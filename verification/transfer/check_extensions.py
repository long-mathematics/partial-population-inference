#!/usr/bin/env python3
"""Independent finite corroboration of Certified Transfer and Generation.

Exact rational algebra checks finite kernels and optimal count couplings.
Enumerations use floating point only where exponential moments are needed.
Grid linear programs check grid duality, NOT continuum certification.
Run: python check_extensions.py
Dependencies: Python >=3.10, numpy, scipy, sympy.
"""
from __future__ import annotations
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
from math import comb, exp, factorial, log, prod
from pathlib import Path
import json
import numpy as np
import sympy as sp
from scipy.optimize import linprog
from scipy.stats import poisson

CHECKS: dict[str, int] = {}
DETAILS: dict[str, object] = {}

def check(name: str, value: bool) -> None:
    if not value:
        raise AssertionError(name)
    CHECKS[name] = CHECKS.get(name, 0) + 1

def choose(n: int, k: int) -> int:
    return comb(n, k) if 0 <= k <= n else 0

def hyper(L: int, n: int) -> sp.Matrix:
    return sp.Matrix(L+1, n+1, lambda A,x:
        sp.Rational(choose(A,x)*choose(L-A,n-x), choose(L,n)))

def binsize(L: int, n: int) -> sp.Matrix:
    return sp.Matrix(L+1,n+1,lambda A,x:
        sp.Rational(choose(n,x))*sp.Rational(A,L)**x *
        (1-sp.Rational(A,L))**(n-x))

def reduce_counts(n: int, m: int) -> sp.Matrix:
    return sp.Matrix(n+1,m+1,lambda x,y:
        sp.Rational(choose(x,y)*choose(n-x,m-y),choose(n,m)))

def verify_kernels() -> None:
    for L in range(2,8):
        for n in range(L+1):
            H=hyper(L,n)
            check('hypergeometric_rank',H.rank()==n+1)
            for j in range(n+1):
                coeff=sp.Matrix([choose(x,j) for x in range(n+1)])
                expected=sp.Matrix([sp.Rational(choose(n,j),choose(L,j))*choose(A,j)
                                    for A in range(L+1)])
                check('hypergeometric_factorial',H*coeff==expected)
            for m in range(n+1):
                K=reduce_counts(n,m)
                check('hypergeometric_reduction',H*K==hyper(L,m))
                check('binomial_size_reduction',binsize(L,n)*K==binsize(L,m))
                check('reduction_stochastic',all(z>=0 for z in K) and
                      K*sp.ones(m+1,1)==sp.ones(n+1,1))
            for m in range(L+2):
                rank=H.row_join(binsize(L,m)).rank()
                check('without_to_with_range_classification',
                      (rank==n+1)==(m<=n or n==L))
        C=binsize(L,2); H=hyper(L,2)
        S=sp.Matrix([[1,0,0],[-sp.Rational(1,2*(L-1)),sp.Rational(L,L-1),
                              -sp.Rational(1,2*(L-1))],[0,0,1]])
        check('cross_design_signed_factor',C*S==H)
        check('cross_design_signed_norm',max(sum(abs(S[i,j]) for j in range(3))
                        for i in range(3))==sp.Rational(L+1,L-1))
        # Positive transfer LP is infeasible, as the analytic endpoint argument says.
        # Variables S[x,y], x,y=0,1,2.
        rows=[]; rhs=[]
        for A in range(L+1):
            for y in range(3):
                row=np.zeros(9)
                for x in range(3): row[3*x+y]=float(C[A,x])
                rows.append(row); rhs.append(float(H[A,y]))
        for x in range(3):
            row=np.zeros(9); row[3*x:3*x+3]=1
            rows.append(row); rhs.append(1.)
        lp=linprog(np.zeros(9),A_eq=rows,b_eq=rhs,bounds=(0,None),method='highs')
        check('cross_design_positive_obstruction_LP',lp.status==2)


def cdf(v: list[F]) -> list[F]:
    out=[]; s=F(0)
    for x in v: s+=x; out.append(s)
    return out

def tail_gap(p: list[F],q: list[F]) -> F:
    size=max(len(p),len(q)); pp=p+[F(0)]*(size-len(p)); qq=q+[F(0)]*(size-len(q))
    return max([F(0)]+[sum(qq[k:])-sum(pp[k:]) for k in range(1,size)])

def inverse_cdf(c: list[F],u: F) -> int:
    return next(i for i,z in enumerate(c) if u<z)

def quantile_coupling(p: list[F],q: list[F]) -> list[list[F]]:
    d=tail_gap(p,q); cp=cdf(p); cq=cdf(q)
    points={F(0),F(1),d,*cp}
    # Breakpoints for V=(U-d) mod 1.
    for z in cq:
        w=z+d
        points.add(w if w<=1 else w-1)
    points=sorted(z for z in points if 0<=z<=1)
    out=[[F(0) for _ in q] for _ in p]
    for lo,hi in zip(points,points[1:]):
        if hi==lo: continue
        u=(lo+hi)/2
        v=u-d if u>d else u+1-d
        i=inverse_cdf(cp,u); j=inverse_cdf(cq,v)
        out[i][j]+=hi-lo
    return out

def verify_count_couplings() -> None:
    rng=np.random.default_rng(20261112)
    for _ in range(120):
        n=int(rng.integers(2,8)); m=int(rng.integers(2,8))
        pi=[int(z) for z in rng.integers(0,7,n)]; qi=[int(z) for z in rng.integers(0,7,m)]
        if sum(pi)==0: pi[0]=1
        if sum(qi)==0: qi[0]=1
        p=[F(z,sum(pi)) for z in pi]; q=[F(z,sum(qi)) for z in qi]
        C=quantile_coupling(p,q); d=tail_gap(p,q)
        check('count_coupling_rows',[sum(row) for row in C]==p)
        check('count_coupling_columns',[sum(C[i][j] for i in range(n)) for j in range(m)]==q)
        failure=sum(C[i][j] for i in range(n) for j in range(m) if j>i)
        check('count_coupling_exact_optimum',failure==d)
        # Independent numerical transport LP for the same finite instance.
        costs=np.array([float(j>i) for i in range(n) for j in range(m)])
        A=[]; b=[]
        for i in range(n):
            row=np.zeros(n*m); row[i*m:(i+1)*m]=1; A.append(row); b.append(float(p[i]))
        for j in range(m):
            row=np.zeros(n*m); row[j::m]=1; A.append(row); b.append(float(q[j]))
        lp=linprog(costs,A_eq=A,b_eq=b,bounds=(0,None),method='highs')
        check('count_coupling_transport_LP',lp.success and abs(lp.fun-float(d))<1e-10)
    DETAILS['random_count_laws']=120


def verify_signed_repairs() -> None:
    for a in [F(0),F(1,5),F(2),F(7,2)]:
        row=[1+a,-a]; gamma=sum(abs(x) for x in row)
        K=[F(1),F(0)]
        check('rectification_row_cost',sum(abs(x-y) for x,y in zip(row,K))/2==(gamma-1)/2)
        for e in product([F(0),F(1,2),F(1)],repeat=2):
            v=sum(s*z for s,z in zip(row,e))
            check('signed_e_nonnegative_repair',v+a>=0)
    # Poisson transform of the algebraic thinning polynomial, symbolic.
    x,t,a,z=sp.symbols('x t a z',positive=True)
    check('poisson_pgf_factor',sp.simplify(sp.exp(t*x*((1-a+a*z)-1))-
                                           sp.exp(t*x*a*(z-1)))==0)


def subsets(seq: tuple) -> list[list[tuple]]:
    return [[tuple(seq[i] for i in I) for I in combinations(range(len(seq)),k)]
             for k in range(len(seq)+1)]

def cover(S: tuple[frozenset,...]) -> frozenset:
    return frozenset().union(*S)

def novelty(S: tuple[frozenset,...],T: tuple[frozenset,...]) -> int:
    return len(cover(T)-cover(S))

def verify_novelty() -> None:
    marks=(frozenset(),frozenset({0}),frozenset({1}),frozenset({0,1}))
    for S in product(marks,repeat=2):
        SS=subsets(S); DS=len(cover(S))
        for T in product(marks,repeat=2):
            TT=subsets(T); Y=novelty(S,T)
            for k in range(3):
                for SI in SS[k]:
                    lost=DS-len(cover(SI))
                    for ell in range(3):
                        for TJ in TT[ell]:
                            small=novelty(SI,TJ)
                            check('pathwise_novelty_sandwich',
                                  small-min(lost,2*ell)<=Y<=small+2*(2-ell))
    # Saturation examples for each coefficient.
    check('novelty_upper_coefficient_attained',novelty((),(frozenset({0,1}),))==2)
    check('novelty_lower_coefficient_attained',
          novelty((),(frozenset({0,1}),))-novelty((frozenset({0,1}),),(frozenset({0,1}),))==2)


def fingerprint(seq: tuple[int,...],d: int) -> np.ndarray:
    cc=Counter(seq)
    return np.array([sum(v==r for v in cc.values()) for r in range(1,d+1)],dtype=float)

def kernel_grid(xs: np.ndarray,t: float,u: float,lam: np.ndarray,eta: float) -> np.ndarray:
    out=np.empty_like(xs,dtype=float)
    for j,x in enumerate(xs):
        if x==0:
            out[j]=t*np.expm1(lam[0])+u*np.expm1(eta)
        else:
            bs=np.array([exp(-t*x)*(t*x)**r/factorial(r) for r in range(1,len(lam)+1)])
            q=exp(-t*x)*(-np.expm1(-u*x))
            out[j]=np.log1p(np.dot(bs,np.expm1(lam))+q*np.expm1(eta))/x
    return out

def Hcoefficient(S: tuple[int,...],m: int,t: float,u: float,lam: np.ndarray,eta: float) -> float:
    SS=subsets(S); n=len(S); D=len(set(S)); total=0.
    for k in range(n+1):
        for ell in range(m+1):
            block=0.
            for SI in SS[k]:
                lost=D-len(set(SI))
                exponent=float(np.dot(lam,fingerprint(SI,len(lam))))
                exponent-=max(eta,0)*(m-ell)+max(-eta,0)*min(lost,ell)
                block+=exp(exponent)
            total+=poisson.pmf(k,t)*poisson.pmf(ell,u)*block/comb(n,k)
    return float(total)

def sequence_law(p: np.ndarray,n: int):
    for S in product(range(len(p)),repeat=n):
        yield S,prod(float(p[i]) for i in S)

def verify_joint_transfer() -> None:
    trials=[]; largest=0.; base_residual=0.
    for p in [np.array([1.]),np.array([.3,.7]),np.array([.5,.3,.2])]:
        for n,m,t,u in [(2,2,1.4,1.1),(3,2,2.2,1.5),(4,1,3.0,.8)]:
            for lam,eta in [(np.array([0.,0.]),0.),(np.array([.7,-.4]),.6),
                            (np.array([-.5,.8]),-.9),(np.array([1.1,-.2]),-1.3)]:
                K=float(np.dot(p,kernel_grid(p,t,u,lam,eta)))
                expected=0.; expected_base=0.
                for S,prS in sequence_law(p,n):
                    H=Hcoefficient(S,m,t,u,lam,eta)
                    baseline=float(poisson.pmf(n,t)*poisson.pmf(m,u)*
                                   exp(np.dot(lam,fingerprint(S,len(lam)))))
                    check('joint_H_dominates_full_count_term',H+1e-13>=baseline)
                    for T,prT in sequence_law(p,m):
                        Y=len(set(T)-set(S))
                        expected+=prS*prT*H*exp(eta*Y-K)
                        expected_base+=prS*prT*baseline*exp(eta*Y-K)
                largest=max(largest,expected)
                check('joint_fixed_count_evalue',expected<=1+1e-11)
                check('joint_full_dominates_single_term',expected+1e-12>=expected_base)
                if eta==0 and np.max(abs(lam))==0:
                    exact=poisson.cdf(n,t)*poisson.cdf(m,u)
                    base_residual=max(base_residual,abs(expected-exact))
                    check('constant_joint_evalue_retained_mass',abs(expected-exact)<2e-12)
                trials.append({'population':p.tolist(),'n':n,'m':m,'t':t,'u':u,
                               'lambda':lam.tolist(),'eta':eta,'expectation':expected})
    DETAILS['joint_evalue_enumerations']={'cases':len(trials),'largest_expectation':largest,
            'constant_case_max_error':base_residual,'cases_detail':trials,
            'interpretation':'Finite-population exact-law enumeration with floating-point exponentials; not a universal coverage proof.'}


def verify_grid_dual() -> None:
    # Construct actual transferred constraints; additionally check manufactured
    # feasible tight constraints so primal/dual signs are tested nontrivially.
    S=(0,0,0,0,1,1); n=len(S); m=3; t=4.; u=2.; alpha=.1
    params=[(np.array([l,0.]),e) for l in [-2.,0.,2.] for e in [-1.,1.]]
    xs=np.unique(np.r_[np.linspace(0,1,161),.3,.7]); J=len(params)
    K=np.array([kernel_grid(xs,t,u,l,e) for l,e in params]); etas=np.array([e for _,e in params])
    HH=np.array([Hcoefficient(S,m,t,u,l,e) for l,e in params]); c=np.log(HH)-log(J/alpha)
    cases=[]
    for tag,cc in [('transferred',c),('manufactured_feasible',.4*K[:,48]+.6*K[:,100]-etas*1.5-.03)]:
        # Maximize sign*y; negative minimum reconstructs L.
        for sign in [1.,-1.]:
            ngrid=len(xs)
            obj=np.r_[np.zeros(ngrid),-sign]
            A=np.column_stack([-K,etas]); b=-cc
            eq=np.r_[np.ones(ngrid),0.][None,:]
            primal=linprog(obj,A_ub=A,b_ub=b,A_eq=eq,b_eq=[1.],
                           bounds=[(0,None)]*ngrid+[(0,m)],method='highs')
            # variables a>=0, z unrestricted, h>=0
            dualobj=np.r_[-cc,1.,m]
            Adual=np.column_stack([K.T,-np.ones(ngrid),np.zeros(ngrid)])
            last=np.r_[-etas,0.,-1.][None,:]
            dual=linprog(dualobj,A_ub=np.vstack([Adual,last]),b_ub=np.r_[np.zeros(ngrid),-sign],
                         bounds=[(0,None)]*J+[(None,None),(0,None)],method='highs')
            check('finite_grid_dual_feasible',primal.success and dual.success)
            gap=abs(-primal.fun-dual.fun)
            check('finite_grid_duality_value',gap<1e-8)
            cases.append({'case':tag,'target_sign':sign,'primal_value':float(-primal.fun),
                          'dual_value':float(dual.fun),'gap':float(gap)})
    DETAILS['grid_dual']={'nodes':len(xs),'cases':cases,
        'warning':'Finite grid corroboration only. A grid maximum is not a certified continuum upper bound.'}


def verify_batch_extraction() -> None:
    labels=range(3); B=2
    marks=[frozenset(C) for k in range(B+1) for C in combinations(labels,k)]
    rng=np.random.default_rng(101)
    for _ in range(40):
        w=[int(z) for z in rng.integers(1,8,len(marks))]; w=[F(z,sum(w)) for z in w]
        q=[sum(v*(i in C) for C,v in zip(marks,w)) for i in labels]
        r=[sum(v*F(int(i in C),B) for C,v in zip(marks,w)) for i in labels]
        null=sum(v*(1-F(len(C),B)) for C,v in zip(marks,w))
        check('batch_proportional_rate',r==[v/B for v in q])
        check('batch_normalization',sum(r)+null==1)
        check('batch_subprobability_mass',sum(r)<=1)
    # Equal marginal rates do not determine the full feature law.
    t=1.; qq=1-exp(-t/2)
    check('batch_prediction_joint_obstruction',2*qq*(1-qq)<4*qq*(1-qq))


def main() -> None:
    verify_kernels(); verify_count_couplings(); verify_signed_repairs()
    verify_novelty(); verify_joint_transfer(); verify_grid_dual(); verify_batch_extraction()
    report={'all_passed':True,'total_checks':sum(CHECKS.values()),'groups':CHECKS,
            'details':DETAILS,'scope':'Exact finite algebra or explicitly marked floating-point enumeration. Proofs and novelty are separate.'}
    path=Path(__file__).with_name('extensions_results.json')
    path.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='details'},indent=2))

if __name__=='__main__':
    main()
