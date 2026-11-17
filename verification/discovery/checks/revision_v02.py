#!/usr/bin/env python3
"""Reproducible checks for manuscript v0.2.

Exact rational checks certify the displayed finite tank rule. Floating-point
checks corroborate transforms and continuum theorems; they are not proofs of
coverage, global optimization, or roundoff error. Run from any directory.
"""
from __future__ import annotations
from fractions import Fraction as F
from math import comb, exp, factorial, sqrt
from pathlib import Path
import importlib.util
import json
import sys
import numpy as np
from scipy.optimize import brentq, linprog, minimize_scalar
from scipy.stats import poisson
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'checks'

def exact_game() -> dict:
    rule = {
        3: [((3,5),F(13,14)),((3,7),F(3,140)),((5,7),F(1,20))],
        4: [((4,6),F(23,84)),((4,8),F(7,15)),((4,9),F(22,105)),((5,10),F(1,20))],
        5: [((5,7),F(137,840)),((5,12),F(211,280)),((6,12),F(1,12))],
        6: [((6,9),F(41,350)),((6,10),F(21,100)),((6,11),F(11,40)),((6,12),F(547,1400)),((7,12),F(1,140))],
        7: [((7,12),F(1))],8: [((8,12),F(1))],9: [((9,12),F(1))],
        10:[((10,12),F(1))],11:[((11,12),F(1))],12:[((12,12),F(1))],
    }
    params=range(3,13); checks=0; rows=[]
    P=lambda n,m:F(comb(m-1,2),comb(n,3)) if m<=n else F(0)
    for m,actions in rule.items():
        assert sum(w for _,w in actions)==1 and all(w>=0 for _,w in actions); checks+=1
        assert all(m<=a<=b<=12 for (a,b),_ in actions);checks+=1
    for n in params:
        cov=sum(P(n,m)*w*(a<=n<=b) for m in params for (a,b),w in rule[m])
        risk=sum(P(n,m)*w*(b-a) for m in params for (a,b),w in rule[m])
        assert cov==F(19,20);checks+=1
        assert risk<=F(697,140);checks+=1
        rows.append({'N':n,'coverage':str(cov),'risk':str(risk)})
    assert [r['risk'] for r in rows[-4:]]==['1369/336','1657/480','167/60','167/80'];checks+=1
    lam={n:F(int(n==7)) for n in params}
    gam={n:F(comb(n,3),35) for n in params}
    v={m:-F(comb(m-1,2),35)*(1 if m<=7 else 13-m) for m in params}
    actions=[None]+[(a,b) for a in params for b in range(a,13)]
    for m in params:
        for I in actions:
            length=0 if I is None else I[1]-I[0]
            rhs=sum(P(n,m)*(lam[n]*length-gam[n]*int(I is not None and I[0]<=n<=I[1])) for n in params)
            assert v[m]<=rhs;checks+=1
    assert sum(lam.values())==1 and sum(gam.values())==F(143,7) and sum(v.values())==-F(101,7);checks+=1
    dual=F(19,20)*sum(gam.values())+sum(v.values())
    assert dual==max(F(r['risk']) for r in rows)==F(697,140);checks+=1
    cert={'version':'0.2-support-clipped','parameters':list(params),'coverage':'19/20','optimal_value':str(dual),
          'rule':{str(m):[{'interval':list(I),'probability':str(w)} for I,w in aa] for m,aa in rule.items()},
          'parameter_results':rows,'dual':{'risk_weights':{str(n):str(w) for n,w in lam.items()},
          'coverage_weights':{str(n):str(w) for n,w in gam.items()},'observation_weights':{str(n):str(w) for n,w in v.items()}},
          'exact_rational_constraints_checked':checks}
    (OUT/'finite_game_certificate_clipped.json').write_text(json.dumps(cert,indent=2)+'\n')
    return {'checks':checks,'optimal_value':str(dual),'risks':rows,'arithmetic':'exact rational'}

def joint_checks() -> dict:
    rng=np.random.default_rng(20261112); errors=[]; directions=0
    counts=np.arange(200)
    for t in (.4,2.,5.):
        for u in (.3,3.):
            for x in (.02,.2,.5,1.):
                probs=poisson.pmf(counts,t*x)
                b=np.array([exp(-t*x)*(t*x)**r/factorial(r) for r in range(1,4)])
                q=exp(-t*x)*(-np.expm1(-u*x)); future=-np.expm1(-u*x)
                for _ in range(5):
                    lam=rng.uniform(-.4,.4,3);eta=float(rng.uniform(-.3,.3));c=rng.uniform(-.6,.6,3)
                    base=np.zeros(200);base[1:4]=lam
                    pred=np.zeros(200);pred[1:4]=c
                    direct_M=np.dot(probs,np.exp(base+eta*x*(counts==0)))
                    form_M=1+np.dot(b,np.expm1(lam))+exp(-t*x)*np.expm1(eta*x)
                    direct_Y=np.dot(probs[1:],np.exp(base[1:]))+probs[0]*(1-future+future*exp(eta))
                    form_Y=1+np.dot(b,np.expm1(lam))+q*np.expm1(eta)
                    direct_E=np.dot(probs[1:],np.exp(base[1:]+eta*pred[1:]))+probs[0]*(1-future+future*exp(-eta))
                    form_E=1+np.dot(b,np.expm1(lam+eta*c))+q*np.expm1(-eta)
                    errors.extend([abs(direct_M-form_M),abs(direct_Y-form_Y),abs(direct_E-form_E)]);directions+=3
    assert max(errors)<3e-13
    # Symbolic first-order expansions give the continuous limits at zero.
    x,t,u,eta,l1,l2,c1,c2=sp.symbols('x t u eta l1 l2 c1 c2',real=True)
    b1=sp.exp(-t*x)*t*x;b2=sp.exp(-t*x)*(t*x)**2/2
    q=sp.exp(-t*x)*(1-sp.exp(-u*x))
    ins=[1+b1*(sp.exp(l1)-1)+b2*(sp.exp(l2)-1)+sp.exp(-t*x)*(sp.exp(eta*x)-1),
         1+b1*(sp.exp(l1)-1)+b2*(sp.exp(l2)-1)+q*(sp.exp(eta)-1),
         1+b1*(sp.exp(l1+eta*c1)-1)+b2*(sp.exp(l2+eta*c2)-1)+q*(sp.exp(-eta)-1)]
    targets=[t*(sp.exp(l1)-1)+eta,t*(sp.exp(l1)-1)+u*(sp.exp(eta)-1),t*(sp.exp(l1+eta*c1)-1)+u*(sp.exp(-eta)-1)]
    for A,target in zip(ins,targets):
        assert A.subs(x,0)==1
        assert sp.simplify(sp.diff(A,x).subs(x,0)-target)==0
    return {'enumerated_comparisons':directions,'max_absolute_error':max(errors),'symbolic_zero_limits':3}

def moment_checks() -> dict:
    p=np.array([.5,.3,.2]);r=np.arange(180);t=1.3;u=1.7;ratio=u/t;errors=[];v2s=[];v4s=[]
    for x in p:
        pr=poisson.pmf(r,t*x);pred=np.zeros(180);pred[1:]=(-1.)**(r[1:]+1)*ratio**r[1:]
        q=exp(-t*x)*(1-exp(-u*x));v2=exp(-t*x)*(exp(t*ratio**2*x)-exp(-u*x));v4=exp(-t*x)*(exp(t*ratio**4*x)-exp(-u*x))
        errors.extend([abs(np.dot(pr,pred**2)+q-v2),abs(np.dot(pr,pred**4)+q-v4)])
        v2s.append(v2);v4s.append(v4)
    V=sum(v2s);K4=sum(a-3*b*b for a,b in zip(v4s,v2s));M4=K4+3*V*V
    alt=sum(v4s)+6*sum(v2s[i]*v2s[j] for i in range(3) for j in range(i+1,3))
    assert max(errors)<1e-11 and abs(M4-alt)<1e-11
    V1=exp(-1)*(exp(4)-exp(-2));M41=exp(-1)*(exp(16)-exp(-2));beta=.05
    examples={'chebyshev':sqrt(V1/beta),'fourth_moment':(M41/beta)**.25,'threshold_beta':V1*V1/M41}
    assert examples['chebyshev']<examples['fourth_moment']
    return {'max_enumeration_error':max(errors),'V':V,'K4':K4,'M4':M4,'counterexample':examples}

def geometry_checks() -> dict:
    spec=importlib.util.spec_from_file_location('rg_v02',ROOT/'originals/recovery_geometry/recovery_geometry.py')
    rg=importlib.util.module_from_spec(spec);sys.modules[spec.name]=rg;spec.loader.exec_module(rg)
    rng=np.random.default_rng(20261113);count=0
    for _ in range(60):
        a,b=.04,.8;t=float(rng.uniform(.4,4.))
        mu=rg.Measure(rng.uniform(a,b,6),rng.dirichlet(np.ones(6)))
        env=rg.rate_envelopes(a,b,t,mu.rate(t),mu.integral(lambda x:x*np.exp(-t*x)))
        for q in (-10.,-3.,-1.,-.1,0.,.5,1.,2.,3.,5.):
            # The retained solver may restrict q's advertised range; use the definition directly.
            def hill(m):
                return float(np.exp(-m.integral(np.log))) if q==1 else float(m.integral(lambda x:x**(q-1))**(1/(1-q)))
            lo,val,hi=map(hill,(env.lower,mu,env.upper))
            assert lo<=val+2e-8 and val<=hi+2e-8;count+=1
    # Intermediate sign-change case; use dense sampling plus bounded optimization
    # of every local bracket, not an unsupported global unimodality assumption.
    p=np.array([.5,.3,.2]);a,b,t,q=.05,.6,1.,2.25
    r0=float(np.dot(p,np.exp(-t*p)));r1=float(np.dot(p*p,np.exp(-t*p)));beta=r1/r0;C=1/r0
    E=lambda x:np.exp(t*x)
    moment=lambda x,y:((y-beta)*E(x)+(beta-x)*E(y))/(y-x)
    c=brentq(lambda y:moment(a,y)-C,beta+1e-10,b)
    d=brentq(lambda x:moment(x,b)-C,a,beta-1e-10)
    def two(x):
        y=c if abs(x-a)<1e-12 else b if abs(x-d)<1e-12 else brentq(lambda y:moment(x,y)-C,beta+1e-10,b+1e-10)
        w=(y-beta)/(y-x)
        return r0*(w*E(x)*x**(q-1)+(1-w)*E(y)*y**(q-1))
    def three(z):
        xs=np.array([a,z,b]);w=np.linalg.solve(np.array([np.ones(3),xs,E(xs)]),[1,beta,C])
        assert np.min(w)>-1e-8
        return r0*np.dot(w,E(xs)*xs**(q-1))
    families={}
    for name,f,interval in [('two',two,(a,d)),('three',three,(d,c))]:
        grid=np.linspace(*interval,401);vals=np.array([f(x) for x in grid]);result={}
        for sign,kind in [(1,'min'),(-1,'max')]:
            candidates=[(grid[0],vals[0]),(grid[-1],vals[-1])]
            for j in range(1,len(grid)-1):
                if sign*vals[j]<=min(sign*vals[j-1],sign*vals[j+1]):
                    sol=minimize_scalar(lambda x:sign*f(x),bounds=(grid[j-1],grid[j+1]),method='bounded',options={'xatol':1e-14})
                    candidates.append((sol.x,f(sol.x)))
            z,val=min(candidates,key=lambda e:sign*e[1]);result[kind]={'node':float(z),'value':float(val)}
        families[name]=result
    xs=np.unique(np.r_[np.linspace(a,b,5001),p]);A=np.array([np.ones_like(xs),np.exp(-t*xs),xs*np.exp(-t*xs)])
    rhs=[1,r0,r1];g=xs**(q-1)
    lp1=linprog(g,A_eq=A,b_eq=rhs,bounds=(0,None),method='highs')
    lp2=linprog(-g,A_eq=A,b_eq=rhs,bounds=(0,None),method='highs')
    lo=min(v['min']['value'] for v in families.values());hi=max(v['max']['value'] for v in families.values())
    assert lp1.success and lp2.success and abs(lp1.fun-lo)<2e-6 and abs(-lp2.fun-hi)<2e-6
    return {'hill_order_comparisons':count,'intermediate_families':families,'grid_min':float(lp1.fun),'grid_max':float(-lp2.fun),'grid_absolute_discrepancy':max(abs(lp1.fun-lo),abs(-lp2.fun-hi)),
            'note':'Numerical corroboration only. Analytic zero-counting proves family exhaustiveness.'}

def truncation_and_perspective() -> dict:
    rows=[]
    for p in (np.array([1.]),np.array([.5,.3,.2]),np.ones(10)/10):
        for t in (.3,2.,8.):
            for R in (1,2,5,10):
                logs=poisson.logcdf(R,t*p);sf=poisson.sf(R,t*p)
                direct=1-float(np.prod(poisson.cdf(R,t*p)))
                kernel=-float(np.expm1(np.dot(p,logs/p)))
                union=float(np.dot(p,sf/p))
                assert abs(direct-kernel)<2e-14 and kernel<=union+2e-14
                rows.append(abs(direct-kernel))
    v,z=sp.symbols('v z',positive=True);K=sp.Function('K')
    identity=sp.diff(v*K(1/v),v,2)-sp.Subs(sp.diff(K(z),z,2),z,1/v)/v**3
    assert sp.simplify(identity)==0
    # One-sided exponential quadratic inequality, including large negative values.
    bennett_checks=0
    for c in (.1,1.,3.):
        for lam in (.05,.5,2.):
            for zz in np.linspace(-10,c,101):
                rhs=1+lam*zz+(zz/c)**2*np.expm1(lam*c)-(zz/c)**2*lam*c
                assert exp(lam*zz)<=rhs+1e-11;bennett_checks+=1
    return {'truncation_comparisons':len(rows),'max_absolute_error':max(rows),'perspective_identity':'symbolically verified',
            'one_sided_exponential_inequality_examples':bennett_checks}

def main() -> None:
    report={'scope':'Exact arithmetic for the finite game; other checks explicitly numerical or symbolic.',
            'finite_game':exact_game(),'joint_transforms':joint_checks(),'moments':moment_checks(),
            'recovery':geometry_checks(),'truncation_perspective':truncation_and_perspective()}
    (OUT/'revision_v02.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':
    main()
