"""Independent exact finite-design audits for the capture–recapture module."""
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations, product
from math import comb, inf, prod
from pathlib import Path
import json
import sympy as sp
from partial_population import capture_recapture as cr

counts=Counter()
def check(condition, family):
    if not condition:
        raise AssertionError(f"failed: {family}")
    counts[family]+=1

def contains(interval,N):
    L,U=interval
    return L<=N and (U is None or N<=U)

def C(n,k):
    return comb(n,k) if 0<=k<=n else 0

def direct_pmf(N,a,b):
    law=Counter(len(set(x)&set(range(a))) for x in combinations(range(N),b))
    return tuple(F(law[j],comb(N,b)) for j in range(min(a,b)+1))

# Explicit labelled subsets are independent of the implementation's formula.
for N in range(1,9):
  for a in range(N+1):
    for b in range(N+1):
      check(cr.hypergeom_pmf(N,a,b)==direct_pmf(N,a,b),'subset_kernel')

for N in range(1,21):
  for a in range(1,N+1):
    for b in range(1,N+1):
      p=cr.hypergeom_pmf(N,a,b)
      check(sum(p)==1,'normalization')
      for r in range(min(a,b)+1):
        lhs=sum(v*prod(range(j-r+1,j+1)) for j,v in enumerate(p))
        rhs=F(prod(range(a-r+1,a+1))*prod(range(b-r+1,b+1)),prod(range(N-r+1,N+1)))
        check(lhs==rhs,'factorial_moments')
      mean,var,bias=cr.chapman_moments(N,a,b)
      check(bias==-(N-a-b)*p[0],'chapman_bias')
      check(var>=0 and mean-N==bias,'chapman_variance')
      for k in range(1,5):
        lhs=sum(v*cr.shifted_reciprocal(a,b,j,k) for j,v in enumerate(p))
        boundary=F(prod(range(b+1,b+k+1))*sum(C(a+k,h)*C(N-a,b+k-h) for h in range(k)),C(N,b))
        check(lhs==prod(range(N+1,N+k+1))-boundary,'reciprocal_shifts')
      if a==b==1:
        check(var==F(4,N)-F(4,N*N),'variance_closed_example')

for a,b in product(range(1,6),repeat=2):
    m=min(a,b)
    matrix=sp.Matrix([[sp.Rational(x.numerator,x.denominator) for x in cr.hypergeom_pmf(N,a,b)] for N in range(max(a,b),max(a,b)+m+1)])
    check(matrix.det()!=0,'range_rank')
    for j in range(1,m+1):
        nmin=a+b-j
        cap=a*b//j+2
        likelihood={N:cr.hypergeom_pmf(N,a,b)[j] for N in range(nmin,cap+1)}
        maximum=max(likelihood.values())
        actual={N for N,p in likelihood.items() if p==maximum}
        expected={a*b//j}
        if a*b%j==0 and a*b//j-1>=nmin:
            expected.add(a*b//j-1)
        check(actual==expected,'integer_mle')

# The implementation searches unbounded N; this audit evaluates the resulting
# exact endpoints on finite truth grids and checks their sharpness witnesses.
for a,b in product(range(1,8),repeat=2):
  for alpha in (F(1,20),F(1,4)):
    intervals=[cr.fixed_interval(a,b,j,alpha) for j in range(min(a,b)+1)]
    for N in range(max(a,b),41):
        p=cr.hypergeom_pmf(N,a,b)
        coverage=sum(v for j,v in enumerate(p) if contains(intervals[j],N))
        check(coverage>=1-alpha,'fixed_coverage')
    for j in range(min(a,b)+1):
        L=cr.fixed_lower_limit(a,b,j,alpha)
        U=cr.fixed_upper_limit(a,b,j,alpha)
        check(sum(cr.hypergeom_pmf(L,a,b)[:j+1])>alpha,'lower_sharpness_witness')
        if L>max(a,b):
            check(sum(cr.hypergeom_pmf(L-1,a,b)[:j+1])<=alpha,'lower_crossing')
        if U is None:
            check(j==0,'infinite_zero_overlap_limit')
        else:
            check(sum(cr.hypergeom_pmf(U,a,b)[j:])>alpha,'upper_sharpness_witness')
            check(sum(cr.hypergeom_pmf(U+1,a,b)[j:])<=alpha,'upper_crossing')
        if j:
            check(L<=cr.fixed_lower_limit(a,b,j-1,alpha),'limit_monotonicity')
            old=cr.fixed_upper_limit(a,b,j-1,alpha)
            check(old is None or U<=old,'limit_monotonicity')
check(cr.fixed_upper_limit(1,1,1,F(1,2))==1,'strict_boundary')
check(cr.fixed_lower_limit(1,1,0,F(1,2))==3,'strict_boundary')

for a,b,c,d in product(range(1,5),repeat=4):
    if min(c,d)>min(a,b):
        try:
            cr.signed_transfer(a,b,c,d)
        except ValueError:
            check(True,'transfer_nonexistence_rejection')
        else:
            check(False,'transfer_nonexistence_rejection')
        continue
    S=cr.signed_transfer(a,b,c,d)
    check(all(sum(row)==1 for row in S),'transfer_row_mass')
    for N in range(max(a,b,c,d),max(a,b,c,d)+8):
        p=cr.hypergeom_pmf(N,a,b)
        out=tuple(sum(p[j]*S[j][y] for j in range(len(p))) for y in range(min(c,d)+1))
        check(out==cr.hypergeom_pmf(N,c,d),'signed_transfer_law')
    if c<=a and d<=b:
        check(all(x>=0 for row in S for x in row),'positive_deletion')
        for j in range(min(a,b)+1):
            first=cr.hypergeom_pmf(a,j,c)
            composed=[sum(v*F(C(h,y)*C(b-h,d-y),C(b,d)) for h,v in enumerate(first)) for y in range(min(c,d)+1)]
            check(tuple(composed)==S[j],'deletion_composition')
S=cr.signed_transfer(2,2,1,3)
check(S==((F(1),F(0)),(F(1,4),F(3,4)),(F(-1,2),F(3,2))),'signed_counterexample')
check(max(sum(abs(x) for x in row) for row in S)==2,'signed_norm')

for N in range(1,21):
  for a in range(1,N+1):
    for r in range(1,a+1):
      p=cr.inverse_hypergeom_pmf(N,a,r)
      check(sum(p.values())==1,'inverse_normalization')
      check(sum(v*t for t,v in p.items())==F(r*(N+1),a+1),'inverse_mean')
      check(sum(v*(cr.sequential_estimate(a,r,t)-N)**2 for t,v in p.items())==cr.sequential_variance(N,a,r),'inverse_variance')
      if N<=8:
        law=Counter(sorted(s)[r-1]+1 for s in combinations(range(N),a))
        check({t:F(v,comb(N,a)) for t,v in law.items()}==p,'order_statistic_enumeration')
      if r==a:
        check(all(v==F(C(t-1,a-1),C(N,a)) for t,v in p.items()),'tank_equivalence')
for a in range(1,6):
  for r in range(1,a+1):
    for N in range(a,21):
      p=cr.inverse_hypergeom_pmf(N,a,r)
      coverage=sum(v for t,v in p.items() if contains(cr.sequential_interval(a,r,t),N))
      check(coverage>=F(19,20),'inverse_coverage')

for N in range(1,16):
  for a in range(1,N+1):
    for b in range(1,8):
      p=cr.binomial_pmf(N,a,b)
      check(sum(p)==1,'binomial_normalization')
      check(sum(v*cr.bailey(a,b,j) for j,v in enumerate(p))==N*(1-(1-F(a,N))**(b+1)),'bailey_bias')
      if b<=4 and N<=8:
        intervals=[cr.binomial_interval(a,b,j) for j in range(b+1)]
        check(sum(v for j,v in enumerate(p) if contains(intervals[j],N))>=F(19,20),'binomial_coverage')
z,p=sp.symbols('z p',positive=True)
for r in range(1,6):
    pgf=(p*z/(1-(1-p)*z))**r
    mean=sp.diff(pgf,z).subs(z,1)
    variance=sp.diff(pgf,z,2).subs(z,1)+mean-mean**2
    check(sp.simplify(mean-r/p)==0,'inverse_binomial_symbolic_mean')
    check(sp.simplify(variance-r*(1-p)/p**2)==0,'inverse_binomial_symbolic_variance')
# Infinite-support coverage: sum until the entire remaining upper tail is
# rejected by the true parameter's survival test, then close that tail exactly.
for N in range(1,7):
  for a in range(1,N+1):
    for r in (1,2):
      prob=F(a,N)**r
      previous_cdf=F(0)
      t=r
      coverage=F(0)
      while 1-previous_cdf>F(1,40):
        interval=cr.inverse_binomial_interval(a,r,t)
        if contains(interval,N):
            coverage+=prob
        previous_cdf+=prob
        prob*=F(t,t-r+1)*(1-F(a,N))
        t+=1
      check(coverage>=F(19,20),'inverse_binomial_exact_tail_coverage')

# Repeated-capture laws from actual labelled subsets, not a count-state kernel.
for N in range(2,6):
  for initial in range(1,N+1):
    sizes=(1,min(2,N),min(2,N))
    for weighted in (False,True):
      law=defaultdict(F)
      def visit(marked,counts_so_far,js,prob,t):
        if t==len(sizes):
          estimate=cr.repeated_estimate(counts_so_far,sizes,js,weight_marked=weighted)
          law[estimate]+=prob
          return
        b=sizes[t]
        for sample in combinations(range(N),b):
          sample=set(sample)
          visit(marked|sample,counts_so_far+(len(marked),),js+(len(marked&sample),),prob/F(comb(N,b)),t+1)
      visit(set(range(initial)),(),(),F(1),0)
      check(dict(law)==cr.repeated_distribution(N,initial,sizes,weight_marked=weighted),'repeated_subset_law')

# Adaptive but predictable betting; every terminal history is enumerated.
for N in range(2,9):
  for initial in range(1,N+1):
    expected=defaultdict(F)
    crossing=F(0)
    def walk(a,js,factors,prob,t):
      global crossing
      if t==4:
        path=cr.capture_evalue_path(N,initial,(2,)*4,js,factors)
        for s,value in enumerate(path): expected[s]+=prob*value
        if max(path)>=4: crossing+=prob
        return
      factor=F(2) if not js or js[-1]==0 else F(1,2)
      for j,v in enumerate(cr.hypergeom_pmf(N,a,2)):
        if v: walk(a+2-j,js+(j,),factors+(factor,),prob*v,t+1)
    walk(initial,(),(),F(1),0)
    check(all(x==1 for x in expected.values()),'repeated_martingale_expectation')
    check(crossing<=F(1,4),'repeated_anytime_crossing')
for N in range(1,10):
  for a in range(1,N+1):
    expected=[F(0)]*N
    crossing=F(0)
    for positions in combinations(range(N),a):
      draws=tuple(int(s in positions) for s in range(N))
      factors=tuple(F(3) if sum(draws[:s])*2<=s else F(1,3) for s in range(N))
      path=cr.draw_evalue_path(N,a,draws,factors)
      prob=F(1,comb(N,a))
      for s,value in enumerate(path): expected[s]+=prob*value
      if max(path)>=4: crossing+=prob
    check(all(x==1 for x in expected),'draw_martingale_expectation')
    check(crossing<=F(1,4),'draw_anytime_crossing')

populations=[((F(1,3),F(2,5)),),((F(1,3),F(2,5)),(F(2,3),F(1,4))),((F(0),F(0)),(F(1),F(1)),(F(1,2),F(1,3)))]
for population in populations:
  for weights in ((F(1),F(2),F(3),F(1,2)),(F(3,2),F(1,3),F(5,2),F(2))):
    direct=F(0)
    for bits in product((0,1),repeat=2*len(population)):
      probability=prod(q[k] if bits[2*i+k] else 1-q[k] for i,q in enumerate(population) for k in range(2))
      score=prod(weights[bits[2*i]+2*bits[2*i+1]] for i in range(len(population)))
      direct+=probability*score
    check(direct==cr.history_transform(population,weights),'history_transform_enumeration')
# Bernoulli two-occasion conditioning removes the homogeneous nuisance rates.
for N in range(1,6):
  raw=defaultdict(F)
  q1,q2=F(1,3),F(2,5)
  for bits in product((0,1),repeat=2*N):
    a=sum(bits[0::2]); b=sum(bits[1::2]); j=sum(bits[2*i]*bits[2*i+1] for i in range(N))
    probability=prod((q1 if bits[2*i] else 1-q1)*(q2 if bits[2*i+1] else 1-q2) for i in range(N))
    raw[a,b,j]+=probability
  for a,b in product(range(N+1),repeat=2):
    mass=sum(raw[a,b,j] for j in range(min(a,b)+1))
    check(tuple(raw[a,b,j]/mass for j in range(min(a,b)+1))==cr.hypergeom_pmf(N,a,b),'bernoulli_conditioning')
for qs in ((F(1,2),F(1,2)),(F(1,4),F(2,3),F(1)),(F(0),F(1,3),F(3,4))):
  for K in range(2,9):
    m=[sum(comb(K,j)*q**j*(1-q)**(K-j) for q in qs) for j in range(3)]
    bound=cr.binomial_moment_lower_bound(K,m[1],m[2])
    check(bound<=m[0],'binomial_moment_bound')
    if len(set(qs))==1: check(bound==m[0],'moment_bound_equality')
for rho in ((F(1,3),F(2,3)),(F(1,4),F(1,2),F(1))):
    law=[]
    for bits in product((0,1),repeat=len(rho)):
        prob=prod(p if bit else 1-p for p,bit in zip(rho,bits))
        estimate=cr.horvitz_thompson(p for p,bit in zip(rho,bits) if bit)
        law.append((prob,estimate))
    check(sum(p*x for p,x in law)==len(rho),'HT_mean')
    check(sum(p*(x-len(rho))**2 for p,x in law)==sum((1-p)/p for p in rho),'HT_variance')

# Dependent fixed-size inclusions: HT is constant and covariance cancels.
for N in range(2,8):
  for b in range(1,N+1):
    rho=F(b,N); pair=F(b*(b-1),N*(N-1))
    variance=N*(1-rho)/rho+N*(N-1)*(pair-rho*rho)/(rho*rho)
    check(variance==0 and cr.horvitz_thompson((rho,)*b)==N,'HT_dependent_design')
check(cr.chao_incidence(2,(4,2))==8,'chao_point_convention')
check(cr.chao_incidence(3,(2,0,1))==inf,'chao_point_convention')
check(cr.chao_incidence(3,(0,0,2))==2,'chao_point_convention')

# Error handling protects design distinctions and exact threshold arithmetic.
invalid=[lambda:cr.hypergeom_pmf(2,3,1),lambda:cr.chapman(1,2,2),lambda:cr.fixed_interval(2,2,1,0.05),lambda:cr.sequential_estimate(1,2,2),lambda:cr.bailey(1,2,3),lambda:cr.repeated_estimate((1,1),(2,2),(0,0)),lambda:cr.capture_evalue_path(3,1,(1,),(0,),(0,)),lambda:cr.history_transform(((F(1,2),),),(1,)),lambda:cr.horvitz_thompson((0,))]
for call in invalid:
    try: call()
    except ValueError: check(True,'input_rejection')
    else: check(False,'input_rejection')
check(cr.lincoln_petersen(2,2,0)==inf,'zero_ratio_convention')
check(cr.bailey(1,3,3)==1,'replacement_repeats')
check(cr.inverse_binomial_estimate(1,3,3)==1,'replacement_repeats')
check(cr.capture_evalue_path(3,2,(2,),(0,),(2,))==(inf,),'impossible_candidate')
check(cr.draw_evalue_path(2,1,(1,1),(2,2))[-1]==inf,'impossible_candidate')
report={'arithmetic':'exact rational; symbolic probability-generating derivatives','checks':dict(sorted(counts.items())),'total_checks':sum(counts.values()),'scope':'Finite examples corroborate the analytic proofs; counts are not measures of proof assurance.'}
Path(__file__).with_name('results.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
