#!/usr/bin/env python3
"""Independent checks of the retained rational game certificate and added lemmas.
Run from any directory. Requires numpy, scipy, sympy for non-certificate checks.
Exact rational verification is distinct from floating-point and symbolic tests.
"""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import importlib.util
import json
import numpy as np
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
cert_path=ROOT/'originals/statistical_calibration/finite_game_certificate.json'
c=json.loads(cert_path.read_text())
N=c['parameters']; obs=c['observations']; acts=c['actions']
k=min(N); alpha=F(c['alpha']); R=F(c['optimal_value'])
D=[[F(x) for x in r] for r in c['decision_probabilities']]
y=[F(x) for x in c['inequality_dual']]; z=[F(x) for x in c['equality_dual']]
P=[[F(comb(m-1,k-1),comb(n,k)) if m<=n else F(0) for m in obs] for n in N]
cost=[F(0) if a is None else F(a[1]-a[0]) for a in acts]
inc=[[int(a is not None and a[0]<=n<=a[1]) for a in acts] for n in N]
checks=0
for row in D:
    assert all(x>=0 for x in row) and sum(row)==1; checks+=1
covs=[]; risks=[]
for j,n in enumerate(N):
    cv=sum(P[j][o]*D[o][a]*inc[j][a] for o in range(len(obs)) for a in range(len(acts)))
    ri=sum(P[j][o]*D[o][a]*cost[a] for o in range(len(obs)) for a in range(len(acts)))
    assert cv>=1-alpha and ri<=R; checks+=2
    covs.append(str(cv)); risks.append(str(ri))
assert all(v<=0 for v in y) and -sum(y[::2])<=1; checks+=2
for o in range(len(obs)):
    for a in range(len(acts)):
        val=z[o]+sum(P[j][o]*(cost[a]*y[2*j]-inc[j][a]*y[2*j+1]) for j in range(len(N)))
        assert val<=0; checks+=1
assert -(1-alpha)*sum(y[1::2])+sum(z)==R==F(697,140); checks+=1
x,t,q,s=sp.symbols('x t q s', positive=True)
def Lop(g): return sp.diff(g,x,3)+2*t*sp.diff(g,x,2)+t*t*sp.diff(g,x)
assert sp.simplify(Lop(sp.exp(-s*x))+s*(t-s)**2*sp.exp(-s*x))==0
assert sp.simplify(Lop(x**(q-1))-(q-1)*x**(q-4)*((t*x+q-2)**2+2-q))==0
assert sp.simplify(Lop(-sp.log(x))+((t*x-1)**2+1)/x**3)==0
spec=importlib.util.spec_from_file_location('rg',ROOT/'originals/recovery_geometry/recovery_geometry.py')
import sys
rg=importlib.util.module_from_spec(spec);sys.modules[spec.name]=rg;spec.loader.exec_module(rg)
rng=np.random.default_rng(20261112); extra=0
for _ in range(50):
    a,b=.03,.8; tt=float(rng.uniform(.3,5))
    mu=rg.Measure(rng.uniform(a,b,5),rng.dirichlet(np.ones(5)))
    e=rg.rate_envelopes(a,b,tt,mu.rate(tt),mu.integral(lambda x:x*np.exp(-tt*x)))
    for qq in (3.,3.1,4.,7.,10.):
        lo,v,hi=[m.hill(qq) for m in (e.lower,mu,e.upper)]
        assert lo<=v+2e-8 and v<=hi+2e-8;extra+=1
p=np.array([.5,.3,.2]); mu=rg.Measure(p,p); tt=.3
E=rg.rate_envelopes(.05,.6,tt,mu.rate(tt),mu.integral(lambda x:x*np.exp(-tt*x)))
reverse=[m.hill(2.25) for m in (E.lower,mu,E.upper)]
assert reverse[0]>reverse[1]>reverse[2]
report={'exact_certificate_checks':checks,'optimal_value':str(R),
        'minimum_exact_coverage':str(min(map(F,covs))),'parameter_risks':risks,
        'symbolic_operator_identities':3,'additional_hill_checks':extra,
        'hill_q_2_25_reversed_order':reverse,
        'interpretation':'Rational verification, symbolic identities, and floating-point examples are separately identified.'}
(ROOT/'checks/independent_audit.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
