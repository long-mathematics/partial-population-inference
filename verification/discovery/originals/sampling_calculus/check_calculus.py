from itertools import combinations
from fractions import Fraction as F
from math import comb

units=[{0,1},{0},{1,2},{2},{0,2}]
w=[F(1),F(2),F(3)]
N=len(units)

def coverage(ids):
    hit=set().union(*(units[i] for i in ids)) if ids else set()
    return sum((w[a] for a in hit), F(0))

def phi(ids,r):
    return sum((w[a] for a in range(3) if sum(a in units[i] for i in ids)==r),F(0))

def avg(n,fn):
    return sum((fn(ids) for ids in combinations(range(N),n)),F(0))/comb(N,n)

Fs=[avg(n,coverage) for n in range(N+1)]
checks=0
for n in range(N):
    for r in range(1,N-n+1):
        delta=sum((-1)**(r-j)*comb(r,j)*Fs[n+j] for j in range(r+1))
        rhs=avg(n+r,lambda ids:phi(ids,r))/comb(n+r,r)
        assert (-1)**(r-1)*delta==rhs, (n,r,delta,rhs)
        checks+=1
for n in range(1,N+1):
  for ids in combinations(range(N),n):
    counts=[sum(a in units[i] for i in ids) for a in range(3)]
    for m in range(n+1):
      rare=sum(w[a]*(1-F(comb(n-counts[a],m) if m<=n-counts[a] else 0,comb(n,m))) for a in range(3))
      direct=sum((coverage(sub) for sub in combinations(ids,m)),F(0))/comb(n,m)
      assert rare==direct
      checks+=1
    for s in [F(1,3),F(2,3),F(4,3),F(7,3)]:
      direct=sum(w[a]*(1-(1-s)**counts[a]) for a in range(3))
      C=[sum(w[a]*comb(counts[a],r) for a in range(3) if counts[a]>=r) for r in range(n+1)]
      collision=sum((-1)**(r+1)*s**r*C[r] for r in range(1,n+1))
      z=s-1
      incidence=coverage(ids)+sum((-1)**(r+1)*z**r*phi(ids,r) for r in range(1,n+1))
      assert direct==collision==incidence
      checks+=1
for k in range(1,5):
  for L in range(k,12):
    for psi in [lambda n:F(n),lambda n:F(n*n),lambda n:F(1,n),lambda n:F(2)**n]:
      expected=F(0)
      for m in range(k,L+1):
        h=psi(m) if m==k else psi(m)+F(m-k,k)*(psi(m)-psi(m-1))
        expected+=h*F(comb(m-1,k-1),comb(L,k))
      assert expected==psi(L)
      checks+=1
print('Exact rational checks passed:', checks)
print('Overlapping-feature mean coverage:', Fs)
