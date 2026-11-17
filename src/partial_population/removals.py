"""Uniform repeated capture with known removals of captured individuals.

N0 is the initial population size. There are no unrecorded losses, recruitment,
tag errors, or unequal catchability. Each observed removal occurs after that
occasion's capture; the next sample is uniform among the remaining individuals.
"""
from fractions import Fraction
from math import inf
from scipy.optimize import brentq
from .capture_recapture import hypergeom_pmf

def removal_states(sizes,overlaps,removed_after,initial_marked=0):
    sizes,overlaps,removed_after=map(tuple,(sizes,overlaps,removed_after))
    if not sizes or len(sizes)!=len(overlaps) or len(removed_after)!=len(sizes)-1:
        raise ValueError('one overlap per catch and one removal count between catches required')
    if type(initial_marked) is not int or initial_marked<0:raise ValueError('invalid initial marked count')
    M,D=initial_marked,0;out=[]
    for t,(C,R) in enumerate(zip(sizes,overlaps)):
        if type(C) is not int or C<1 or type(R) is not int or not 0<=R<=min(M,C):
            raise ValueError('invalid capture or recapture count')
        out.append((M,D,C,R))
        if t<len(removed_after):
            d=removed_after[t]
            if type(d) is not int or not 0<=d<=C:raise ValueError('removals must be captured individuals')
            M+=C-R-d;D+=d
    return tuple(out)

def removal_evalue(N0,states,factor):
    """Exact rational constant-tilt test martingale at the final occasion."""
    if type(N0) is not int or N0<1:raise ValueError('invalid initial population size')
    z=Fraction(factor)
    if z<=0:raise ValueError('positive factor required')
    value=Fraction(1)
    for M,D,C,R in states:
        if N0<D+M+C-R:return inf
        p=hypergeom_pmf(N0-D,M,C)
        if not p[R]:return inf
        value*=z**R/sum((prob*z**r for r,prob in enumerate(p)),Fraction(0))
    return value

def removal_interval(sizes,overlaps,removed_after,*,initial_marked=0,alpha=Fraction(1,20),factor=Fraction(2)):
    """Exact two-tail e-value interval for N0, with an unbounded integer search.

Constant factors z and 1/z are fixed before observing outcomes. Their
normalizers are monotone in N0, so doubling plus bisection is exact. None is
an infinite upper endpoint. The coverage claim is under removal_states' model.
"""
    states=removal_states(sizes,overlaps,removed_after,initial_marked)
    alpha,z=Fraction(alpha),Fraction(factor)
    if not 0<alpha<1 or z<=1:raise ValueError('require 0<alpha<1 and factor>1')
    support=max(D+M+C-R for M,D,C,R in states);threshold=2/alpha
    def accepted_low(N):return removal_evalue(N,states,1/z)<=threshold
    lo=support
    if not accepted_low(lo):
        high=max(lo+1,2*lo)
        while not accepted_low(high):high*=2
        left=lo
        while left+1<high:
            mid=(left+high)//2
            if accepted_low(mid):high=mid
            else:left=mid
        lo=high
    # E_z increases towards z^(sum R); a finite upper bound exists only if
    # that limit strictly exceeds the threshold (equality is still accepted).
    if z**sum(s[3] for s in states)<=threshold:hi=None
    elif removal_evalue(support,states,z)>threshold:hi=support-1
    else:
        left=support;right=max(support+1,2*support)
        while removal_evalue(right,states,z)<=threshold:right*=2
        while left+1<right:
            mid=(left+right)//2
            if removal_evalue(mid,states,z)<=threshold:left=mid
            else:right=mid
        hi=left
    return [support,None] if hi is not None and lo>hi else [lo,hi]

def removal_estimates(sizes,overlaps,removed_after,initial_marked=0):
    """Historical constant-N ratios and roots of removal-aware unbiased scores.

The score has expectation zero at N0; its root is not claimed unbiased or
exactly calibrated. Marked-pool updates use all recorded removals in both.
"""
    states=removal_states(sizes,overlaps,removed_after,initial_marked)
    support=max(D+M+C-R for M,D,C,R in states);out={}
    for label,weighted in [('schnabel',False),('schumacher_eschmeyer',True)]:
        denom=sum((M if weighted else 1)*R for M,D,C,R in states)
        num=sum((M if weighted else 1)*M*C for M,D,C,R in states)
        if not denom:raise ValueError('positive recapture information required')
        out[label+'_constant_population']=float(Fraction(num,denom))
        def score(N):return sum((M if weighted else 1)*M*C/(N-D) for M,D,C,R in states)-denom
        high=max(support+1,2*support)
        while score(high)>0:high*=2
        out[label+'_initial_population_score_root']=(float(support) if score(support)<=0 else float(brentq(score,support,high)))
    return out
