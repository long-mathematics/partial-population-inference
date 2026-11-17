"""Explicit approximate calibration, separate from continuum certification.

An estimator, its variance/bias model, and its calibration are independent
choices. The generic interface permits other estimators and studentizations;
it does not upgrade a caller's assumptions into a coverage theorem.
"""
import math
from scipy.stats import norm

def normal_interval(estimate,variance,*,target,assumptions,alpha=.05,
                    support=None,bias_bound=0.):
    if not target or not assumptions:
        raise ValueError('state the target and approximation assumptions')
    if not all(math.isfinite(x) for x in (estimate,variance,alpha,bias_bound)):
        raise ValueError('finite numerical inputs required')
    if variance<0 or bias_bound<0 or not 0<alpha<1:
        raise ValueError('invalid variance, bias bound or level')
    radius=float(norm.ppf(1-alpha/2))*math.sqrt(variance)+bias_bound
    raw=[estimate-radius,estimate+radius];interval=raw.copy()
    if support is not None:
        lo,hi=support
        if not lo<=hi:raise ValueError('invalid support')
        interval=[max(lo,raw[0]),min(hi,raw[1])]
        if interval[0]>interval[1]:raise ValueError('normal interval misses declared support')
    return {'estimate':estimate,'variance':variance,'raw_interval':raw,'interval':interval,
            'target':target,'alpha':alpha,'support':support,'bias_bound':bias_bound,
            'calibration':'normal approximation; no general finite-sample guarantee',
            'assumptions':assumptions}

def good_turing_normal(n,phi1,phi2,*,variance_model='fixed_n',alpha=.05,clip=True):
    if type(n) is not int or n<1 or any(type(v) is not int or v<0 for v in (phi1,phi2)) or phi1+2*phi2>n:
        raise ValueError('invalid known sample size or fingerprint')
    if variance_model not in ('fixed_n','poisson'):
        raise ValueError('variance_model must be fixed_n or poisson')
    v=(phi1+2*phi2-(phi1**2/n if variance_model=='fixed_n' else 0))/n**2
    assumptions=('Fixed n iid species draws; an applicable triangular-array missing-mass CLT, '
                 'nondegenerate studentization, negligible standardized bias, and no dominant rare-event contribution.'
                 if variance_model=='fixed_n' else
                 'Independent Poisson species counts at known effort n; a Lindeberg-type CLT and consistent studentization. '
                 'Using this variance at fixed count is an approximation, not an exact transfer.')
    result=normal_interval(phi1/n,max(0.,v),target='realized missing mass',assumptions=assumptions,
                           alpha=alpha,support=[0.,1.] if clip else None)
    result['variance_model']=variance_model
    result['degenerate_variance']=v==0
    return result
