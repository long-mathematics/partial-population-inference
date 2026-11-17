"""Exact small-design capture–recapture calculations.

All probabilities, finite estimates, and confidence thresholds use Fraction.
None denotes an infinite confidence endpoint; inf denotes an undefined ratio
estimate at zero recaptures. Labels carry identity only, not numerical ranks.
These routines favor exactness over speed for large populations/sample sizes.
"""
from fractions import Fraction
from math import comb, inf, prod


def _integer(x, name, minimum=0):
    if isinstance(x, bool) or not isinstance(x, int) or x < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return x


def _rational(x, name):
    if isinstance(x, bool) or not isinstance(x, (int, Fraction)):
        raise ValueError(f"{name} must be an int or Fraction, not a float")
    return Fraction(x)


def _alpha(alpha):
    alpha = _rational(alpha, "alpha")
    if not 0 < alpha < 1:
        raise ValueError("alpha must lie strictly between zero and one")
    return alpha


def choose(n, k):
    """Binomial coefficient, zero outside 0 <= k <= n."""
    return comb(n, k) if 0 <= k <= n else 0


def falling(n, r):
    return prod(range(n-r+1, n+1))


def rising(n, r):
    return prod(range(n, n+r))


def _fixed_data(a, b, j):
    _integer(a, "a", 1)
    _integer(b, "b", 1)
    _integer(j, "j")
    if j > min(a, b):
        raise ValueError("overlap exceeds a sample size")


def hypergeom_pmf(N, a, b):
    """P(J=j), j=0,...,min(a,b), including impossible zero cells."""
    _integer(N, "N", 1)
    _integer(a, "a")
    _integer(b, "b")
    if max(a, b) > N:
        raise ValueError("sample/marked count exceeds population")
    return tuple(Fraction(choose(a,j)*choose(N-a,b-j), choose(N,b))
                 for j in range(min(a,b)+1))


def lincoln_petersen(a, b, j):
    _fixed_data(a, b, j)
    return Fraction(a*b, j) if j else inf


def chapman(a, b, j):
    _fixed_data(a, b, j)
    return Fraction((a+1)*(b+1), j+1)-1


def chapman_moments(N, a, b):
    """Exact mean, variance, and bias; no plug-in variance approximation."""
    _fixed_data(a, b, 0)
    p = hypergeom_pmf(N, a, b)
    mean = sum(v*chapman(a,b,j) for j,v in enumerate(p))
    variance = sum(v*(chapman(a,b,j)-mean)**2 for j,v in enumerate(p))
    return mean, variance, mean-N


def shifted_reciprocal(a, b, j, k=1):
    _fixed_data(a, b, j)
    _integer(k, "k", 1)
    return Fraction(rising(a+1,k)*rising(b+1,k), rising(j+1,k))


def binomial_pmf(N, a, b):
    """With-replacement marked-count law, j=0,...,b."""
    _integer(N, "N", 1)
    _integer(a, "a", 1)
    _integer(b, "b", 1)
    if a > N:
        raise ValueError("marked count exceeds population")
    p = Fraction(a,N)
    return tuple(comb(b,j)*p**j*(1-p)**(b-j) for j in range(b+1))


def bailey(a, b, j):
    """Binomial second-sample correction; repeats of marked units count."""
    _integer(a, "a", 1)
    _integer(b, "b", 1)
    _integer(j, "j")
    if j > b:
        raise ValueError("recaptures exceed draws")
    return Fraction(a*(b+1), j+1)


def _first_true(start, predicate):
    if predicate(start):
        return start
    lo, hi = start, max(start+1, 2*start)
    while not predicate(hi):
        lo, hi = hi, 2*hi
    while hi-lo > 1:
        mid = (lo+hi)//2
        if predicate(mid):
            hi = mid
        else:
            lo = mid
    return hi


def fixed_lower_limit(a, b, j, alpha=Fraction(1,20)):
    """Largest nonincreasing one-sided lower limit at error alpha."""
    _fixed_data(a,b,j)
    alpha = _alpha(alpha)
    return _first_true(a+b-j,
                       lambda N: sum(hypergeom_pmf(N,a,b)[:j+1]) > alpha)


def fixed_upper_limit(a, b, j, alpha=Fraction(1,20)):
    """Smallest nonincreasing one-sided upper limit; None at j=0."""
    _fixed_data(a,b,j)
    alpha = _alpha(alpha)
    if j == 0:
        return None
    return _first_true(a+b-j,
                       lambda N: sum(hypergeom_pmf(N,a,b)[j:]) <= alpha)-1


def fixed_interval(a, b, j, alpha=Fraction(1,20)):
    """Equal-tail interval with coverage >= 1-alpha on unbounded integer N.

    Strict acceptance inequalities are compared as rational numbers. A pair
    with lower > upper represents the empty interval; no endpoint is clipped
    to an arbitrary search ceiling.
    """
    alpha = _alpha(alpha)
    return (fixed_lower_limit(a,b,j,alpha/2),
            fixed_upper_limit(a,b,j,alpha/2))


def signed_transfer(a, b, c, d):
    """Unique normalized signed kernel (a,b)->(c,d) for unknown N.

    The common parameter space is every N >= max(a,b,c,d). A transfer exists
    iff min(c,d) <= min(a,b). Nonnegative entries are needed for calibration.
    """
    for name, x in zip(("a","b","c","d"),(a,b,c,d)):
        _integer(x,name,1)
    m, h = min(a,b), min(c,d)
    if h > m:
        raise ValueError("no signed transfer on the unbounded parameter space")
    ratios = [Fraction(falling(c,r)*falling(d,r),
                       falling(a,r)*falling(b,r)) for r in range(h+1)]
    return tuple(tuple(sum((-1)**(r-y)*choose(r,y)*choose(j,r)*ratios[r]
                           for r in range(y,h+1))
                       for y in range(h+1)) for j in range(m+1))


def _sequential_data(a, r, t):
    _integer(a,"a",1)
    _integer(r,"r",1)
    _integer(t,"t",r)
    if r > a:
        raise ValueError("without-replacement recapture target exceeds marked count")


def inverse_hypergeom_pmf(N, a, r):
    _sequential_data(a,r,r)
    _integer(N,"N",a)
    return {t: Fraction(choose(t-1,r-1)*choose(N-t,a-r),choose(N,a))
            for t in range(r,N-a+r+1)}


def sequential_estimate(a, r, t):
    _sequential_data(a,r,t)
    return Fraction((a+1)*t,r)-1


def sequential_variance(N, a, r):
    _sequential_data(a,r,r)
    _integer(N,"N",a)
    return Fraction((N+1)*(N-a)*(a+1-r),r*(a+2))


def _sequential_cdf(N, a, r, t):
    if t < r:
        return Fraction(0)
    if t >= N:
        return Fraction(1)
    return sum(hypergeom_pmf(N,a,t)[r:])


def sequential_interval(a, r, t, alpha=Fraction(1,20)):
    """Exact equal-tail interval for a preassigned r, without replacement."""
    _sequential_data(a,r,t)
    alpha = _alpha(alpha)
    start = a+t-r
    lower = _first_true(start,
                        lambda N: 1-_sequential_cdf(N,a,r,t-1) > alpha/2)
    upper = _first_true(start,
                        lambda N: _sequential_cdf(N,a,r,t) <= alpha/2)-1
    return lower, upper


def inverse_binomial_estimate(a, r, t):
    """With-replacement stopping at r marked draws; r may exceed a."""
    _integer(a,"a",1)
    _integer(r,"r",1)
    _integer(t,"t",r)
    return Fraction(a*t,r)


def repeated_estimate(marked, sizes, overlaps, *, weight_marked=False):
    """Schnabel (default) or Schumacher–Eschmeyer (weight_marked=True).

    Each occasion returns its sample and marks all newly captured individuals.
    The supplied marked counts must obey the observed update relation.
    """
    marked, sizes, overlaps = tuple(marked), tuple(sizes), tuple(overlaps)
    if not marked or not len(marked)==len(sizes)==len(overlaps):
        raise ValueError("nonempty sequences must have equal lengths")
    numerator, denominator = 0, 0
    for t,(a,b,j) in enumerate(zip(marked,sizes,overlaps)):
        _integer(a,"marked")
        _integer(b,"size",1)
        _integer(j,"overlap")
        if j > min(a,b):
            raise ValueError("invalid overlap")
        if t and a != marked[t-1]+sizes[t-1]-overlaps[t-1]:
            raise ValueError("marked counts violate the mark-all update")
        w = a if weight_marked else 1
        numerator += w*a*b
        denominator += w*j
    if numerator == 0:
        raise ValueError("no marked-exposure information")
    return Fraction(numerator,denominator) if denominator else inf


def capture_evalue_path(N, initial_marked, sizes, overlaps, factors):
    """Exact repeated-occasion exponential test-martingale values.

    factors[t]=exp(lambda_t)>0 must be chosen before occasion t using only
    past observations. Passing retrospectively chosen factors is not valid.
    Impossible candidate histories return inf from their first zero cell.
    """
    _integer(N,"N",1)
    _integer(initial_marked,"initial_marked")
    sizes, overlaps, factors = tuple(sizes), tuple(overlaps), tuple(factors)
    if not len(sizes)==len(overlaps)==len(factors):
        raise ValueError("sequences must have equal lengths")
    for b,j,z in zip(sizes,overlaps,factors):
        _integer(b,"size",1)
        _integer(j,"overlap")
        if j > b:
            raise ValueError("overlap exceeds sample size")
        if _rational(z,"factor") <= 0:
            raise ValueError("factors must be positive")
    a, value, result = initial_marked, Fraction(1), []
    for b,j,z in zip(sizes,overlaps,factors):
        z = Fraction(z)
        if value == inf or max(a,b)>N or j>min(a,b):
            value = inf
        else:
            p = hypergeom_pmf(N,a,b)
            if not p[j]:
                value = inf
            else:
                normalizer = sum(v*z**k for k,v in enumerate(p))
                value *= z**j/normalizer
        result.append(value)
        a += b-j
    return tuple(result)


def draw_evalue_path(N, a, draws, factors):
    """Exact test martingale for sequential draws WITHOUT replacement.

    Each draw is a marked/unmarked bit. Predictable positive factors permit
    arbitrary data-dependent stopping; fixed-r intervals do not.
    """
    _integer(N,"N",1)
    _integer(a,"a",1)
    draws, factors = tuple(draws), tuple(factors)
    if len(draws)!=len(factors):
        raise ValueError("sequences must have equal lengths")
    for bit,z in zip(draws,factors):
        if type(bit) is not int or bit not in (0,1):
            raise ValueError("draws must be integer zero/one bits")
        if _rational(z,"factor") <= 0:
            raise ValueError("factors must be positive")
    successes, value, result = 0, Fraction(1), []
    for s,(bit,z) in enumerate(zip(draws,factors)):
        z = Fraction(z)
        if value==inf or a>N or s>=N or not 0<=a-successes<=N-s:
            value = inf
        else:
            p = Fraction(a-successes,N-s)
            chance = p if bit else 1-p
            value = value*z**bit/(1-p+p*z) if chance else inf
        result.append(value)
        successes += bit
    return tuple(result)


def horvitz_thompson(observed_inclusion_probabilities):
    """Population-total estimate from known positive inclusion probabilities."""
    probabilities = tuple(_rational(p,"inclusion probability")
                          for p in observed_inclusion_probabilities)
    if any(not 0<p<=1 for p in probabilities):
        raise ValueError("observed inclusion probabilities must lie in (0,1]")
    return sum((1/p for p in probabilities), Fraction(0))


def history_probabilities(catchabilities):
    """Independent-occasion history probabilities; bit k is occasion k."""
    q = tuple(_rational(p,"catchability") for p in catchabilities)
    if not q or any(not 0<=p<=1 for p in q):
        raise ValueError("a nonempty catchability vector in [0,1] is required")
    return tuple(prod(p if h & (1<<k) else 1-p for k,p in enumerate(q))
                 for h in range(1<<len(q)))


def history_transform(catchabilities, weights):
    """Exact E product_h weights[h]**C_h for fixed individual probabilities.

    weights[0] can calibrate realized unseen count. This product is conditional
    on the fixed finite catchability population, not a random-effects mixture.
    """
    catchabilities = tuple(tuple(q) for q in catchabilities)
    weights = tuple(_rational(w,"weight") for w in weights)
    if not weights or any(w<=0 for w in weights):
        raise ValueError("weights must be positive")
    value = Fraction(1)
    for q in catchabilities:
        probabilities = history_probabilities(q)
        if len(probabilities)!=len(weights):
            raise ValueError("weights must cover all histories of every individual")
        value *= sum(p*w for p,w in zip(probabilities,weights))
    return value


def binomial_moment_lower_bound(K, mean_singletons, mean_doubletons):
    """Lower bound on expected unseen count from POPULATION moments.

    Replacing these moments with observed counts does not give a confidence
    bound. K is the number of independent, equal-probability occasions per
    individual; probabilities may differ across individuals.
    """
    _integer(K,"K",2)
    m1 = _rational(mean_singletons,"mean_singletons")
    m2 = _rational(mean_doubletons,"mean_doubletons")
    if min(m1,m2)<0 or (m2==0 and m1!=0):
        raise ValueError("incompatible nonnegative moments")
    return Fraction(K-1,2*K)*m1*m1/m2 if m2 else Fraction(0)


def repeated_distribution(N, initial_marked, sizes, *, weight_marked=False):
    """Exact law of a repeated-capture ratio for a preassigned size sequence.

    Dynamic programming retains (marked, numerator, denominator). The returned
    law explicitly retains any atom at inf; no finite-risk claim discards it.
    """
    _integer(N,"N",1)
    _integer(initial_marked,"initial_marked")
    if initial_marked>N:
        raise ValueError("marked count exceeds population")
    sizes = tuple(sizes)
    if not sizes:
        raise ValueError("sizes must be nonempty")
    for b in sizes:
        _integer(b,"size",1)
        if b>N:
            raise ValueError("sample size exceeds population")
    states = {(initial_marked,0,0): Fraction(1)}
    for b in sizes:
        updated = {}
        for (a,num,den), probability in states.items():
            w = a if weight_marked else 1
            for j,p in enumerate(hypergeom_pmf(N,a,b)):
                if p:
                    key = (a+b-j,num+w*a*b,den+w*j)
                    updated[key] = updated.get(key,Fraction(0))+probability*p
        states = updated
    law = {}
    for (_,num,den),p in states.items():
        if num==0:
            raise ValueError("no marked-exposure information")
        estimate = Fraction(num,den) if den else inf
        law[estimate] = law.get(estimate,Fraction(0))+p
    return law


def binomial_interval(a, b, j, alpha=Fraction(1,20)):
    """Exact equal-tail interval for a fixed with-replacement second sample."""
    bailey(a,b,j)  # validate this design's observation range
    alpha = _alpha(alpha)
    lower = _first_true(a,lambda N: sum(binomial_pmf(N,a,b)[:j+1])>alpha/2)
    upper = (None if j==0 else
             _first_true(a,lambda N: sum(binomial_pmf(N,a,b)[j:])<=alpha/2)-1)
    return lower, upper


def _inverse_binomial_cdf(N, a, r, t):
    return sum(binomial_pmf(N,a,t)[r:]) if t>=r else Fraction(0)


def inverse_binomial_interval(a, r, t, alpha=Fraction(1,20)):
    """Exact equal-tail interval at a preassigned marked-draw target r."""
    inverse_binomial_estimate(a,r,t)
    alpha = _alpha(alpha)
    lower = _first_true(a,lambda N: 1-_inverse_binomial_cdf(N,a,r,t-1)>alpha/2)
    upper = _first_true(a,lambda N: _inverse_binomial_cdf(N,a,r,t)<=alpha/2)-1
    return lower, upper


def chao_incidence(K, frequencies):
    """Finite-occasion Chao-type POINT estimate, not a confidence bound.

    frequencies[j-1] counts units seen on exactly j of K occasions. Zero
    doubletons with positive singletons returns inf. If both are zero, S is
    an explicit point-estimate convention and does not identify the unseen.
    """
    _integer(K,"K",2)
    frequencies = tuple(frequencies)
    if len(frequencies)!=K:
        raise ValueError("provide frequencies for counts 1 through K")
    for count in frequencies:
        _integer(count,"frequency")
    S=sum(frequencies)
    f1,f2=frequencies[:2]
    if f2==0:
        return inf if f1 else Fraction(S)
    return S+Fraction((K-1)*f1*f1,2*K*f2)
