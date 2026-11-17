"""Log-domain deletion calibrators for ordinary iid, single-label draws.

These are numerical research prototypes, not interval-arithmetic certificates.
Directions and mixture weights must be fixed independently of the observations.
Counts describe observed species only: every count must be positive.
"""
from math import comb, exp, isfinite, log, log1p, lgamma
from numbers import Integral

NEG_INF = float("-inf")


def logadd(a, b):
    """Stable log(exp(a)+exp(b)), including zero terms."""
    if a == NEG_INF:
        return b
    if b == NEG_INF:
        return a
    hi, lo = max(a, b), min(a, b)
    return hi + log1p(exp(lo - hi))


def logsum(values):
    ans = NEG_INF
    for value in values:
        ans = logadd(ans, value)
    return ans


def log_poisson_mass(k, effort):
    if not isinstance(k, Integral) or k < 0:
        raise ValueError("k must be a nonnegative integer")
    if not isfinite(effort) or effort <= 0:
        raise ValueError("effort must be finite and positive")
    return -effort + k * log(effort) - lgamma(k + 1)


def _validate(counts, direction):
    counts, direction = tuple(counts), tuple(direction)
    if any(not isinstance(x, Integral) or x <= 0 for x in counts):
        raise ValueError("observed counts must be positive integers")
    if any(not isfinite(v) for v in direction):
        raise ValueError("direction must be finite")
    return counts, direction


def _score(j, direction):
    return direction[j - 1] if 1 <= j <= len(direction) else 0.0


def fingerprint_score(counts, direction):
    counts, direction = _validate(counts, direction)
    return sum(_score(x, direction) for x in counts)


def log_subset_coefficients(counts, direction):
    """Log coefficients A[k] = sum_{|I|=k} exp(lambda . Phi(S_I)).

    Sequential convolution uses O(n^2) arithmetic operations and O(n) memory.
    direction[r-1] is lambda_r; lambda_0 and undeclared directions are zero.
    """
    counts, direction = _validate(counts, direction)
    coeff = [0.0]
    for x in counts:
        factor = [log(comb(x, j)) + _score(j, direction)
                  for j in range(x + 1)]
        out = [NEG_INF] * (len(coeff) + x)
        for k, a in enumerate(coeff):
            for j, b in enumerate(factor):
                out[k + j] = logadd(out[k + j], a + b)
        coeff = out
    return coeff


def log_lost_species_coefficients(counts, direction):
    """Log A[k][r], additionally tracking r species lost by deletion.

    A is the coefficient array of product_i
    (w + sum_{j=1}^{X_i} binom(X_i,j) exp(lambda_j) z^j).
    Dense bounds: O(D*n^2) operations, O(n*D) memory.
    """
    counts, direction = _validate(counts, direction)
    coeff = [[0.0]]
    for species_index, x in enumerate(counts):
        out = [[NEG_INF] * (species_index + 2)
               for _ in range(len(coeff) + x)]
        factor = [log(comb(x, j)) + _score(j, direction)
                  for j in range(1, x + 1)]
        for k, row in enumerate(coeff):
            for r, a in enumerate(row):
                if a == NEG_INF:
                    continue
                out[k][r + 1] = logadd(out[k][r + 1], a)
                for j, b in enumerate(factor, 1):
                    out[k + j][r] = logadd(out[k + j][r], a + b)
        coeff = out
    return coeff


def log_deletion_score(counts, direction, effort):
    """Log H for dominated Poisson deletion; this is not yet an e-value.

    Subtract the candidate population's Poisson cumulant transform to
    obtain a log e-value. Do not renormalize the retained Poisson weights.
    """
    counts, direction = _validate(counts, direction)
    n = sum(counts)
    coeff = log_subset_coefficients(counts, direction)
    return logsum(log_poisson_mass(k, effort) + a - log(comb(n, k))
                  for k, a in enumerate(coeff))


def log_prediction_deletion_score(counts, direction, eta, future_size,
                                  effort, future_effort):
    """Log H for the transfer draft's fixed-sample novelty sandwich, B=1.

    The final candidate log e-value is H + eta*y - K_{lambda,eta}(nu).
    Its validity is joint/unconditional, not conditional on the past sample.
    """
    counts, direction = _validate(counts, direction)
    if not isinstance(future_size, Integral) or future_size < 0:
        raise ValueError("future_size must be a nonnegative integer")
    if not isfinite(eta):
        raise ValueError("eta must be finite")
    n, m = sum(counts), future_size
    if eta >= 0:
        future = logsum(log_poisson_mass(ell, future_effort)
                        - eta * (m - ell) for ell in range(m + 1))
        return log_deletion_score(counts, direction, effort) + future
    coeff = log_lost_species_coefficients(counts, direction)
    future = [logsum(log_poisson_mass(ell, future_effort)
                     + eta * min(r, ell) for ell in range(m + 1))
              for r in range(len(counts) + 1)]
    return logsum(log_poisson_mass(k, effort) - log(comb(n, k))
                  + a + future[r]
                  for k, row in enumerate(coeff)
                  for r, a in enumerate(row) if a != NEG_INF)


def optimized_full_count_log_e(log_score, sample_sizes, effort_grid,
                               cumulant):
    """Optimize deterministic normalizers for a COMMON fixed-count score.

    cumulant(efforts) supplies K at a fixed candidate population and direction.
    Returns (log_e, selected_efforts). A finite, predeclared, nonempty grid
    makes the optimization exact over that grid, up to floating-point error.

    This does NOT license maximizing log_deletion_score over effort:
    the latter changes the data score, so the common-score proof does not apply.
    """
    sample_sizes = tuple(sample_sizes)
    if not sample_sizes or any(not isinstance(n, Integral) or n < 0
                               for n in sample_sizes):
        raise ValueError("sample sizes must be nonnegative integers")
    if not isfinite(log_score):
        raise ValueError("log_score must be finite")
    best, selected = float("inf"), None
    for efforts in effort_grid:
        efforts = tuple(efforts)
        if len(efforts) != len(sample_sizes):
            raise ValueError("each effort tuple must match the sample sizes")
        penalty = -sum(log_poisson_mass(n, t)
                       for n, t in zip(sample_sizes, efforts))
        value = cumulant(efforts) + penalty
        if not isfinite(value):
            raise ValueError("cumulant and normalizer must be finite")
        if value < best:
            best, selected = value, efforts
    if selected is None:
        raise ValueError("effort grid must be nonempty")
    return log_score - best, selected
