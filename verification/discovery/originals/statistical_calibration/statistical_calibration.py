#!/usr/bin/env python3
"""Finite-sample checks for statistical calibration and optimality.

Requires Python >=3.10, numpy, scipy. Run this file to write verification.json.

The fingerprint distribution is enumerated exactly up to floating-point
arithmetic: all per-species counts outside {1,...,R} are aggregated into one
category, not truncated. No Monte Carlo is used.

The confidence-set LP optimizes ONLY the explicitly declared finite experiment.
It is not a claim of minimax optimality for unrestricted populations.
The exponential calibrators are fixed before seeing data.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction
from math import comb, log
from pathlib import Path
from typing import Iterable
import json

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import brentq, linprog
from scipy.special import logsumexp, xlogy
from scipy.stats import poisson

Array = NDArray[np.float64]


def checked_p(p: Iterable[float]) -> Array:
    out = np.asarray(list(p), dtype=float)
    if out.ndim != 1 or len(out) == 0 or np.any(out <= 0):
        raise ValueError('p must be a nonempty positive probability vector')
    if not np.isclose(out.sum(), 1., rtol=0., atol=1e-12):
        raise ValueError('probabilities must sum to one')
    return out


def occupancy_probabilities(p: Iterable[float], t: float, R: int) -> Array:
    p = checked_p(p)
    if t <= 0 or R < 1:
        raise ValueError('require t>0 and R>=1')
    return poisson.pmf(np.arange(1, R+1)[None, :], t*p[:, None])


def fingerprint_distribution(p: Iterable[float], t: float, R: int = 2
                             ) -> tuple[Array, Array, Array]:
    """Return states, exact state probabilities (numerically), and E Phi."""
    b = occupancy_probabilities(p, t, R)
    distribution = {(0,)*R: 1.}
    for row in b:
        probabilities = np.r_[1-row.sum(), row]
        new: dict[tuple[int, ...], float] = {}
        for state, mass in distribution.items():
            for category, probability in enumerate(probabilities):
                if probability == 0:
                    continue
                s = list(state)
                if category:
                    s[category-1] += 1
                key = tuple(s)
                new[key] = new.get(key, 0.) + mass*float(probability)
        distribution = new
    states = np.asarray(list(distribution), dtype=float)
    weights = np.asarray(list(distribution.values()), dtype=float)
    if abs(weights.sum()-1) > 1e-11:
        raise ArithmeticError('enumerated distribution does not sum to one')
    return states, weights, b.sum(axis=0)


def exact_log_mgf(p: Iterable[float], t: float, lambdas: Array) -> float:
    """Log E exp(lambda dot Phi), using independence across species."""
    lambdas = np.asarray(lambdas, dtype=float)
    b = occupancy_probabilities(p, t, len(lambdas))
    return float(np.log1p(b @ np.expm1(lambdas)).sum())


def calibration_kernel(x: Array, t: float, lambdas: Array) -> Array:
    """kappa_lambda(x); integrates against sum_i p_i delta_{p_i}."""
    x = np.asarray(x, dtype=float)
    if np.any(x <= 0):
        raise ValueError('kernel routine requires x>0')
    lambdas = np.asarray(lambdas, dtype=float)
    b = poisson.pmf(np.arange(1, len(lambdas)+1)[None, :], t*x[:, None])
    return np.log1p(b @ np.expm1(lambdas))/x


def poisson_deviance(observed: Array, means: Array) -> Array:
    """h(z,m)=z log(z/m)-z+m, with h(0,m)=m."""
    observed, means = np.broadcast_arrays(observed, means)
    with np.errstate(divide='ignore', invalid='ignore'):
        result = xlogy(observed, observed/means)-observed+means
    return np.where((observed == 0) & (means == 0), 0., result)


def mean_confidence_interval(observed: int, alpha: float) -> tuple[float, float]:
    """Two-sided >=1-alpha interval for a Poisson-binomial mean.

    This is Chernoff inversion, not the exact shortest confidence interval.
    """
    if observed < 0 or not (0 < alpha < 1):
        raise ValueError('require nonnegative observation and 0<alpha<1')
    L = log(2/alpha)
    if observed == 0:
        return 0., L
    def equation(mean: float) -> float:
        return observed*(log(observed)-log(mean))-observed+mean-L
    lower = brentq(equation, np.finfo(float).tiny, float(observed))
    upper_bracket = float(observed+L+1)
    while equation(upper_bracket) < 0:
        upper_bracket *= 2
    upper = brentq(equation, float(observed), upper_bracket)
    return lower, upper


def mixture_log_evalues(states: Array, cgfs: Array, prior: Array,
                        directions: Array) -> Array:
    if np.any(prior <= 0) or abs(prior.sum()-1) > 1e-12:
        raise ValueError('mixture weights must be positive and sum to one')
    return logsumexp(states @ directions.T-cgfs[None, :]
                     + np.log(prior)[None, :], axis=1)


def tank_upper(maximum: int, sample_size: int,
               alpha: Fraction = Fraction(1, 20)) -> int:
    """Smallest nondecreasing nonrandomized (1-alpha) upper limit.

    Uses exact integer comparisons: retain N when F_N(maximum)>alpha.
    """
    if not (1 <= sample_size <= maximum and 0 < alpha < 1):
        raise ValueError('invalid tank parameters')
    numerator = comb(maximum, sample_size)*alpha.denominator
    def accepted(N: int) -> bool:
        return numerator > alpha.numerator*comb(N, sample_size)
    lo, hi = maximum, maximum+1
    while accepted(hi):
        hi *= 2
    while hi-lo > 1:
        mid = (lo+hi)//2
        if accepted(mid):
            lo = mid
        else:
            hi = mid
    return lo


def finite_tank_confidence_lp(k: int = 3, Nmax: int = 12,
                              alpha: float = .05) -> dict:
    """Exact finite-game formulation for minimax expected interval length.

    Parameter values N=k,...,Nmax; observations M=k,...,Nmax; actions are
    all closed intervals with endpoints in that parameter grid, plus empty.
    Randomization is permitted. LP numerical errors are reported.
    """
    parameters = np.arange(k, Nmax+1)
    observations = parameters.copy()
    actions = [None] + [(int(a), int(b)) for a in parameters
                       for b in parameters if b >= a]
    J, O, A = len(parameters), len(observations), len(actions)
    P = np.zeros((J, O))
    for j, N in enumerate(parameters):
        for o, m in enumerate(observations):
            if m <= N:
                P[j, o] = comb(int(m-1), k-1)/comb(int(N), k)
    coverage = np.array([[0. if I is None else float(I[0] <= N <= I[1])
                         for I in actions] for N in parameters])
    costs = np.array([0. if I is None else float(I[1]-I[0]) for I in actions])
    nvars = O*A+1
    c = np.zeros(nvars); c[-1] = 1.
    ub, rhs = [], []
    for j in range(J):
        riskrow = np.zeros(nvars)
        riskrow[:-1] = (P[j, :, None]*costs[None, :]).ravel()
        riskrow[-1] = -1.
        ub.append(riskrow); rhs.append(0.)
        covrow = np.zeros(nvars)
        covrow[:-1] = -(P[j, :, None]*coverage[j, None, :]).ravel()
        ub.append(covrow); rhs.append(-(1-alpha))
    eq = np.zeros((O, nvars))
    for o in range(O):
        eq[o, o*A:(o+1)*A] = 1.
    ub, rhs = np.array(ub), np.array(rhs)
    solution = linprog(c, A_ub=ub, b_ub=rhs, A_eq=eq, b_eq=np.ones(O),
                       bounds=(0, None), method='highs')
    if not solution.success:
        raise RuntimeError(solution.message)
    rule = solution.x[:-1].reshape(O, A)
    achieved_coverage = np.array([(P[j, :, None]*rule*coverage[j]).sum()
                                 for j in range(J)])
    risks = np.array([(P[j, :, None]*rule*costs).sum() for j in range(J)])
    dual = float(rhs @ solution.ineqlin.marginals
                 + np.ones(O) @ solution.eqlin.marginals)
    assert achieved_coverage.min() >= 1-alpha-2e-9
    assert risks.max() <= solution.fun+2e-9
    assert abs(dual-solution.fun) < 2e-9
    def rational(value: float) -> Fraction:
        if abs(value) < 1e-10:
            return Fraction(0)
        return Fraction(float(value)).limit_denominator(10**8)
    qr = [[rational(v) for v in row] for row in rule]
    yr = [rational(v) for v in solution.ineqlin.marginals]
    zr = [rational(v) for v in solution.eqlin.marginals]
    RR = rational(solution.fun)
    ar = rational(alpha)
    Pr = [[Fraction(comb(int(m-1), k-1), comb(int(N), k)) if m <= N
           else Fraction(0) for m in observations] for N in parameters]
    for row in qr:
        assert all(v >= 0 for v in row) and sum(row) == 1
    for j, N in enumerate(parameters):
        cov = sum(Pr[j][o]*qr[o][a]
                  for o in range(O) for a,I in enumerate(actions)
                  if I is not None and I[0] <= N <= I[1])
        risk = sum(Pr[j][o]*qr[o][a]*int(costs[a])
                   for o in range(O) for a in range(A))
        assert cov >= 1-ar and risk <= RR
    assert all(v <= 0 for v in yr)
    assert -sum(yr[2*j] for j in range(J)) <= 1
    for o in range(O):
        for a,I in enumerate(actions):
            lhs = zr[o] + sum(Pr[j][o]*(int(costs[a])*yr[2*j]
                        - int(I is not None and I[0] <= N <= I[1])*yr[2*j+1])
                        for j,N in enumerate(parameters))
            assert lhs <= 0
    rational_dual = -(1-ar)*sum(yr[2*j+1] for j in range(J))+sum(zr)
    assert rational_dual == RR
    certificate = {'description': 'Exact rational primal/dual certificate for the finite interval game',
        'parameters': parameters.tolist(), 'observations': observations.tolist(),
        'actions': actions, 'alpha': str(ar), 'optimal_value': str(RR),
        'decision_probabilities': [[str(v) for v in row] for row in qr],
        'inequality_dual': [str(v) for v in yr], 'equality_dual': [str(v) for v in zr]}
    Path(__file__).with_name('finite_game_certificate.json').write_text(
        json.dumps(certificate, indent=2)+'\n')
    return {'parameter_range': [k, Nmax], 'sample_size': k, 'alpha': alpha,
            'exact_rational_optimum': str(RR),
            'rational_primal_dual_verified': True,
            'randomization_allowed': True, 'empty_intervals_allowed': True,
            'minimum_coverage': float(achieved_coverage.min()),
            'optimal_worst_expected_length': float(solution.fun),
            'dual_value': dual, 'numerical_duality_gap': abs(dual-solution.fun)}


def run_checks() -> dict:
    rng = np.random.default_rng(20261112)
    R, alpha = 2, .05
    directions = np.array([[a, b] for a in (-2., -.5, 0., .5, 2.)
                           for b in (-2., -.5, 0., .5, 2.)
                           if a or b])
    prior = np.full(len(directions), 1/len(directions))
    cases = []
    errors = {'mgf_identity': 0., 'measure_kernel_identity': 0.,
              'variance_identity': 0., 'evalue_mean_identity': 0.}
    counts = {'populations_and_efforts': 0, 'mgf_directions': 0,
              'integer_tank_coverage': 0, 'tank_minimality_witnesses': 0,
              'two_point_minimax': 0}
    populations = [[.5, .3, .2], [.25]*4, [1/10]*10,
                   [.7, .1, .1, .1], [1.]]
    populations += [rng.dirichlet(np.ones(n)).tolist() for n in (3, 5, 8)]
    for p in populations:
        p = checked_p(p)
        for t in (.5, 3., 15.):
            states, probability, means = fingerprint_distribution(p, t, R)
            cgfs = np.array([exact_log_mgf(p, t, lam) for lam in directions])
            envelopes = np.expm1(directions) @ means
            for lam, cgf, envelope in zip(directions, cgfs, envelopes):
                actual = float(logsumexp(np.log(probability)+states @ lam))
                errors['mgf_identity'] = max(errors['mgf_identity'], abs(actual-cgf))
                represented = float(p @ calibration_kernel(p, t, lam))
                errors['measure_kernel_identity'] = max(
                    errors['measure_kernel_identity'], abs(represented-cgf))
                assert cgf <= envelope+2e-11
                counts['mgf_directions'] += 1
            logE = mixture_log_evalues(states, cgfs, prior, directions)
            logE_bound = mixture_log_evalues(states, envelopes, prior, directions)
            exact_e_mean = float(probability @ np.exp(logE))
            bound_e_mean = float(probability @ np.exp(logE_bound))
            errors['evalue_mean_identity'] = max(errors['evalue_mean_identity'],
                                               abs(exact_e_mean-1))
            assert bound_e_mean <= 1+2e-11
            assert np.all(logE >= logE_bound-2e-11)
            cov_exact = float(probability[logE < log(1/alpha)].sum())
            cov_bound = float(probability[logE_bound < log(1/alpha)].sum())
            h = poisson_deviance(states, means)
            cov_box = float(probability[np.all(h <= log(2*R/alpha), axis=1)].sum())
            assert min(cov_exact, cov_bound, cov_box) >= 1-alpha-2e-11
            coefficients = np.array([.4, -.7])
            statistic = states @ coefficients
            actual_var = float(probability @ (statistic-probability@statistic)**2)
            b = occupancy_probabilities(p, t, R)
            formula_var = float(np.sum(b @ coefficients**2-(b@coefficients)**2))
            errors['variance_identity'] = max(errors['variance_identity'],
                                              abs(actual_var-formula_var))
            cases.append({'p': p.tolist(), 't': t, 'fingerprint_states': len(states),
                          'exact_mgf_mixture_coverage': cov_exact,
                          'poisson_envelope_mixture_coverage': cov_bound,
                          'deviance_box_coverage': cov_box,
                          'poisson_envelope_evalue_expectation': bound_e_mean})
            counts['populations_and_efforts'] += 1
    max_tank_error = Fraction(0)
    for k in range(1, 7):
        for N in range(k, 81):
            error_numerator = sum(comb(m-1, k-1) for m in range(k, N+1)
                                  if tank_upper(m, k) < N)
            error = Fraction(error_numerator, comb(N, k))
            assert error <= Fraction(1,20)
            max_tank_error = max(max_tank_error, error)
            counts['integer_tank_coverage'] += 1
        for m in range(k, 41):
            U = tank_upper(m, k)
            assert Fraction(comb(m,k), comb(U,k)) > Fraction(1,20)
            assert Fraction(comb(m,k), comb(U+1,k)) <= Fraction(1,20)
            counts['tank_minimality_witnesses'] += 1
    two_point_max_gap = 0.
    for common_probability in np.linspace(.01, 1., 31):
        r = float(common_probability)
        # Estimate 1+d on the common observation; estimate 2 otherwise.
        sol = linprog([0., 1.], A_ub=[[1., -1.], [-r, -1.]],
                      b_ub=[0., -r], bounds=[(0.,1.), (0.,None)], method='highs')
        if not sol.success:
            raise RuntimeError(sol.message)
        gap = abs(sol.fun-r/(1+r))
        assert gap < 1e-10
        two_point_max_gap = max(two_point_max_gap, gap)
        counts['two_point_minimax'] += 1
    from itertools import product
    missing_and_prediction_checks = 0
    maximum_missing_variance_error = 0.
    maximum_prediction_mgf_error = 0.
    for pp in ([.5,.3,.2], [.25]*4):
        pp = np.array(pp)
        for time in (.5, 3., 12.):
            bb = occupancy_probabilities(pp, time, 2)
            b0 = np.exp(-time*pp)
            # Exact GT-minus-realized-missing-mass error: -p, 1/t, or 0.
            vals = np.column_stack((-pp, np.full(len(pp), 1/time), np.zeros(len(pp))))
            probs = np.column_stack((b0, bb[:,0], 1-b0-bb[:,0]))
            mean, second = 0., 0.
            for state in product(range(3), repeat=len(pp)):
                mass = float(np.prod([probs[i,j] for i,j in enumerate(state)]))
                value = sum(vals[i,j] for i,j in enumerate(state))
                mean += mass*value; second += mass*value*value
            expected_variance = (bb[:,0].sum()+2*bb[:,1].sum())/time**2
            maximum_missing_variance_error = max(maximum_missing_variance_error,
                                                  abs(second-expected_variance))
            assert abs(mean) < 1e-12 and abs(second-expected_variance) < 1e-12
            future = 1.5*time
            target_chance = b0*(-np.expm1(-future*pp))
            coeff = np.array([.8,-.6])
            # Error values -1 (a new future species), c1, c2, or 0.
            error_vals = np.array([-1., coeff[0], coeff[1], 0.])
            error_probs = np.column_stack((target_chance, bb[:,0], bb[:,1],
                                          1-target_chance-bb.sum(axis=1)))
            for lam in (-1., .5, 1.):
                direct = 0.
                for state in product(range(4), repeat=len(pp)):
                    mass = float(np.prod([error_probs[i,j] for i,j in enumerate(state)]))
                    value = sum(error_vals[j] for j in state)
                    direct += mass*np.exp(lam*value)
                formula = float(np.prod(1+bb@np.expm1(lam*coeff)
                                        + target_chance*np.expm1(-lam)))
                error = abs(direct-formula)
                assert error < 2e-12
                maximum_prediction_mgf_error = max(maximum_prediction_mgf_error,error)
                missing_and_prediction_checks += 1
    errors['missing_mass_variance'] = maximum_missing_variance_error
    errors['prediction_error_mgf'] = maximum_prediction_mgf_error
    counts['prediction_error_mgf_enumerations'] = missing_and_prediction_checks
    assert max(errors.values()) < 1e-10
    return {'scope': 'Exact finite-distribution enumeration, with floating-point arithmetic; no Monte Carlo.',
            'counts': counts, 'maximum_numerical_identity_errors': errors,
            'minimum_coverage_across_24_test_cases': {
                key: min(case[key] for case in cases) for key in
                ('exact_mgf_mixture_coverage', 'poisson_envelope_mixture_coverage',
                 'deviance_box_coverage')},
            'tank_largest_exact_noncoverage': str(max_tank_error),
            'two_point_maximum_lp_gap': two_point_max_gap,
            'finite_confidence_game': finite_tank_confidence_lp(),
            'cases': cases}


if __name__ == '__main__':
    report = run_checks()
    destination = Path(__file__).with_name('verification.json')
    destination.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k != 'cases'}, indent=2))
    print(f'Full report: {destination}')
