#!/usr/bin/env python3
"""Checks for discovery-rate recovery geometry.

Requirements: Python 3.10+, numpy, scipy.
Run: python recovery_geometry.py

All moment inputs below are exact-model quantities, not noisy observations.
The two-node bounds optimize over probability mixing measures; they do NOT
assert integer population realizability or confidence coverage.
Numerical tests corroborate the accompanying analytic arguments, not replace them.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import json
from typing import Callable

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import brentq, linprog

Array = NDArray[np.float64]


@dataclass(frozen=True)
class Measure:
    nodes: Array
    weights: Array

    def integral(self, function: Callable[[Array], Array]) -> float:
        return float(np.dot(self.weights, function(self.nodes)))

    def rate(self, time: float) -> float:
        return self.integral(lambda x: np.exp(-time * x))

    def hill(self, order: float) -> float:
        if abs(order - 1.0) < 1e-12:
            return float(np.exp(self.integral(lambda x: -np.log(x))))
        power_sum = self.integral(lambda x: x ** (order - 1.0))
        return power_sum ** (1.0 / (1.0 - order))


@dataclass(frozen=True)
class Envelopes:
    lower: Measure
    upper: Measure
    tilted_lower: Measure
    tilted_upper: Measure


def rate_envelopes(a: float, b: float, time: float,
                   rate: float, negative_derivative: float) -> Envelopes:
    """Sharp relaxed envelopes given mass=1, r(time), and -r'(time).

    The lower envelope is supported at {d,b}; the upper at {a,c}.
    Valid for 0<a<b<=1 and time>0. This double-precision implementation
    is intended for moderately scaled inputs, not near-singular data.
    """
    if not (0.0 < a < b <= 1.0 and time > 0.0 and rate > 0.0):
        raise ValueError("Require 0<a<b<=1, time>0, and rate>0.")
    if time * b > 500.0:
        raise ValueError("Use a high-precision/scaled solver for this time range.")
    beta = negative_derivative / rate
    tolerance = 2e-11
    if not (a - tolerance <= beta <= b + tolerance):
        raise ValueError("Tilted mean is outside the specified interval.")
    beta = float(np.clip(beta, a, b))
    target = 1.0 / rate
    low = np.exp(time * beta)
    high = ((b - beta) * np.exp(time * a)
            + (beta - a) * np.exp(time * b)) / (b - a)
    scale = max(1.0, target, high)
    if target < low - tolerance * scale or target > high + tolerance * scale:
        raise ValueError("The rate moments are infeasible in the relaxed class.")

    if abs(target - low) <= tolerance * scale:
        rho = Measure(np.array([beta]), np.array([1.0]))
        mu = Measure(np.array([beta]), np.array([1.0]))
        return Envelopes(mu, mu, rho, rho)
    if abs(target - high) <= tolerance * scale:
        nodes = np.array([a, b])
        weights = np.array([b - beta, beta - a]) / (b - a)
        rho = Measure(nodes, weights)
        mu = Measure(nodes, rate * np.exp(time * nodes) * weights)
        return Envelopes(mu, mu, rho, rho)

    def left_moment(c: float) -> float:
        return ((c - beta) * np.exp(time * a)
                + (beta - a) * np.exp(time * c)) / (c - a)

    def right_moment(d: float) -> float:
        return ((b - beta) * np.exp(time * d)
                + (beta - d) * np.exp(time * b)) / (b - d)

    c = brentq(lambda x: left_moment(x) - target, beta, b,
               xtol=5e-15, rtol=1e-14)
    d = brentq(lambda x: right_moment(x) - target, a, beta,
               xtol=5e-15, rtol=1e-14)
    rho_upper = Measure(np.array([a, c]),
                        np.array([c - beta, beta - a]) / (c - a))
    rho_lower = Measure(np.array([d, b]),
                        np.array([b - beta, beta - d]) / (b - d))

    def untilt(rho: Measure) -> Measure:
        return Measure(rho.nodes,
                       rate * np.exp(time * rho.nodes) * rho.weights)

    return Envelopes(untilt(rho_lower), untilt(rho_upper),
                     rho_lower, rho_upper)


def check_minimax_duality(delta: float) -> float:
    """Finite-grid primal modulus vs. affine minimax dual, infinity-norm noise."""
    x = np.linspace(0.1, 1.0, 31)
    n = len(x)
    K = np.column_stack((x, x*x))
    g = np.exp(-3*x)
    difference = np.column_stack((K.T, -K.T))
    equality = np.zeros((2, 2*n))
    equality[0, :n] = 1.0
    equality[1, n:] = 1.0
    primal = linprog(
        np.r_[-g, g],
        A_ub=np.vstack((difference, -difference)),
        b_ub=np.full(4, 2*delta), A_eq=equality, b_eq=np.ones(2),
        bounds=(0.0, None), method="highs")
    if not primal.success:
        raise RuntimeError(primal.message)
    # Variables: c_plus[2], c_minus[2], residual_upper, residual_lower.
    rows1 = np.column_stack((-K, K, -np.ones(n), np.zeros(n)))
    rows2 = np.column_stack((K, -K, np.zeros(n), np.ones(n)))
    dual = linprog(
        np.r_[np.full(4, delta), .5, -.5],
        A_ub=np.vstack((rows1, rows2)), b_ub=np.r_[-g, g],
        bounds=[(0., None)]*4 + [(None, None)]*2, method="highs")
    if not dual.success:
        raise RuntimeError(dual.message)
    gap = abs(-.5*primal.fun - dual.fun)
    assert gap < 2e-8, gap
    return float(gap)


def main() -> None:
    rng = np.random.default_rng(20261111)
    tests = {"moment_checks": 0, "rate_sandwiches": 0,
             "hill_sandwiches": 0, "lp_endpoint_checks": 0,
             "minimax_duality_checks": 0, "exact_rational_checks": 0}
    maximum_moment_error = 0.0
    maximum_lp_gap = 0.0

    for case in range(100):
        a, b = float(rng.uniform(.02, .15)), float(rng.uniform(.6, 1.0))
        time = float(rng.uniform(.5, 6.0))
        n = int(rng.integers(3, 9))
        original = Measure(rng.uniform(a + .05*(b-a), b - .05*(b-a), n),
                           rng.dirichlet(np.full(n, 2.0)))
        rate = original.rate(time)
        derivative = original.integral(lambda x: x*np.exp(-time*x))
        envelopes = rate_envelopes(a, b, time, rate, derivative)
        for measure in (envelopes.lower, envelopes.upper):
            errors = [abs(np.sum(measure.weights)-1),
                      abs(measure.rate(time)-rate),
                      abs(measure.integral(lambda x:x*np.exp(-time*x))-derivative)]
            maximum_moment_error = max(maximum_moment_error, *errors)
            assert max(errors) < 1e-9
            tests["moment_checks"] += 3
        for s in np.r_[0., time, np.geomspace(.001, 100., 45)]:
            lower = envelopes.lower.rate(float(s))
            value = original.rate(float(s))
            upper = envelopes.upper.rate(float(s))
            assert lower <= value + 2e-10 and value <= upper + 2e-10
            tests["rate_sandwiches"] += 1
        for q in np.linspace(0., 2., 17):
            lo, value, hi = (m.hill(float(q)) for m in
                             (envelopes.lower, original, envelopes.upper))
            assert lo <= value + 2e-8 and value <= hi + 2e-8
            tests["hill_sandwiches"] += 1

        if case < 12:
            grid = np.unique(np.r_[np.linspace(a, b, 201),
                                   envelopes.lower.nodes, envelopes.upper.nodes])
            A = np.vstack((np.ones_like(grid), np.exp(-time*grid),
                           grid*np.exp(-time*grid)))
            rhs = np.array([1., rate, derivative])
            for kernel in (lambda x: np.exp(-(time+2)*x),
                           lambda x: 1/x, lambda x: -np.log(x)):
                g = kernel(grid)
                for sign, measure in ((1., envelopes.lower),
                                      (-1., envelopes.upper)):
                    answer = linprog(sign*g, A_eq=A, b_eq=rhs,
                                     bounds=(0., None), method="highs",
                                     options={"primal_feasibility_tolerance":1e-9,
                                              "dual_feasibility_tolerance":1e-9})
                    if not answer.success:
                        raise RuntimeError(answer.message)
                    gap = abs(sign*answer.fun - measure.integral(kernel))
                    maximum_lp_gap = max(maximum_lp_gap, gap)
                    assert gap < 5e-7, gap
                    tests["lp_endpoint_checks"] += 1

    duality_gaps = []
    for delta in (0., .001, .01, .03, .1):
        duality_gaps.append(check_minimax_duality(delta))
        tests["minimax_duality_checks"] += 1

    p_exact = [Fraction(1,2), Fraction(1,4), Fraction(1,4)]
    assert sum(p_exact) == 1
    assert sum(x*x for x in p_exact) == Fraction(3,8)
    assert 1 / Fraction(3,8) == Fraction(8,3)
    assert 6-8*Fraction(3,8) == 3
    tests["exact_rational_checks"] = 4

    p = np.array([.5, .3, .2])
    original = Measure(p, p)
    time = 3.
    rate = original.rate(time)
    derivative = original.integral(lambda x:x*np.exp(-time*x))
    env = rate_envelopes(.05, .6, time, rate, derivative)
    kernels = {
        "future_discoveries_u2": lambda x: np.exp(-time*x)*(-np.expm1(-2*x))/x,
        "unseen_richness": lambda x: np.exp(-time*x)/x,
        "support_size": lambda x: 1/x,
        "entropy": lambda x: -np.log(x)}
    example = {name: [m.integral(kernel) for m in (env.lower, original, env.upper)]
               for name, kernel in kernels.items()}
    example["hill_2"] = [m.hill(2.) for m in (env.lower, original, env.upper)]
    example["chao_lower"] = rate*rate/derivative
    example["rate"] = rate
    example["negative_derivative"] = derivative
    example["lower_nodes"] = env.lower.nodes.tolist()
    example["lower_weights"] = env.lower.weights.tolist()
    example["upper_nodes"] = env.upper.nodes.tolist()
    example["upper_weights"] = env.upper.weights.tolist()
    print(json.dumps({"tests": tests, "total_checks": sum(tests.values()),
                      "maximum_moment_error": float(maximum_moment_error),
                      "maximum_lp_gap": float(maximum_lp_gap),
                      "maximum_minimax_duality_gap": max(duality_gaps),
                      "example_lower_true_upper": example}, indent=2))


if __name__ == "__main__":
    main()
