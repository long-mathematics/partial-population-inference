"""Independent finite-enumeration checks for the new research note."""
import itertools
import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path

from partial_population.deletion import (
    fingerprint_score, log_deletion_score, log_lost_species_coefficients,
    log_poisson_mass, log_prediction_deletion_score, log_subset_coefficients,
    optimized_full_count_log_e,
)

HERE = Path(__file__).resolve().parent
checks = Counter()
max_error = 0.0


def close(a, b, group, tol=2e-10):
    global max_error
    error = abs(a - b) / max(1.0, abs(a), abs(b))
    max_error = max(max_error, error)
    assert error < tol, (group, a, b, error)
    checks[group] += 1


def partitions(n, smallest=1):
    if n == 0:
        yield ()
    for x in range(smallest, n + 1):
        for rest in partitions(n - x, x):
            yield (x,) + rest


def brute_coeff(counts, weights):
    labels = [i for i, x in enumerate(counts) for _ in range(x)]
    ans = [[Fraction(0) for _ in range(len(counts) + 1)]
           for _ in range(len(labels) + 1)]
    for mask in range(1 << len(labels)):
        c = Counter(labels[j] for j in range(len(labels)) if mask >> j & 1)
        score = math.prod(weights.get(v, Fraction(1)) for v in c.values())
        ans[sum(c.values())][len(counts) - len(c)] += score
    return ans


def poisson_cgf(p, direction, t, eta=0.0, u=0.0, target="novelty"):
    """Independent finite-population product, not numerical integration."""
    ans = 0.0
    for x in p:
        empty = math.exp(-t * x)
        val = 1.0 + sum(empty * (t * x) ** r / math.factorial(r)
                        * math.expm1(lam)
                        for r, lam in enumerate(direction, 1))
        if target == "novelty":
            val += empty * (-math.expm1(-u * x)) * math.expm1(eta)
        elif target == "missing_mass":
            val += empty * math.expm1(eta * x)
        else:
            raise ValueError(target)
        ans += math.log(val)
    return ans


def outcomes(p, n):
    """Enumerate ordered draws, accumulating equal count vectors."""
    result = Counter()
    for word in itertools.product(range(len(p)), repeat=n):
        c = tuple(word.count(i) for i in range(len(p)))
        result[c] += math.prod(p[i] for i in word)
    return result.items()


def main():
    weights = {1: Fraction(2), 2: Fraction(1, 3), 3: Fraction(3, 2)}
    direction = tuple(math.log(float(weights[r])) for r in range(1, 4))
    for n in range(8):
        for counts in partitions(n):
            exact = brute_coeff(counts, weights)
            coefficients = log_subset_coefficients(counts, direction)
            joint = log_lost_species_coefficients(counts, direction)
            for k, row in enumerate(exact):
                close(math.exp(coefficients[k]), float(sum(row)),
                      "univariate_vs_exact_subsets")
                for r, val in enumerate(row):
                    close(math.exp(joint[k][r]), float(val),
                          "bivariate_vs_exact_subsets")
            for t in (0.7, 2.4):
                expected = sum(math.exp(log_poisson_mass(k, t))
                               * float(sum(row)) / math.comb(n, k)
                               for k, row in enumerate(exact))
                close(math.exp(log_deletion_score(counts, direction, t)),
                      expected, "population_deletion_vs_exact_subsets")
                for eta in (-0.8, 0.0, 0.7):
                    m, u = 3, 1.6
                    expected = sum(
                        math.exp(log_poisson_mass(k, t)
                                 + log_poisson_mass(ell, u)
                                 - max(eta, 0) * (m - ell)
                                 - max(-eta, 0) * min(r, ell))
                        * float(val) / math.comb(n, k)
                        for k, row in enumerate(exact)
                        for r, val in enumerate(row)
                        for ell in range(m + 1))
                    got = log_prediction_deletion_score(
                        counts, direction, eta, m, t, u)
                    close(math.exp(got), expected,
                          "prediction_deletion_vs_exact_subsets")
                    full = (log_poisson_mass(n, t)
                            + log_poisson_mass(m, u)
                            + fingerprint_score(counts, direction))
                    assert got >= full - 1e-12
                    checks["full_count_domination"] += 1

    populations = [(Fraction(1, 2), Fraction(1, 3), Fraction(1, 6)),
                   (Fraction(1, 4), Fraction(3, 4))]
    grid = tuple(itertools.product((0.4, 1.2, 2.8, 5.0), repeat=2))
    worst_e = 0.0
    for rational_p in populations:
        p = tuple(map(float, rational_p))
        for n, m in ((1, 1), (3, 2), (4, 3)):
            for eta in (-1.1, 0.0, 0.9):
                for target in ("novelty", "missing_mass"):
                    total_opt, total_delete = 0.0, 0.0
                    for counts, prob in outcomes(rational_p, n):
                        observed = tuple(c for c in counts if c)
                        score = fingerprint_score(observed, direction)
                        if target == "missing_mass":
                            y = sum(p[i] for i, c in enumerate(counts) if c == 0)
                            choices = [(y, Fraction(1))]
                            sizes, efforts_grid = (n,), [(t,) for t in (0.4, 1.2, 2.8, 5.0)]
                            def cgf(e):
                                return poisson_cgf(p, direction, e[0], eta,
                                                   target=target)
                        else:
                            choices = [(sum(c == 0 and z > 0
                                            for c, z in zip(counts, future)), pr)
                                       for future, pr in outcomes(rational_p, m)]
                            sizes, efforts_grid = (n, m), grid
                            def cgf(e):
                                return poisson_cgf(p, direction, e[0], eta, e[1])
                        for y, pr in choices:
                            log_e, _ = optimized_full_count_log_e(
                                score + eta * y, sizes, efforts_grid, cgf)
                            total_opt += float(prob * pr) * math.exp(log_e)
                            if target == "novelty":
                                log_e_delete = log_prediction_deletion_score(
                                    observed, direction, eta, m, 1.2, 2.8)
                                total_delete += float(prob * pr) * math.exp(
                                    log_e_delete + eta * y - cgf((1.2, 2.8)))
                    assert total_opt <= 1 + 2e-12, (target, total_opt)
                    checks["optimized_fixed_count_e_expectation"] += 1
                    worst_e = max(worst_e, total_opt)
                    if target == "novelty":
                        assert total_delete <= 1 + 2e-12, total_delete
                        checks["joint_deletion_e_expectation"] += 1

    # Population relaxation: convexity survives effort minimization.
    # Here measures are mixtures of two actual abundance measures.
    p, q = (0.5, 0.3, 0.2), (0.7, 0.3)
    score, n = 0.7, 4
    one_grid = [(t,) for t in (0.4, 1.2, 2.8, 4.0, 5.0)]
    endpoints = [optimized_full_count_log_e(
        score, (n,), one_grid,
        lambda e, pp=pp: poisson_cgf(pp, direction, e[0]))[0]
        for pp in (p, q)]
    for w in (0.1, 0.3, 0.5, 0.9):
        val, _ = optimized_full_count_log_e(
            score, (n,), one_grid,
            lambda e: w * poisson_cgf(p, direction, e[0])
            + (1 - w) * poisson_cgf(q, direction, e[0]))
        assert val <= w * endpoints[0] + (1 - w) * endpoints[1] + 1e-12
        checks["log_convexity_on_measure_mixtures"] += 1

    # Search for an explicit improvement over the conventional effort t=n.
    p, n = (0.5, 0.3, 0.2), 4
    example_direction = (2.0, -1.0)
    effort_grid = [(j / 20,) for j in range(1, 201)]
    cgf = lambda e: poisson_cgf(p, example_direction, e[0])
    optimized, selected = optimized_full_count_log_e(
        0.0, (n,), effort_grid, cgf)
    baseline = -cgf((float(n),)) + log_poisson_mass(n, float(n))
    assert optimized >= baseline - 1e-12
    result = {
        "all_passed": True, "checks": dict(checks),
        "total_checks": sum(checks.values()),
        "max_relative_error": max_error,
        "max_tested_optimized_e_expectation": worst_e,
        "effort_optimization_example": {
            "population": p, "n": n, "direction": example_direction,
            "grid": "0.05, 0.10, ..., 10.00",
            "selected_effort": selected[0],
            "log_e_gain_over_t_equals_n": optimized - baseline,
        },
        "interpretation": "Finite numerical corroboration; no global optimality or floating-point certification.",
    }
    (HERE / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
