#!/usr/bin/env python3
"""Numerical corroboration of the full Good--Toulmin prediction variance.

This finite-series comparison is supplemental, not a proof or a rigorous
truncation-error certificate. The exact identity is proved in the manuscript.
"""
from pathlib import Path
import json
import numpy as np
from scipy.stats import poisson


def main() -> None:
    counts = np.arange(1, 201, dtype=float)
    populations = (np.array([0.5, 0.3, 0.2]), np.array([1.0]), np.full(10, 0.1))
    cases = 0
    max_error = 0.0
    for p in populations:
        for t in (0.5, 3.0, 10.0):
            mass = poisson.pmf(counts[None, :], t * p[:, None])
            for ratio in (0.25, 1.0, 2.0):
                # Predictor and future-discovery contributions are disjoint
                # for each species; the resulting per-species error has mean zero.
                second_moment = np.sum(mass * ratio ** (2.0 * counts))
                future_mean = np.sum(np.exp(-t * p) * (-np.expm1(-ratio * t * p)))
                series_variance = second_moment + future_mean
                closed_variance = np.sum(np.exp(-t * p) *
                    (np.exp(ratio * ratio * t * p) - np.exp(-ratio * t * p)))
                error = abs(series_variance - closed_variance) / max(1.0, abs(closed_variance))
                if not np.isfinite(error) or error >= 1e-11:
                    raise AssertionError((p.tolist(), t, ratio, error))
                max_error = max(max_error, float(error))
                cases += 1
    result = {
        'numerical_variance_checks': cases,
        'max_scaled_error': max_error,
        'note': 'Finite series through r=200 compared with the analytic generating-function identity. Numerical corroboration, not an error-bound proof.'
    }
    Path(__file__).with_suffix('.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
