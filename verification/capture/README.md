# Capture–recapture verification

Run `make check-capture` from the repository root. The suite uses exact fractions and symbolic generating-function derivatives; it writes [results.json](results.json).

It checks the formulas against independently enumerated labelled subsets, marked positions, repeated-capture histories, and independent Bernoulli histories. The checks cover reciprocal moments, Chapman and Bailey bias, signed transfers, confidence coverage and endpoint witnesses, sequential variance, repeated-estimator distributions, and martingale crossing probabilities. The analytic proofs are in [the capture–recapture manuscript](../../manuscripts/capture_recapture.tex).

## Using the exact routines

With the repository's `src` directory on `PYTHONPATH`:

```python
from partial_population.capture_recapture import (
    chapman, fixed_interval, sequential_estimate, sequential_interval,
)

# Mark 20, then sample 25 without replacement, observing 5 recaptures.
assert chapman(20, 25, 5) == 90
assert fixed_interval(20, 25, 5) == (56, 274)

# A different experiment: stop without replacement at the fifth recapture.
# Here 25 is the resulting random sample size.
assert sequential_estimate(20, 5, 25) == 104
assert sequential_interval(20, 5, 25) == (61, 274)
```

Intervals default to at least 95% coverage under the specified design. Supply an exact `Fraction` to change the error level. `None` is an infinite upper endpoint; it occurs at zero overlap in the fixed-size design. The endpoint search uses monotone bracketing over unbounded integer population sizes, not a numerical upper cap.

`bailey` and `binomial_interval` use a second sample **with replacement**. `inverse_binomial_estimate` and `inverse_binomial_interval` stop at a predetermined number of marked draws with replacement. Repeated appearances count in these designs.

`capture_evalue_path` supplies candidate-wise calibration for repeated capture occasions, with all newly encountered individuals marked and returned. `draw_evalue_path` supplies calibration for sequential binary draws without replacement. Betting factors must be chosen predictably, before seeing the observation they score. Running maximum inversion gives simultaneous coverage; a fixed-recapture-count interval does not itself justify a data-selected stopping count.

`repeated_distribution` computes the exact law of the Schnabel or Schumacher–Eschmeyer ratio by dynamic programming. It retains any atom at infinity. `history_transform` conditions on a fixed finite catchability population. `horvitz_thompson` requires known positive inclusion probabilities. `chao_incidence` is a point estimate with an explicit zero-doubleton convention; it is not a confidence bound.

These are exact research routines for moderate designs. They do not implement a production continuum optimizer for heterogeneous population measures or guarantee shortest two-sided intervals.
