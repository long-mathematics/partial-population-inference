# Certified projection

The implementation projects a fixed eight-direction joint e-value region over
all probability measures on [0,1]. It uses fixed Poisson efforts t=n and u=m,
then conditions on the specified sample sizes. It does not optimize effort or
impose integer realizability of the relaxed population measure.

CLARABEL/CVXPY proposes endpoint weights. FLINT/Arb verifies the kernel on every
cell of a partition covering the entire continuum, including the removable
singularity at zero. An entropy dual inequality bounds each endpoint. Rational
dual weights sum exactly to 10^12. Outward binary64 rounding is checked before
serialization and again on replay. A separately checked primal witness, when
available, brackets the optimum of the relaxed region. Its gap is not a bound
on the difference from the integer-realizable problem.

Install `requirements-empirical.txt`, then run:

```sh
PYTHONPATH=src python scripts/project_interval.py --n 1000 --fingerprint 50 --output build/example
PYTHONPATH=src python scripts/project_interval.py --verify build/example
PYTHONPATH=src python scripts/project_interval.py --n 1000 --future 1000 --fingerprint 50 20 10 --output build/future
make check-projection
make check-empirical
PYTHONPATH=src python scripts/project_interval.py --n 1000 --fingerprint 50 20 --normal fixed_n --output build/normal
PYTHONPATH=src python scripts/project_interval.py --n 1000 --fingerprint 50 20 --normal poisson --output build/poisson-normal
```

Normal calibration is a separate, selectable choice. Both commands reproduce
the Good–Turing center and their stated studentizers exactly, preserve the
untruncated interval, and record approximation assumptions. They do not claim
finite-sample certification. `normal_interval` accepts other point estimates,
variances and bias allowances with a declared target and assumptions.
`--score-weights` supplies a predeclared linear score to the certified path;
data-dependent fitting requires separate justification.

`removals.py` implements the hypergeometric kernel with known removals of
captured animals. It updates both available population and marked pool,
provides historical ratios and removal-aware score roots, and exactly inverts
two rational test martingales for initial population size. `check_choices.py`
enumerates its sampling laws, reproduces normal choices, and verifies a
nondefault score through the continuum certificate path.

Missing mass uses the known n and Phi_1; realized future discoveries use Phi_1,
Phi_2, Phi_3 and a future sample size chosen independently of the data. Coverage
is unconditional over the initial and future iid samples, not conditional on
the observed fingerprint. Custom directions must be chosen before seeing the
sample. Multiple horizons here have marginal 95% prediction coverage each,
not simultaneous 95% coverage.

Solver failure uses a certified one-direction bound. Empty reported intersections
use full target support. Failed kernel enclosure raises an error instead of
returning an uncertified interval. A caller supplying its own kernel table must
run `verify_kernel_table` before trusting it; `verify_interval` assumes that
shared-table check has already passed.

`check_projection.py` enumerates exact rational fixed-size laws, tests zero
limits, forces optimizer failure, and rejects forged certificates. General
coverage follows from the analytic proof in the empirical manuscript, not from
these finite examples or Monte Carlo. Replay trusts the kernel implementation,
FLINT/Arb and Python arithmetic; it is not proof-assistant verification.
