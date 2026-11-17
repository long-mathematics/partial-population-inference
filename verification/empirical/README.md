# Paired empirical study

`classical_results.json` holds the data transformations and model-based real-data
results. `coverage_width.json` contains 64 coverage/width summary rows from
200 paired replicates in each of eight population/sample-size scenarios.
Initial and independent future counts are multinomial. The four populations
are uniform on 100 species, a normalized geometric sequence on 100 species,
a mass-0.8 common species with 200 mass-0.001 species, and the normalized
Rothamsted observed abundances. The last is a specified reference population,
not an estimate of real-world coverage.

All methods see the same draws. Targets are realized missing mass and the
realized number of new species in a future sample of size m=n. Intervals have
nominal 95% coverage. Benchmarks are equal-budget single-direction calibration,
full target support, and (for missing mass) explicitly asymptotic Good–Turing
normal intervals with fixed-count and Poisson variance choices. No total-richness or
point-prediction benchmark is mislabeled as a matching interval benchmark.

Compressed files preserve every unique certificate, shared continuum table,
and all 3,200 target cases. `check_artifacts.py` checks all source hashes,
recomputes the complete continuum enclosures and endpoint certificates without
running the optimizer, regenerates the random draws from seed 20261117, and
checks all intervals, coverage counts, exact binomial Monte Carlo intervals,
mean widths, and their Monte Carlo standard errors.

Floating-point simulation truths and descriptive summaries are compared with
relative and absolute tolerances of 1e-12 for portability across CPU builds.
Integer counts and coverage decisions are checked exactly; the separate
continuum-certificate verifier retains its outward interval checks.

It also reproduces the normal intervals for all abundance datasets, the
known-removal sunfish estimates and interval, and the point-summary sensitivity
to disputed taxi and golf-tee source versions. The normal variants illustrate
calibration choices within the framework, not a mandatory conservative rule.

```sh
python scripts/extract_schnabel.py
PYTHONPATH=src python scripts/classical_analysis.py
PYTHONPATH=src python scripts/coverage_width.py
PYTHONPATH=src python verification/empirical/check_artifacts.py
python scripts/empirical_tables.py
```

The simulation regeneration can take several minutes. Routine CI replays the
stored certificates and cases instead. Conic proposals can differ across solver
versions; their validity is established by replay. Gzip container timestamps
are not part of the mathematical certificate.
