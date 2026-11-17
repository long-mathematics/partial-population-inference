# Partial-Population Inference

Sampling kernels, recovery geometry, and finite-sample calibration for partially observed populations.

This project studies missing mass, species discovery, population size, and inference across sampling designs. Its aim is a common mathematical framework for familiar estimators and new confidence procedures: identify the observation kernel, determine which population quantities it recovers, and track what remains valid when the experiment changes.

**Status:** active research drafts with reproducible computational checks.

## Classical estimators and applications

| Estimator or method | Inferential target | Treatment in this repository |
|---|---|---|
| **Good–Turing** | Missing mass and discovery rate | Singleton identities, the fixed-sample size shift, and the exact Poissonized relation; joint calibration for realized missing mass. |
| **Good–Toulmin** | New species in a future sample | Discovery-curve extrapolation, prediction-error moments, and long-horizon instability and repairs. |
| **Chao-type richness bounds** | Unseen species and total richness | The underlying population moment inequality and sharper moment-based recovery bounds. The observed-count ratio alone is not a confidence bound. |
| **German tank estimator** | Population size from sampled serial numbers | The classical minimum-variance unbiased estimator, general triangular inversion, exact upper confidence limits, and a finite decision game. |
| **Lincoln–Petersen** | Population size from two capture samples and their overlap | Exact overlap moments, reciprocal-factorial recovery, integer likelihood maximizers, and exact confidence limits. |
| **Chapman** | Bias-corrected two-sample capture–recapture estimation | Exact reciprocal-shift derivation, boundary bias, finite-design moments, and separately calibrated confidence intervals. |
| **Bailey** | Population size with replacement in the second sample | Binomial correction, exact bias, and exact confidence intervals. |
| **Schnabel; Schumacher–Eschmeyer** | Population size across repeated capture occasions | Conditional estimating equations, exact finite-design distributions, and separate test-martingale calibration. |
| **Horvitz–Thompson** | Population totals with known inclusion probabilities | Unbiased weighting, variance, and its relation to capture-history models. |
| **Inverse-sampling estimators** | Population size when sampling stops at a prescribed recapture count | With- and without-replacement laws, unbiased estimates, exact intervals, and the German tank equivalence. |

The current theory also treats collision-based power sums, Hill diversity, and multiple-survey discovery patterns. Closed-population tag–retag (capture–recapture) is developed in the [capture–recapture paper](manuscripts/capture_recapture.tex). Broader heterogeneity and inverse coupon collection remain in the [research roadmap](research/ROADMAP.md).

## Mathematical framework

- **Sampling and transfer.** Positive kernels transport decisions and calibration; dominated kernels transport calibration with a loss of power; signed kernels transport expectations and reveal extrapolation instability.
- **Recovery geometry.** The abundance measure turns discovery curves and other population targets into linear functionals. Moment constraints yield simultaneous envelopes, while integer realizability distinguishes actual populations from relaxed measures.
- **Calibration and computation.** Independent Poisson species counts give an additive cumulant identity. Joint e-values calibrate observed, realized, and future quantities; dominated deletion carries calibration to fixed-size samples. Polynomial algorithms compute the deletion averages, and common-score conditioning permits deterministic optimization over Poisson effort.

The assumptions depend on the experiment and target. Some richness and envelope results require a known positive lower abundance bound; missing-mass and finite-horizon calibration do not always require one. The German tank model is a separate finite experiment within the framework.

## Read the papers

| Manuscript | Download | Main contents |
|---|---|---|
| [Discovery-rate inversion](manuscripts/discovery_rate_inversion.tex) | [PDF](https://github.com/long-mathematics/partial-population-inference/raw/refs/heads/main/pdf/discovery_rate_inversion.pdf) | Sampling identities; moment recovery and simultaneous envelopes; Hill-order anomaly; integer realizability; joint calibration; tank inference. |
| [Transfer and generation](manuscripts/transfer_generation.tex) | [PDF](https://github.com/long-mathematics/partial-population-inference/raw/refs/heads/main/pdf/transfer_generation.pdf) | Positive, signed, and dominated transfers; random-effort deficiency; fixed-size calibration and prediction; batch extraction. |
| [Additive cumulants and fixed-sample computation](research/master_identity_and_deletion.tex) | [PDF](https://github.com/long-mathematics/partial-population-inference/raw/refs/heads/main/pdf/master_identity_and_deletion.pdf) | Integrability conditions; polynomial deletion algorithms; threshold costs; effort-optimized full-count calibration. |
| [Capture–recapture](manuscripts/capture_recapture.tex) | [PDF](https://github.com/long-mathematics/partial-population-inference/raw/refs/heads/main/pdf/capture_recapture.pdf) | Classical estimators; exact recovery and transfers; fixed-size and sequential intervals; repeated captures; independent-history calibration. |

Each paper develops the mathematics first and places its detailed comparison with prior work near the end, with brief citations at the points of use. See also the [literature audit](audit/LITERATURE.md), [mathematical and computational audit](audit/AUDIT.md), and [source provenance](provenance/SOURCES.md).

## Reproduce

Use Python 3.12 or later:

```sh
python -m pip install -r requirements.txt
make check
```

To run an individual verification suite, use `make check-discovery`, `make check-transfer`, `make check-new`, or `make check-capture`. The suites include exact rational certificates, symbolic identities, and numerical checks.

To build all four PDFs, install a TeX distribution with `latexmk` and the packages declared in the sources, then run:

```sh
make pdf
```

PDFs are written to `build/`, and the build checks the logs for errors, unresolved references, and layout warnings. Run `make publish-pdfs` to refresh the downloadable copies in [pdf/](pdf/README.md). CI checks their source hashes and independently rebuilds all four papers; it also retains PDFs as build artifacts.

[src/partial_population/deletion.py](src/partial_population/deletion.py) implements log-domain population and prediction deletion scores and common-score effort optimization. The population polynomial takes O(n²) arithmetic operations and O(n) storage. This is research code; certified projection and numerical error bounds remain [development goals](research/ROADMAP.md).

[src/partial_population/capture_recapture.py](src/partial_population/capture_recapture.py) supplies exact rational estimators, tail inversion with unbounded population search, signed transfers, and sequential e-values. See the [capture verification guide](verification/capture/README.md) and [audit](audit/CAPTURE_RECAPTURE.md) for the covered designs and limitations.

## Author and citation

Christopher D. Long. Please cite the relevant manuscript and its version; repository citation metadata is in [CITATION.cff](CITATION.cff). Original code and text are released under the [MIT license](LICENSE).
