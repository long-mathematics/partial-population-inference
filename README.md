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
| **Lincoln–Petersen** | Population size from two capture samples and their overlap | The tag–retag hypergeometric experiment is stated; explicit treatment of this estimator remains to be incorporated. |
| **Chapman** | Bias-corrected two-sample capture–recapture estimation | A planned specialization of the same tag–retag experiment; its derivation and calibration are not yet developed here. |

The current theory also treats collision-based power sums, Hill diversity, and multiple-survey discovery patterns. Tag–retag (capture–recapture) and inverse coupon collection are explicit continuation targets in the [research roadmap](research/ROADMAP.md).

## Mathematical framework

- **Sampling and transfer.** Positive kernels transport decisions and calibration; dominated kernels transport calibration with a loss of power; signed kernels transport expectations and reveal extrapolation instability.
- **Recovery geometry.** The abundance measure turns discovery curves and other population targets into linear functionals. Moment constraints yield simultaneous envelopes, while integer realizability distinguishes actual populations from relaxed measures.
- **Calibration and computation.** Independent Poisson species counts give an additive cumulant identity. Joint e-values calibrate observed, realized, and future quantities; dominated deletion carries calibration to fixed-size samples. Polynomial algorithms compute the deletion averages, and common-score conditioning permits deterministic optimization over Poisson effort.

The assumptions depend on the experiment and target. Some richness and envelope results require a known positive lower abundance bound; missing-mass and finite-horizon calibration do not always require one. The German tank model is a separate finite experiment within the framework.

## Read the papers

| Manuscript | Main contents |
|---|---|
| [Discovery-rate inversion](manuscripts/discovery_rate_inversion.tex) | Sampling identities; moment recovery and simultaneous envelopes; Hill-order anomaly; integer realizability; joint calibration; exact tank inference. |
| [Transfer and generation](manuscripts/transfer_generation.tex) | Positive, signed, and dominated transfers; random-effort deficiency; fixed-size calibration and prediction; batch extraction. |
| [Additive cumulants and fixed-sample computation](research/master_identity_and_deletion.tex) | Master identity with integrability conditions; polynomial deletion algorithms; threshold costs; effort-optimized full-count calibration. |

The manuscript bibliographies contain the references and comparisons with existing literature. See also the [literature audit](audit/LITERATURE.md), [mathematical and computational audit](audit/AUDIT.md), and [source provenance](provenance/SOURCES.md).

## Reproduce

Use Python 3.12 or later:

```sh
python -m pip install -r requirements.txt
make check
```

To run an individual verification suite, use `make check-discovery`, `make check-transfer`, or `make check-new`. The suites include exact rational certificates, symbolic identities, and numerical checks.

To build all three PDFs, install a TeX distribution with `latexmk` and the packages declared in the sources, then run:

```sh
make pdf
```

PDFs are written to `build/`. The [GitHub Actions workflow](https://github.com/long-mathematics/partial-population-inference/actions/workflows/verify.yml) runs the verification suites and publishes compiled PDFs as build artifacts.

[src/partial_population/deletion.py](src/partial_population/deletion.py) implements log-domain population and prediction deletion scores and common-score effort optimization. The population polynomial takes O(n²) arithmetic operations and O(n) storage. This is research code; certified projection and numerical error bounds remain [development goals](research/ROADMAP.md).

## Author and citation

Christopher D. Long. Please cite the relevant manuscript and its version; repository citation metadata is in [CITATION.cff](CITATION.cff). Original code and text are released under the [MIT license](LICENSE).
