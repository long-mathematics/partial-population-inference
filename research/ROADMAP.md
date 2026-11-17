# Research continuation

The central objective remains a transfer-and-generation theory: state which experiment and target are given, generate calibrated statistics from the Poisson additive transform where applicable, and identify what survives transfer.

## 1. Certified projected intervals

The [classical-data manuscript](../manuscripts/classical_data.tex) and
[projection implementation](../src/partial_population/projection.py) now
provide continuum-certified outer intervals for a specified eight-direction,
fixed-effort region. Entropy duality and interval arithmetic certify endpoints;
feasible relaxed measures provide explicit optimization-gap bounds. Saved
certificates replay independently of the conic optimizer. Integer-realizable
population optimization remains open.

Evaluate a preassigned mixture of effort-optimized full-count and fixed-effort deletion calibration. The continuation theorem optimizes normalizers for a common score; arbitrary data-selected directions or maxima of full-deletion e-values require separate justification.

Completed first deliverable: reproducible missing-mass and fixed-horizon
novelty intervals, finite-law checks, saved continuum certificates and
optimization-gap reports where witnesses are available. Improve efficiency
without confusing a stronger theoretical region with the region actually
implemented. Normal approximations and custom scores are selectable choices
with their own assumptions, rather than being excluded by the framework.

## 2. Matched comparisons and classical data

The first paired study and multi-dataset manuscript are now implemented:
four specified populations, two sample sizes, 200 replicates each, 64
coverage/width summaries, and five abundance datasets plus classical capture
series. The [audit](../audit/EMPIRICAL.md) records assumptions and source
discrepancies. Known sunfish removals and disputed taxi/golf-tee frequencies
are handled explicitly. The historical serial-number record still lacks an
authenticated raw sample and sampling rule; no substitute is invented.

Compare richness estimates on a common minimum-mass class with Chao and Wu–Yang as appropriate. Compare prediction with smoothed Good–Toulmin and Eriksson under their stated regimes. Treat individual unseen probability and maximum unseen probability as distinct targets when discussing Painsky and Colombi et al. Do not compare point-estimator error bounds with confidence-interval width without constructing a matching inferential task.

## 3. Capture–recapture extensions and inverse coupon

The [capture–recapture paper](../manuscripts/capture_recapture.tex) now develops the closed-population uniform fixed-size, with-replacement, repeated-capture, and preassigned-recapture experiments. It includes the classical estimators, exact bias and recovery identities, signed-transfer classification, exact confidence limits, and separate martingales for arbitrary stopping. Independent heterogeneous histories connect to the counting-measure cumulant and moment framework. The [audit](../audit/CAPTURE_RECAPTURE.md) records the scope.

Next, seek a sharp bounded-heterogeneity inverse-coupon interval with a finite-parameter or moment dual certificate. Detection heterogeneity, tag loss, erroneous matches, spatial detection, and open populations each require explicit observation models. General capture-history calibration does not itself produce a certified continuum optimizer or a new sharp interval for all of those models.

The older private tag–retag work has not been recovered here; the current development is self-contained and does not claim to reproduce it.

## 4. Structural limits and manuscript consolidation

Preserve the established theorem numbering while integrating specializations through explicit cross-references. Put detailed prior-art comparisons near the end of each paper, retaining citations at use. A combined exposition should proceed from experiments and transfer kernels to the additive transform, measure geometry, calibration, and dual certificates. The tank game belongs as a finite decision example; it is not a Poisson species-factorization corollary.

Connect bounded-batch fluctuation work after recovering its current version. Overlap changes the needed joint incidence data, and one-label extraction loses information about the original batches.

Before submission, replace literature-search impressions by explicit theorem-level comparisons. The moment-envelope and common-score conditioning arguments may have close antecedents; no priority conclusion follows from the searches completed so far.
