# Certified computation and classical-data audit

Date: 2026-11-17. Baseline: repository commit `35f6902`.

## Implemented scope

- Joint-region outer endpoints for realized missing mass and realized future
  discoveries under fixed-size iid sampling. Eight predeclared directions use
  fixed Poisson efforts n and m and explicit conditioning penalties.
- Complete [0,1] kernel coverage with FLINT/Arb, including the zero limit;
  rational entropy-dual certificates; checked outward rounding; feasible
  relaxed-measure witnesses and endpoint gap reports.
- Separate estimation and calibration choices: fixed-count and Poisson
  Good–Turing normal intervals, generic normal calibration with declared
  assumptions and a bias allowance, and custom predeclared linear scores.
- Known-removal capture kernels, historical ratio reproduction, initial-N
  score roots, and an exact rational two-tail martingale interval. The sunfish
  example includes all three known deaths and does not interpret a final
  placeholder as additional removal.
- Five abundance datasets, trout, two full-history datasets, all fifteen
  textbook frequency series, and explicit source-version sensitivity for
  taxicab A, taxicab B, and golf tees. Source data are hashed and attributed.
- A paired study with 1,600 initial/future sample pairs, 3,200 target cases,
  64 summaries, 1,438 distinct simulation certificates, 12 real-data
  certificates, and 16 shared continuum tables.

## Mathematical checks

The empirical manuscript proves the entropy-dual endpoint inequality and
explains fixed-count conditioning. Optimizer success is neither necessary nor
sufficient for acceptance: replay uses the stored rational weights and entire
interval cells. A failed optimizer uses a certified individual direction;
failed kernel enclosure raises an error. Empty output intersections use full
support. A missing primal witness means the optimization gap is unknown,
not that coverage has failed. Witnesses concern the relaxed measure problem,
not integer realizability.

The normal-approximation derivation retains the fixed-count variance correction
and exact finite-n bias. It states the need for a suitable CLT, nondegenerate
and consistent studentization, and negligible standardized bias. Gaussian
conditioning motivates a formula but is not itself a finite-count proof.
The generic calibration interface cannot certify an unsupported user assumption.

Known-removal products have conditional mean one under the declared uniform
remaining-population kernel. Constant tilts make their values monotone in N0;
integer doubling/bisection and the infinite-N limit avoid an arbitrary cutoff.
Exact finite-law checks corroborate the construction. Nonlinear score roots
are not mislabeled as unbiased estimates.

## Findings and limits

The certified simulation intervals cover all 200 targets in each scenario;
they are substantially wider than the normal choices. Joint versus split
improvement is modest. The 18 failed simulation endpoint proposals remain in
the results through certified fallback. Two real-data lower endpoints also
use that fallback. Replay reports 2,880 verified witnesses out of 2,900
endpoints. These data do not establish optimality or state-of-the-art dominance.

The manuscript preserves truncated butterfly and Shakespeare tables, the BCI
diameter/census distinction, the unresolved cottontail Schnabel discrepancy,
the wood-mouse missing columns, and disputed singleton frequencies. It does
not claim that ecological or linguistic observations have been validated as
iid samples. Normal calibration, exact calibration, observation-model choice,
and uncertainty about source records are separate axes of analysis.

No authenticated raw wartime serial sample with a selection rule was recovered.
An archival locator and aggregate production comparisons do not supply that
sampling law. A known coarsening can be composed with the sampling kernel;
unknown records or mechanisms can be retained as nuisance uncertainty.

Effort-optimized regions, integer-population optimization, dependent-text
holdouts, all modern prediction benchmarks, and broad open-population capture
models remain extensions. This is research code using trusted FLINT/Arb and
Python, not a proof-assistant formalization.

## Verification and provenance

`make check` runs all original suites plus finite-law/tamper/failure tests,
removal/normal/custom-score tests, and replay of every stored empirical
certificate, simulated target, interval and summary. `make pdf` rejects errors,
unresolved references and overfull/underfull boxes. The PDF manifest binds
source, generated tables/figure, analysis JSON, and published PDF hashes.
Visual page review complements log checks.

The unavailable earlier Work archive was not recovered. This implementation
was reconstructed from the published theory and preserved sources. The new
certificate counts and results are its own; they are not attributed to the
unavailable implementation. Binary files are transferred by ordinary Git,
not embedded as base64 in conversation/tool arguments.
