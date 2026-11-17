# Partial-Population Inference

Research on missing mass, species discovery, population size, and the transfer of inference between sampling experiments.

The organizing question is when sampling kernels turn additive observations into recoverable population functionals and valid confidence regions. Poisson independence supplies linear cumulant functionals of the abundance measure. Positive and dominated transfers carry calibration to other experiments; signed transfers carry expectations and expose instability. The German tank problem supplies an exact finite decision example. Tag–retag and inverse-coupon inference are part of the research agenda, with their experiment-specific assumptions still to be developed here.

**Status:** research drafts and computational verification; not peer-reviewed publications. This repository does not establish publication priority, universal minimaxity, or formal proof-assistant verification.

## Manuscripts and current results

| File | Version | Content |
|---|---|---|
| [Discovery-rate inversion](manuscripts/discovery_rate_inversion.tex) | 0.2, November 12, 2026 | Sampling calculus; moment recovery; simultaneous envelopes; Hill-order anomaly; integer realizability; joint calibration; exact tank game |
| [Transfer and generation](manuscripts/transfer_generation.tex) | 0.1, November 12, 2026 | Positive/signed transfers; random-effort deficiency; dominated deletion; fixed-size prediction; batch extraction |
| [Additive cumulants and fixed-sample computation](research/master_identity_and_deletion.tex) | 0.1, November 16, 2026 | Explicit integrability conditions; polynomial deletion algorithms; threshold costs; effort-optimized full-count calibration |

The first two files preserve the most recent versions recovered from prior work. The supplied discovery attachment was version 0.1; the recovered version 0.2 already implements most of the review suggestions. See [source provenance](provenance/SOURCES.md), the [current audit](audit/AUDIT.md), and the [literature audit](audit/LITERATURE.md).

The continuation note proves a useful normalization result. When conditioning different Poisson efforts gives the same fixed-count data score, one can minimize its deterministic cumulant normalizer over effort without splitting the error level. The resulting e-value remains convex in the candidate measure. This does not justify maximizing arbitrary deletion e-values or choosing data-dependent score directions.

## Reproduce

Use Python 3.12 or later. Install the recorded numerical environment:

~~~sh
python -m pip install -r requirements.txt
make check
~~~

The three check targets run independently:

~~~sh
make check-discovery
make check-transfer
make check-new
~~~

Discovery verification distinguishes exact rational certificates, symbolic identities, and floating-point examples. Transfer verification reports 8,350 checks. The new deletion and normalization verifier reports 2,150 checks, using direct subset enumeration and finite-population expectations. Counts overlap in purpose and are not a measure of mathematical assurance.

To compile all three documents, install a TeX distribution with latexmk, Latin Modern, amsmath, microtype, cleveref, and the other packages declared in the sources:

~~~sh
make pdf
~~~

PDFs are written under build/. The baseline papers describe their original bundle paths; the corresponding code is now under verification/discovery/ and verification/transfer/.

## Code

[src/partial_population/deletion.py](src/partial_population/deletion.py) implements log-domain population and prediction deletion scores, plus effort optimization for a common fixed-count score. This is a research prototype. The population polynomial takes O(n²) arithmetic operations and O(n) storage; the two-sided prediction polynomial also tracks lost species. Log-domain arithmetic reduces numerical overflow but does not supply rigorous error bounds.

Directions, mixture weights, sampling windows, and population restrictions must satisfy the conditions in the papers. A lower abundance bound is needed for some recovery targets, especially unrestricted richness, but not for every calibration kernel. Relaxed positive measures and integer-realizable populations are different feasible sets.

## Continuing the research

The next concrete milestones are [listed here](research/ROADMAP.md): certified projected intervals, matched statistical benchmarks, and an explicit tag–retag/inverse-coupon experiment layer. No external dataset has yet been used to validate these methods.

Author: Christopher D. Long. Original code and text are released under the [MIT license](LICENSE). Cited literature retains its own rights; this repository does not redistribute third-party papers.
