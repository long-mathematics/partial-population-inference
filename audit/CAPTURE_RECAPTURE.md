# Capture–recapture audit — November 17, 2026

This is an internal proof and computational audit of the self-contained capture–recapture module. It does not claim independent referee review, formal proof-assistant verification, or recovery of the older private tag–retag manuscript.

## Model and theorem audit

| Component | Checked argument and scope |
|---|---|
| Unknown-size overlap kernel | Uniform recapture subset conditional on a marked set. Numeric label magnitudes are excluded. The unknown population-size parameterization differs from the known-population/unknown-marked-count theorem in the transfer paper. |
| Reciprocal-factorial recovery | Ordered-tuple inclusion probabilities give the moments. A triangular statistic basis and polynomial independence after clearing denominators give the full range and uniqueness. No finite-valued unbiased estimator of unbounded N exists in this finite overlap experiment. |
| Lincoln–Petersen likelihood | Consecutive likelihood ratios give the integer maximizer, including ties and support truncation. Zero overlap has no finite maximizer. Infinite-value conventions retain their consequences for risk. |
| Chapman | Vandermonde summation gives the exact residual bias −(N−a−b)P(J=0). The generalized shifted-reciprocal identity retains every omitted boundary term. The estimate is at least the observed union. Exact variance is evaluated from the finite law. |
| Bailey | Binomial second-stage sampling gives a different correction and residual bias. Its observation range allows repeated marked sightings and can exceed the number marked. |
| Signed transfer | Infinite-parameter range inclusion is necessary and sufficient. Coefficient extraction gives a unique normalized matrix. Positive deletion is checked separately; the (2,2)→(1,3) matrix has a negative entry and variation norm 2. No finite-grid classification is inferred. |
| Fixed-size limits | Coupling by adding unmarked units proves stochastic order. Discrete tail inequalities give coverage; explicit tail events witness pointwise one-sided optimality in the nonincreasing class. Strict acceptance inequalities are retained. Zero overlap requires an infinite upper endpoint. |
| Repeated captures | Conditional hypergeometric means give Schnabel and Schumacher–Eschmeyer estimating functions. A zero-mean estimating function is not asserted to give an unbiased ratio. Dynamic programming retains zero denominators. |
| Sequential calibration | Predictable directions and exact conditional normalizers give test martingales. Bounded stopping and a crossing argument establish simultaneous coverage. Impossible candidate prefixes may be rejected. Retrospective direction choice is excluded. |
| Inverse sampling | Uniform marked positions give the negative-hypergeometric law. Uniform gap compositions give exact mean and variance. Triangular support proves completeness; the conditional binary history proves sufficiency. At the last marked recapture the statistic is exactly the tank maximum. |
| Replacement inverse sampling | Negative-binomial generating functions give unbiasedness and variance. No completeness claim is made for the discrete parameter subfamily. Exact tails are used without truncating the infinite support. |
| Independent histories | The log transform is linear in a fixed-population counting measure. The normalization mass is unknown N. The proof distinguishes conditional fixed catchabilities from marginal random effects, whose transform is different. |
| Horvitz–Thompson | Inclusion weighting is unbiased with known positive probabilities. Pair inclusion probabilities enter the general variance; independence is not imposed for the mean identity. |
| Chao-type recovery | Cauchy–Schwarz applied to the exact binomial kernel retains (K−1)/K; the Poisson bound has a different coefficient. Observed frequency ratios are explicitly point estimates, not confidence bounds. |
| Heterogeneity obstruction | Zero catchability is exactly invisible; arbitrarily small positive catchability yields a coupling proof against uniformly honest, always-finite upper population bounds. |

## Computational verification

[The verifier](../verification/capture/check_capture_recapture.py) passes 54,416 exact or symbolic checks. [The machine-readable report](../verification/capture/results.json) records the families. Counts overlap in mathematical purpose and are not a measure of proof assurance.

Independent constructions include labelled subsets for the overlap law, uniform marked positions for inverse sampling, labelled repeated-capture samples for the ratio distributions, and full Bernoulli histories for the cumulant transform and conditioning identity. Coverage checks use the computed unbounded-search intervals, including strict-threshold equality examples. The negative-binomial coverage check closes the entire remaining tail analytically; it does not discard a truncation remainder. Predictable adaptive betting is tested by exhaustive finite-history expectation and crossing calculations.

The original discovery, transfer, and cumulant verification suites are retained. The integration preserves the established theorem numbering: new cross-references and displays in the existing development are unnumbered, and the new prior-work sections follow the existing mathematical sections.

## Literature audit and positioning

- **Chapman (1951):** monograph metadata verified through the digitized catalog record; the full original monograph was not available for page-level inspection. The correction is also documented in Dettloff's accessible author manuscript. Every identity used here is proved directly.
- **Bailey (1951, 1952):** publisher/journal records identify the two recapture papers. The 1951 publisher title is *On estimating the size of mobile populations from recapture data*. The binomial correction and replacement assumption are also documented by the original `recapr` implementation documentation. The exact expectation is independently derived.
- **Schnabel (1938):** publisher metadata verified, DOI 10.1080/00029890.1938.11990818. **Schumacher–Eschmeyer (1943):** attribution and formula checked against the accessible Dettloff manuscript and fisheries documentation; the original 1943 article was not inspected in full.
- **Horvitz–Thompson (1952):** original article metadata and public scan located, DOI 10.1080/01621459.1952.10483446. The inclusion and pair-inclusion variance arguments are supplied directly.
- **Chao (1987):** full primary paper inspected, especially its model, Poisson approximation, moment inequalities, and explicit Harris comparison. Its approximate Poisson argument is distinguished from the exact finite-occasion binomial calculation here.
- **Lloyd (1994):** publisher abstract and metadata verify the martingale-estimation precedent, DOI 10.1093/biomet/81.2.305.
- **Chao–Yip–Lee–Chu (2001):** publisher and author-institution records explicitly describe a unified estimating-function/sample-coverage approach. Its full text was not accessible in this pass; we do not claim a theorem-by-theorem exclusion of overlap with our algebraic results.
- **Pérez-Salvador (2026):** publisher text verifies the inverse-hypergeometric estimator and variance, approximate intervals, and the separate stability stopping proposal; *Sequential Analysis* 45(3), 417–437, DOI 10.1080/07474946.2026.2620792. We attribute the estimator and do not infer optional-stopping validity from its fixed-recapture-count calculation.
- **Dettloff (2023):** publisher and NOAA records give *Fisheries Research* 265, **106756**, DOI 10.1016/j.fishres.2023.106756. The accessible author manuscript gives the classical formulas and finite-design comparisons. Its simulation conclusions are not promoted to universal optimality theorems.

The detailed mathematical comparisons appear in the new paper's “Relation to prior work” section. The other three papers now use the same placement, with necessary citations retained at use. No first-unification or first-appearance claim is made for the classical formulas, moment arguments, or explicit algebraic consequences.

## Build and distribution

All four LaTeX papers are compiled with halt-on-error. The log checker rejects unresolved citations/references, warnings, and overfull or underfull boxes. Downloadable PDFs are committed under [pdf/](../pdf/README.md), with source/PDF hashes in a manifest checked by CI. CI also rebuilds the PDFs independently.

## Remaining scope

The proofs cover the declared closed-population models and the fixed-catchability history model. They do not solve open-population estimation, tag loss or match error, arbitrary behavioral dependence, certified continuum projection, or the bounded-heterogeneity inverse-coupon research target. Exact arithmetic is intended for moderate designs and does not make the current implementation a large-scale statistical package.
