# Audit and research continuation — November 16, 2026

## Version finding

The attached discovery draft was version 0.1. The prior-context revision and its complete verification bundle contain version 0.2, dated November 12, 2026. The transfer bundle also contains that same discovery version. Both recovered copies agree byte-for-byte.

The attached transfer draft agrees byte-for-byte with the recovered version 0.1. Searches of the available prior-file records found no later matching transfer version. This identifies the latest recovered sources, not a guarantee that no unindexed private copy exists.

The repository preserves the supplied baseline in its first commit and the recovered versions in the next commit. Current manuscript sources are unmodified recovered baselines. New results appear in a separate continuation note.

## Disposition of the supplied review

| Item | Latest-source status / audit conclusion |
|---|---|
| Tank upper endpoint monotonicity | Already proved in discovery v0.2 by nested acceptance sets. |
| Even-block extremizer orientation | v0.2 supplies a direct Hermite-sign proof and cites Karlin–Studden and Pinelis. Pinelis Proposition 2 does contain the relevant extremal-measure statement. |
| Negative Hill orders | Already included for a positive lower support endpoint; the statement covers q ≤ 2 and q ≥ 3. |
| Intermediate Hill orders | v0.2 already proves the two exhaustive one-parameter contact families in the single-sign-change case. Numerical grid agreement is corroboration, not the proof. |
| Two-moment deviance region | Clarified in v0.2; reduction is restricted to constraints depending only on those moment coordinates. |
| Clipped tank intervals | v0.2 uses the clipped rule. All coverages are 19/20, the four improved risks match the review, and the exact optimum remains 697/140. |
| Joint cumulant calibration | Already central to v0.2; a fixed true-containing candidate class and integrable kernels are required. It gives joint, not conditional, target coverage. |
| Missing-mass variance citation | v0.2 correctly distinguishes Poisson variance from Esty's fixed-sample correction. |
| Long-horizon Good–Toulmin | v0.2 includes both fourth moments and truncation. The fourth-moment radius can be much worse at a fixed error level; the explicit counterexample reruns. |
| Chernoff order exchange | v0.2 proves equality over the compact convex relaxation after the perspective substitution. The nonconvex integer class gets only the conservative inequality. |
| Lower abundance assumption | Discussed up front for relevant recovery tasks. Continuous calibration kernels can extend to zero without identifying unrestricted richness. |
| Transfer computational cost | The continuation note supplies the O(n²) population algorithm and a lost-species bivariate extension for prediction. |
| Statistical cost of deletion | The continuation note makes the full-count log-threshold bound explicit. This is a deterministic region comparison, not a universal power theorem. |

No fatal contradiction was found in the proof and scope checks performed here. This is an internal audit, not an independent referee report or a formal verification of every theorem.

## New proof-level work

The continuation note contains:

1. A countable-species additive cumulant theorem with explicit sufficient exponential-integrability conditions and justified differentiation.
2. Polynomial coefficient algorithms for population calibration and the ordinary-draw two-sided novelty correction.
3. An explicit threshold cost of minus log Poisson(n; t), asymptotic to one half log(2πn) at t=n; two-window costs add.
4. A common-score normalization theorem: optimize deterministic upper bounds over effort before normalizing the shared fixed-count statistic. Concavity of the optimized cumulant normalizer preserves convex joint confidence regions.

The last result permits a maximum only within the special family of deterministic multiples of the same random score. It does not validate generic e-value maximization or a data-selected direction. No literature-priority claim is made.

## Verification performed

All six discovery programs passed: the retained sampling, recovery, and calibration programs, the independent finite-game audit, prediction variance checks, and the v0.2 revision checker.

The clipped finite game passed 603 exact rational checks. The revision checker reproduced 360 joint-transform comparisons, three symbolic zero-rate limits, 600 Hill comparisons, and 36 truncation comparisons. The intermediate-Hill grid discrepancy is about 2.74×10⁻⁹; this is a discretization comparison, not a continuum optimization error bound.

The transfer suite passed 8,350 checks: 2,182 retained exact checks and 6,168 extension checks. Their arithmetic types are recorded by the scripts.

The new verifier passed 2,150 checks. It compares polynomial coefficients against rational indexed-subset enumeration, and evaluates the normalized statistics over finite fixed-size sample laws. Maximum relative discrepancy in those coefficient comparisons was below 10⁻¹⁵. This does not certify roundoff for larger inputs or prove coverage outside the tested cases; the analytic proofs supply the general statements.

Environment: Python 3.12.14; NumPy 2.3.5; SciPy 1.17.0; SymPy 1.14.0. JSON reports are retained alongside their scripts. Historical and clipped tank certificates have distinct filenames and should not be substituted for one another.

## Document build

All three LaTeX sources compiled successfully: discovery 35 pages, transfer 23 pages, and continuation 7 pages. The final logs contain no unresolved citations, references, or overfull boxes. The continuation pages and representative baseline pages were visually inspected. GitHub Actions is configured, but no remote CI run is claimed before the repository is created.

## Remaining limits

The work does not yet supply a certified production continuum optimizer, a real-data width study, or the recovered tag–retag manuscript. The repository's name describes the program's scope, not completed sharp solutions for every population problem.

The earlier drafts remain separately versioned. Full consolidation, theorem-level novelty review, and target-matched statistical comparisons are the next research steps.
