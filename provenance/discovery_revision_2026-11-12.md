# Revision 0.2 — changes and audit disposition

**Manuscript:** Shape-Constrained Discovery-Rate Inversion  
**Author:** Christopher D. Long  
**Revision date:** November 12, 2026
**Compiled manuscript:** 35 pages; 25 bibliography entries; no unresolved references or LaTeX warnings.

This is a substantive revision of v0.1, not a claim of independent peer review, formal verification, or universal statistical optimality. The original draft was preserved separately. Classical sources and the specialized results proved in this manuscript are distinguished throughout.

## Main restructuring

Section 7 now centers on a joint population–target calibration theorem. A fixed finite mixture of exact exponential identities yields a convex region for the unknown abundance measure and an additive realized or future target at a single error level. The candidate population class is explicitly fixed, nonrandom, and true-containing. A previously estimated confidence region cannot be intersected at no coverage cost.

Section 8 gives the three separate kernels: fingerprint with realized missing mass; fingerprint with actual future discoveries; and fingerprint with finite-coefficient prediction error. Countable populations and the limits at probability zero are justified explicitly. Direct future-count prediction remains exponentially calibratable even when the full long-horizon Good–Toulmin error is not. These are unconditional guarantees for predeclared targets, not conditional or automatic horizon-uniform prediction guarantees.

## Mathematical improvements

- The Hill-envelope range is extended to all finite real orders in `(-infinity,2] union [3,infinity)` on a positive compact rate interval. The sufficient condition `t*a >= 1/4` gives the common ordering for every finite real order.
- For the exceptional crossing case `2<q<3`, a contact-count theorem proves that two one-parameter families exhaust the relaxed extrema: two-node measures and endpoint–interior–endpoint measures. The proof distinguishes exhaustive analytical reduction from numerical global certification and from integer-realizable population optimization.
- The multi-time even-block envelope theorem now has an explicit real-rooted Hermite-remainder sign lemma. Karlin–Studden and the exact parts of Pinelis Proposition 2 are cited for principal representations; the endpoint orientation is proved rather than inferred from naming conventions.
- Projected Chernoff optimization now has an equality theorem under compact convex assumptions. The change `v=1/lambda` makes the objective convex in `v`, so Sion's theorem applies. Actual-population nonconvexity and infinite exponential transforms are expressly excluded.
- Long-horizon alternatives now include fourth moments and an exact truncation-failure transform. The manuscript gives a counterexample to uniform fourth-moment improvement, explains the nonlinear fourth-moment projection, and includes the error budget needed for a data-dependent truncation gate.

## Corrections and clarifications

- The tank upper limit is explicitly shown to be nondecreasing using the extended distribution function, not incorrectly claiming that the sets restricted by `N>=m` are nested.
- The finite tank rule is clipped to the observed feasible range `[M,12]`. Coverage remains exactly `19/20`; the four final risks are `1369/336`, `1657/480`, `167/60`, and `167/80`. The minimax value remains `697/140`, attained at `N=5,7`; the same analytic dual certifies it. The revision does not claim admissibility merely from clipping.
- Missing-mass Bernstein intervals use `1/t` for the positive error tail and `b` for the negative error tail. The proof uses a valid one-sided exponential inequality. Taking `b=1` needs no added assumption.
- The two-moment confidence region is explicitly the `d=2` deviance region with threshold `log(4/alpha)`, intersected with deterministic feasibility.
- The introduction states target-specific assumptions. A known positive lower probability is used for compact envelope and realizability results, not imposed on all calibration results.
- Esty's fixed-sample studentization includes the `-Phi_1^2/n^3` correction. It is not identified with the exact Poissonized variance.

## Reference changes

The bibliography has 25 entries (seven additions). Added Karlin–Studden, Esty, Zhang–Zhang, Wu–Yang, Sion, Boucheron–Lugosi–Massart, and the published 2026 general-theory article by Polyanskiy–Wu. The combined Polyanskiy–Wu preprint remains separately cited for applications. Eriksson is pinned to v4 of September 9, 2026 and its uniqueness claim is stated only within its synchronous-expeditions/atemporality class. O'Neill's issue and pagination are included; Daley–Smith is pinned to v3. See `reference_audit.md` for what was checked and the access limitations.

## Verification

All six programs in `run_checks.py` passed. Legacy suites remain separately identifiable. The current checker writes `checks/revision_v02.json` and the exact `checks/finite_game_certificate_clipped.json`; it does not silently overwrite the old untrimmed certificate.

The current checker verifies 603 exact rational constraints for the displayed tank procedure; 360 enumerated joint-transform comparisons; three symbolic zero-rate limits; fourth-moment identities and the counterexample; 600 Hill-order comparisons including negative orders; the intermediate-order families against finite-grid linear programs; 36 truncation calculations; the perspective identity; and 909 numerical instances of the one-sided exponential inequality. Numerical checks are corroboration, not proof, and the family optimizer is not advertised as interval-certified.

## Scope retained

No new minimax-rate, first-unification, conditional-coverage, or universal-estimator-generator claim is made. A certified transfer-and-generation theorem remains the next central target; the discussion distinguishes that program from results proved in this revision.
