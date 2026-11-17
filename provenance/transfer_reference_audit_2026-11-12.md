# Reference and source-access audit

Checked November 12, 2026. References are primary author, publisher, or journal sources. URLs below are supplied as a research bibliography, not as evidence that every historical full text was available. No external full-text copyrighted papers are redistributed in the bundle.

## Directly relevant prior work

- **Blackwell (1953), Equivalent Comparisons of Experiments**, Annals of Mathematical Statistics 24(2), 265–272, DOI 10.1214/aoms/1177729032. Publisher landing page and author bibliographic record checked. The publisher page did not expose the historical full text through this retrieval. Cited for classical experiment comparison; the finite theorem here is proved independently.
  https://doi.org/10.1214/aoms/1177729032
- **Le Cam (1964), Sufficiency and Approximate Sufficiency**, Annals of Mathematical Statistics 35(4), 1419–1455, DOI 10.1214/aoms/1177700372. Publisher metadata and author's Berkeley bibliography checked. Historical full-text extraction was unavailable. The transfer inequalities used in the new manuscript are proved there.
  https://doi.org/10.1214/aoms/1177700372
  https://www.stat.berkeley.edu/~rice/LeCam/biblio.html
- **Torgersen (1991), Comparison of Statistical Experiments**, Cambridge University Press, DOI 10.1017/CBO9780511666353. Publisher book record checked. Used as foundational comparison-theory context, not as an asserted full-book line-by-line audit.
  https://doi.org/10.1017/CBO9780511666353
- **Kühn and Kunze (2022), Feller generators with measurable lower order terms**, Positivity 26, article 85, DOI 10.1007/s11117-022-00948-4. Publisher full-text record and author arXiv PDF checked. Section 2 explicitly discusses signed kernel operators and bounded-pointwise/weak continuity. The new manuscript proves its needed representation directly.
  https://link.springer.com/article/10.1007/s11117-022-00948-4
  https://arxiv.org/abs/2108.02542
- **Liu and Pego (2016), On generating functions of Hausdorff moment sequences**, Transactions AMS 368, 8499–8518, DOI 10.1090/tran/6618. Section 1 of the author version explicitly states the Hausdorff complete-monotonicity criterion. Version 4, November 12, 2025, includes corrections to other results; the manuscript uses only the classical characterization in Section 1.
  https://arxiv.org/html/1401.8052v4
  https://doi.org/10.1090/tran/6618
- **Last and Penrose, Lectures on the Poisson Process**, Cambridge University Press, DOI 10.1017/9781316104477. Author PDF and publisher record checked; marking, independence, and thinning results in Chapter 5 checked. The bibliography uses the publisher's 2017 publication date (the author PDF also carries 2018 copyright information).
  https://stoch.math.kit.edu/img/Last/lastpenrose2017.pdf
  https://doi.org/10.1017/9781316104477
- **Strassen (1965), The existence of probability measures with given marginals**, Annals of Mathematical Statistics 36(2), 423–439, DOI 10.1214/aoms/1177700153. Historical publisher record checked; full-text extraction was not available. The special one-sided count-coupling formula in the new manuscript has its own complete quantile proof, and is cited only as part of classical coupling theory.
  https://doi.org/10.1214/aoms/1177700153
- **Sion (1958), On general minimax theorems**, Pacific Journal of Mathematics 8(1), 171–176, DOI 10.2140/pjm.1958.8.171. Full journal PDF and the one-compact-side corollary on printed page 174 inspected. Applied only with a compact convex primal moment body and continuous affine payoffs.
  https://msp.org/pjm/1958/8-1/pjm-v8-n1-p14-p.pdf

## Poissonization, amplification, and calibration

- **Low and Zhou (2007), A complement to Le Cam's theorem**, Annals of Statistics 35(3), 1146–1165, DOI 10.1214/009053607000000091. Author arXiv PDF checked. Important precedent for the dependence of Poissonization equivalence on restrictions of the population class.
  https://arxiv.org/abs/0708.2233
- **Axelrod, Garg, Sharan, and Valiant (2020), Sample Amplification: Increasing Dataset Size even when Learning is Impossible**, ICML/PMLR 119, 442–451. Publisher metadata, abstract, and model scope checked. Explicitly concerns amplification on restricted distribution classes, including bounded support; cited to avoid presenting experiment-based sample generation as a new general framework.
  https://proceedings.mlr.press/v119/axelrod20a.html
- **Axelrod, Garg, Han, Sharan, and Valiant (2024), On the Statistical Complexity of Sample Amplification**, Annals of Statistics 52(6), 2767–2790, DOI 10.1214/24-AOS2444. Journal metadata verified against IMS and author records. Author version 2201.04315v2 inspected, including the explicit Le Cam comparison formulation and Example A.16/Lemma A.17 on Poissonized discrete sampling. These are substantial nearby results; the new manuscript does not claim to supersede their rates.
  https://arxiv.org/html/2201.04315v2
  https://doi.org/10.1214/24-AOS2444
- **Polyanskiy and Wu, Dualizing Le Cam's method for functional estimation, with applications to estimating the unseens**, arXiv:1902.05616v2, January 4, 2021 (first version 2019). Author record and text checked. Cited as a preprint in the consulted version, not assigned an unverified journal publication.
  https://arxiv.org/abs/1902.05616v2
- **Wasserman, Ramdas, and Balakrishnan (2020), Universal inference**, PNAS 117(29), 16880–16890, DOI 10.1073/pnas.1922664117. Full public article/metadata and author preprint record checked. Used for nonnegative mean-one calibration and confidence inversion.
  https://pmc.ncbi.nlm.nih.gov/articles/PMC7382245/
  https://arxiv.org/abs/1912.11436
- **Howard, Ramdas, McAuliffe, and Sekhon (2021), Time-uniform, nonparametric, nonasymptotic confidence sequences**, Annals of Statistics 49(2), 1055–1080, DOI 10.1214/20-AOS1991. Author arXiv record and companion's verified publication details checked. Used for the distinction between fixed-time e-values and sequentially valid procedures.
  https://arxiv.org/abs/1810.08240
- **Eriksson (2026), The Unseen Species Problem Revisited**, arXiv:2602.08769v4, September 9, 2026. Version-specific text checked, including bounded-arity set observations and prediction scope. No singleton-independence formula is imported into the batch model.
  https://arxiv.org/html/2602.08769v4

## Internal companion

- **Long, Shape-Constrained Discovery-Rate Inversion**, version 0.2, November 12, 2026. Source inspected and compiled to recover exact theorem numbers and page locations. The unchanged source/PDF are included in `companion/`, with hashes in `SHA256SUMS.txt`. It is an internal research draft, not cited as published prior literature.

## Novelty and completeness of the search

The searches covered experiment comparison, signed and normal kernel operators, random sample-size comparison, Poissonization, sample amplification, and bounded-batch discovery. No priority claim is made for the exact unrestricted random-effort formula or the specialized fixed-sample joint calibration. In particular, an exhaustive comparison with every result in the historical comparison-of-experiments monographs was not possible from available full texts. The primary mathematical results in the new manuscript are therefore proved self-contained, with their precise scope, rather than advertised as historically first.
