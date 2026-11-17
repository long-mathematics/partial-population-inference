# Recovered source provenance

The November 16 audit compared the attached files with the prior-context sources and full verification bundles. The most recent matching discovery source found was version 0.2, dated November 12. The transfer source remains version 0.1 of the same date.

The discovery v0.2 bundle and the companion copied into the transfer bundle agree byte-for-byte. The attached discovery source is the earlier version 0.1. The attached and recovered transfer sources agree byte-for-byte.

## Source digests

| Source | Bytes | SHA-256 |
|---|---:|---|
| supplied discovery v0.1 | 79075 | 7805a13ad318d14fbb91c796db1c92e202fd90c7d53f8b0904e2b4963870137a |
| recovered discovery v0.2 | 108935 | 542c9fbdf8fcf199c5b2c59dbe1dfa7f1db11a5df63ea5219e8b4a4813ca558c |
| supplied transfer v0.1 | 71871 | 78678dfef9d2ac02293c631c4f69cb6cf49c40cb90e8c6074191c2a3b7a213d2 |
| recovered transfer v0.1 | 71871 | 78678dfef9d2ac02293c631c4f69cb6cf49c40cb90e8c6074191c2a3b7a213d2 |

## Retained material

- The discovery bundle supplies the three original verification programs, the independent audit, the Good–Toulmin variance checker, the revision checker, and the historical and clipped game certificates.
- The transfer bundle supplies its retained finite checks, extension checks, and runner.
- Historical reference audits and the discovery revision memo are retained here under filenames dated November 12.
- The recovered versions are retained in git history. The November 17 integration updates cross-references and the placement of prior-work discussions while preserving established theorem numbering. Original verification-bundle paths correspond to `verification/discovery/` and `verification/transfer/`.
- Original standalone PDFs and duplicate source copies were omitted from the source tree; the PDFs can be rebuilt.

The first local commit preserves the supplied sources. The second records recovery of the newer discovery draft. Subsequent additions are the present audit, repository tooling, and a separate continuation note.

Library identifiers are omitted from the research repository because they are account-specific retrieval metadata. Source content hashes provide a portable identity check.
