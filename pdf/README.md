# Download the papers

These PDFs are compiled from the LaTeX sources in this repository.

- [Discovery-rate inversion](https://github.com/long-mathematics/partial-population-inference/raw/refs/heads/main/pdf/discovery_rate_inversion.pdf)
- [Transfer and generation](https://github.com/long-mathematics/partial-population-inference/raw/refs/heads/main/pdf/transfer_generation.pdf)
- [Additive cumulants and fixed-sample computation](https://github.com/long-mathematics/partial-population-inference/raw/refs/heads/main/pdf/master_identity_and_deletion.pdf)
- [Capture–recapture](https://github.com/long-mathematics/partial-population-inference/raw/refs/heads/main/pdf/capture_recapture.pdf)
- [Classical data: discovery, recapture, and serial numbers](https://github.com/long-mathematics/partial-population-inference/raw/refs/heads/main/pdf/classical_data.pdf)

Run `make publish-pdfs` from the repository root to rebuild and refresh the downloads. The source and PDF hashes in [source_manifest.json](source_manifest.json) are checked by `make check-pdfs` and CI. The manifest verifies which source files accompany the committed PDFs; it is not a claim of byte-identical output across TeX installations.
