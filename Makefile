PYTHON ?= python
LATEXMK ?= latexmk
.PHONY: check check-discovery check-transfer check-new pdf
check: check-discovery check-transfer check-new
check-discovery:
	$(PYTHON) verification/discovery/run_checks.py
check-transfer:
	$(PYTHON) verification/transfer/run_all.py
check-new:
	PYTHONPATH=src $(PYTHON) verification/new/check_deletion.py
pdf:
	mkdir -p build
	$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error -outdir=build manuscripts/discovery_rate_inversion.tex
	$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error -outdir=build manuscripts/transfer_generation.tex
	$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error -outdir=build research/master_identity_and_deletion.tex
