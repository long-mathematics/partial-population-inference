PYTHON ?= python
LATEXMK ?= latexmk
.PHONY: check check-discovery check-transfer check-new check-capture check-projection check-empirical check-pdfs pdf publish-pdfs
check: check-discovery check-transfer check-new check-capture check-projection check-empirical check-pdfs
check-discovery:
	$(PYTHON) verification/discovery/run_checks.py
check-transfer:
	$(PYTHON) verification/transfer/run_all.py
check-new:
	PYTHONPATH=src $(PYTHON) verification/new/check_deletion.py
check-capture:
	PYTHONPATH=src $(PYTHON) verification/capture/check_capture_recapture.py
check-projection:
	PYTHONPATH=src $(PYTHON) verification/projection/check_projection.py
	PYTHONPATH=src $(PYTHON) verification/projection/check_choices.py
check-empirical:
	PYTHONPATH=src $(PYTHON) verification/empirical/check_artifacts.py
	$(PYTHON) scripts/empirical_tables.py --check
check-pdfs:
	$(PYTHON) scripts/pdfs.py
pdf:
	mkdir -p build
	$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error -outdir=build manuscripts/discovery_rate_inversion.tex
	$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error -outdir=build manuscripts/transfer_generation.tex
	$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error -outdir=build research/master_identity_and_deletion.tex
	$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error -outdir=build manuscripts/capture_recapture.tex
	$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error -outdir=build manuscripts/classical_data.tex
	$(PYTHON) scripts/check_latex_logs.py
publish-pdfs: pdf
	$(PYTHON) scripts/pdfs.py --publish
