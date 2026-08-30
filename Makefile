.PHONY: setup pipeline dashboard

VENV := .venv
PYTHON := $(VENV)/bin/python
PIP := $(VENV)/bin/pip
STREAMLIT := $(VENV)/bin/streamlit

# Creates an isolated virtualenv and installs all dependencies into it.
# Using a venv (rather than installing into the system/base Python) avoids
# PEP 668 "externally managed environment" failures on some Codespaces base
# images, and keeps this project's dependency versions isolated.
setup:
	python3 -m venv $(VENV)
	$(PIP) install --upgrade pip
	$(PIP) install -r requirements.txt

# Runs the full pipeline end-to-end with no manual steps:
# initializes + loads the database (Part 1), then generates the Part 2-4
# tables and the Part 3 boxplot into output/.
pipeline:
	$(PYTHON) load_data.py
	$(PYTHON) run_analysis.py

# Starts the interactive dashboard. Binds to 0.0.0.0 so the port forwards
# correctly under GitHub Codespaces; Codespaces will prompt to open the
# forwarded port 8501 in the browser.
dashboard:
	$(STREAMLIT) run dashboard.py --server.port 8501 --server.address 0.0.0.0
