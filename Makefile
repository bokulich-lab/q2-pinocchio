# Makefile for setting up LongQC, its dependencies, and additional Python development tasks

# Define the LongQC repository URL
LONGQC_REPO := https://github.com/yfukasawa/LongQC.git

# Define the LongQC minimap2-coverage directory
MINIMAP2_COVERAGE_DIR := LongQC/minimap2-coverage

# Define the LongQC version
LONGQC_VERSION := 1.2.1

# Define the LongQC flags for Apple Silicon (M1/M2)
LONGQC_FLAGS := arm_neon=1 aarch64=1

# Python executable
PYTHON ?= python

# Define targets
.PHONY: install_dependencies install_longqc_mac install_longqc lint test test-cov install dev clean distclean

# New 'all' target that does not include the installation of dependencies
all: install_longqc_mac install_longqc

install_dependencies:
	# Install LongQC dependencies using conda
	conda install numpy scipy matplotlib scikit-learn pandas jinja2 h5py
	conda install -c bioconda pysam
	conda install -c bioconda edlib
	conda install -c bioconda python-edlib

install_longqc_mac:
	# For Mac users: Install argp-standalone using homebrew
	brew install argp-standalone

install_longqc:
	# Clone LongQC repository and build minimap2-coverage
	if [ -d "LongQC" ]; then rm -rf "LongQC"; fi
	git clone $(LONGQC_REPO) && \
	cd "$(MINIMAP2_COVERAGE_DIR)" && make
	cd ..

lint:
	q2lint
	flake8

test:
	py.test

test-cov:
	coverage run -m pytest
	coverage xml

install:
	$(PYTHON) setup.py install

dev:
	# Only install the necessary tools for development without reinstalling dependencies
	pip install pre-commit
	pip install -e .
	pre-commit install

clean: distclean
	# Remove cloned LongQC repository and built minimap2-coverage
	rm -rf $(INSTALL_PATH)/LongQC

distclean:
	# Additional clean commands for development environment
