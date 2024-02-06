# Makefile for setting up LongQC and its dependencies


# Define the LongQC repository URL
LONGQC_REPO := https://github.com/yfukasawa/LongQC.git

# Define the LongQC minimap2-coverage directory
MINIMAP2_COVERAGE_DIR := LongQC/minimap2-coverage

# Define the LongQC version
LONGQC_VERSION := 1.2.1

# Define the LongQC flags for Apple Silicon (M1/M2)
LONGQC_FLAGS := arm_neon=1 aarch64=1

# Define targets
.PHONY: all install_dependencies install_longqc_mac install_longqc

all: install_dependencies install_longqc_mac install_longqc

install_dependencies:
	# Install LongQC dependencies using conda
	conda install numpy scipy matplotlib scikit-learn pandas jinja2 h5py
	conda install -c bioconda pysam edlib python-edlib

install_longqc_mac:
	# For Mac users: Install argp-standalone using homebrew
	brew install argp-standalone

install_longqc:
	# Clone LongQC repository and build minimap2-coverage
	if [ -d "LongQC" ]; then rm -rf "LongQC"; fi
	git clone $(LONGQC_REPO) && \
	cd "$(MINIMAP2_COVERAGE_DIR)" && make
	cd ..

.PHONY: clean

clean:
	# Remove cloned LongQC repository and built minimap2-coverage
	rm -rf $(INSTALL_PATH)/LongQC

