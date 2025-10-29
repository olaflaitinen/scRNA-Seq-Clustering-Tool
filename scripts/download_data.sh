#!/bin/bash

###############################################################################
# Data Download Script for scRNA-Seq Analysis
#
# This script downloads the 3k PBMC dataset from 10x Genomics and prepares
# it for analysis with the scRNA-Seq clustering workflow.
#
# Usage: bash scripts/download_data.sh
###############################################################################

set -e  # Exit on error

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Print functions
print_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Configuration
DATA_DIR="data/pbmc3k"
DATA_URL="http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
ARCHIVE_NAME="pbmc3k_filtered_gene_bc_matrices.tar.gz"
EXPECTED_PATH="${DATA_DIR}/filtered_gene_bc_matrices/hg19"

###############################################################################
# Main download logic
###############################################################################

print_info "Starting data download for scRNA-Seq analysis..."
echo

# Check if data already exists
if [ -d "$EXPECTED_PATH" ]; then
    print_warn "Data directory already exists at: $EXPECTED_PATH"
    read -p "Do you want to re-download and overwrite? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        print_info "Keeping existing data. Exiting."
        exit 0
    fi
    print_info "Removing existing data..."
    rm -rf "$DATA_DIR"
fi

# Create data directory
print_info "Creating data directory: $DATA_DIR"
mkdir -p "$DATA_DIR"

# Download data
print_info "Downloading PBMC 3k dataset from 10x Genomics..."
print_info "Source: $DATA_URL"
echo

if command -v curl &> /dev/null; then
    curl -L -o "${DATA_DIR}/${ARCHIVE_NAME}" "$DATA_URL" --progress-bar
elif command -v wget &> /dev/null; then
    wget -O "${DATA_DIR}/${ARCHIVE_NAME}" "$DATA_URL"
else
    print_error "Neither curl nor wget is available. Please install one of them."
    exit 1
fi

# Check if download was successful
if [ ! -f "${DATA_DIR}/${ARCHIVE_NAME}" ]; then
    print_error "Download failed. Archive not found."
    exit 1
fi

print_info "Download complete."
echo

# Extract archive
print_info "Extracting archive..."
tar -xzf "${DATA_DIR}/${ARCHIVE_NAME}" -C "$DATA_DIR/"

# Verify extraction
if [ -d "$EXPECTED_PATH" ]; then
    print_info "Extraction successful."

    # Display file information
    echo
    print_info "Dataset contents:"
    echo "  - Location: $EXPECTED_PATH"
    echo "  - Files:"
    ls -lh "$EXPECTED_PATH" | tail -n +2 | awk '{print "    - " $9 " (" $5 ")"}'

    # Clean up archive
    echo
    print_info "Cleaning up archive file..."
    rm "${DATA_DIR}/${ARCHIVE_NAME}"

    echo
    print_info "✓ Data download and extraction complete!"
    print_info "You can now run the analysis notebook:"
    echo "  jupyter-lab notebooks/scrna_analysis_workflow.ipynb"
else
    print_error "Extraction failed. Expected path not found: $EXPECTED_PATH"
    exit 1
fi

echo
print_info "Done."
