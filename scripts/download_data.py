#!/usr/bin/env python3
"""
Data Download Script for scRNA-Seq Analysis

This script downloads the 3k PBMC dataset from 10x Genomics and prepares
it for analysis with the scRNA-Seq clustering workflow.

Usage:
    python scripts/download_data.py
"""

import os
import sys
import tarfile
import urllib.request
from pathlib import Path
from typing import Optional


# Configuration
DATA_DIR = Path("data/pbmc3k")
DATA_URL = "http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
ARCHIVE_NAME = "pbmc3k_filtered_gene_bc_matrices.tar.gz"
EXPECTED_PATH = DATA_DIR / "filtered_gene_bc_matrices" / "hg19"


# ANSI color codes
class Colors:
    RED = '\033[0;31m'
    GREEN = '\033[0;32m'
    YELLOW = '\033[1;33m'
    NC = '\033[0m'  # No Color


def print_info(message: str) -> None:
    """Print info message in green."""
    print(f"{Colors.GREEN}[INFO]{Colors.NC} {message}")


def print_warn(message: str) -> None:
    """Print warning message in yellow."""
    print(f"{Colors.YELLOW}[WARN]{Colors.NC} {message}")


def print_error(message: str) -> None:
    """Print error message in red."""
    print(f"{Colors.RED}[ERROR]{Colors.NC} {message}")


def download_progress(block_num: int, block_size: int, total_size: int) -> None:
    """Show download progress."""
    downloaded = block_num * block_size
    percent = min(downloaded * 100.0 / total_size, 100)
    bar_length = 50
    filled = int(bar_length * percent / 100)
    bar = '█' * filled + '-' * (bar_length - filled)

    sys.stdout.write(f'\r  Progress: |{bar}| {percent:.1f}%')
    sys.stdout.flush()

    if downloaded >= total_size:
        print()  # New line after completion


def check_existing_data() -> bool:
    """
    Check if data already exists.

    Returns:
        True if should proceed with download, False otherwise
    """
    if EXPECTED_PATH.exists():
        print_warn(f"Data directory already exists at: {EXPECTED_PATH}")
        response = input("Do you want to re-download and overwrite? (y/N): ").strip().lower()

        if response in ['y', 'yes']:
            print_info("Removing existing data...")
            import shutil
            shutil.rmtree(DATA_DIR)
            return True
        else:
            print_info("Keeping existing data. Exiting.")
            return False

    return True


def create_data_directory() -> None:
    """Create data directory if it doesn't exist."""
    print_info(f"Creating data directory: {DATA_DIR}")
    DATA_DIR.mkdir(parents=True, exist_ok=True)


def download_dataset() -> Optional[Path]:
    """
    Download the PBMC 3k dataset.

    Returns:
        Path to downloaded archive, or None if download failed
    """
    archive_path = DATA_DIR / ARCHIVE_NAME

    print_info("Downloading PBMC 3k dataset from 10x Genomics...")
    print_info(f"Source: {DATA_URL}")
    print()

    try:
        urllib.request.urlretrieve(DATA_URL, archive_path, download_progress)
        print_info("Download complete.")
        return archive_path
    except Exception as e:
        print_error(f"Download failed: {e}")
        return None


def extract_archive(archive_path: Path) -> bool:
    """
    Extract the downloaded archive.

    Args:
        archive_path: Path to the archive file

    Returns:
        True if extraction successful, False otherwise
    """
    print()
    print_info("Extracting archive...")

    try:
        with tarfile.open(archive_path, 'r:gz') as tar:
            tar.extractall(path=DATA_DIR)

        if EXPECTED_PATH.exists():
            print_info("Extraction successful.")
            return True
        else:
            print_error(f"Extraction failed. Expected path not found: {EXPECTED_PATH}")
            return False

    except Exception as e:
        print_error(f"Extraction failed: {e}")
        return False


def display_dataset_info() -> None:
    """Display information about the extracted dataset."""
    print()
    print_info("Dataset contents:")
    print(f"  - Location: {EXPECTED_PATH}")
    print("  - Files:")

    for file_path in sorted(EXPECTED_PATH.iterdir()):
        size = file_path.stat().st_size
        size_mb = size / (1024 * 1024)
        print(f"    - {file_path.name} ({size_mb:.2f} MB)")


def cleanup_archive(archive_path: Path) -> None:
    """Remove the downloaded archive file."""
    print()
    print_info("Cleaning up archive file...")
    archive_path.unlink()


def main() -> int:
    """
    Main function to orchestrate the download process.

    Returns:
        Exit code (0 for success, 1 for failure)
    """
    print_info("Starting data download for scRNA-Seq analysis...")
    print()

    # Check if data already exists
    if not check_existing_data():
        return 0

    # Create data directory
    create_data_directory()

    # Download dataset
    archive_path = download_dataset()
    if archive_path is None or not archive_path.exists():
        print_error("Download failed. Archive not found.")
        return 1

    # Extract archive
    if not extract_archive(archive_path):
        return 1

    # Display information
    display_dataset_info()

    # Clean up
    cleanup_archive(archive_path)

    # Success message
    print()
    print_info("✓ Data download and extraction complete!")
    print_info("You can now run the analysis notebook:")
    print("  jupyter-lab notebooks/scrna_analysis_workflow.ipynb")
    print()
    print_info("Done.")

    return 0


if __name__ == "__main__":
    sys.exit(main())
