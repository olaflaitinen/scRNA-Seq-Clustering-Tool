# Data Directory

This directory is used to store the scRNA-Seq datasets and processed results.

## Structure

```
data/
├── pbmc3k/                           # Raw data from 10x Genomics
│   └── filtered_gene_bc_matrices/
│       └── hg19/
│           ├── matrix.mtx            # Expression matrix
│           ├── genes.tsv             # Gene annotations
│           └── barcodes.tsv          # Cell barcodes
└── pbmc3k_processed.h5ad             # Processed AnnData object (output)
```

## Downloading Data

The 3k PBMC dataset from 10x Genomics can be downloaded using the provided helper scripts:

### Option 1: Using Python script (recommended)
```bash
python scripts/download_data.py
```

### Option 2: Using Bash script
```bash
bash scripts/download_data.sh
```

### Option 3: Manual download
```bash
# Create directory
mkdir -p data/pbmc3k/

# Download the dataset
curl -o data/pbmc3k_filtered_gene_bc_matrices.tar.gz \
  http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz

# Extract
tar -xzf data/pbmc3k_filtered_gene_bc_matrices.tar.gz -C data/pbmc3k/

# Clean up
rm data/pbmc3k_filtered_gene_bc_matrices.tar.gz
```

## Dataset Information

**Dataset:** 3k PBMCs from a Healthy Donor
**Source:** 10x Genomics
**Technology:** Chromium Single Cell 3' v1
**Cells:** ~3,000 peripheral blood mononuclear cells
**Genes:** ~32,000 genes (human genome hg19)

## Expected Cell Types

This PBMC dataset typically contains the following cell types:
- CD4+ T cells
- CD8+ T cells
- B cells
- NK cells
- Monocytes (CD14+ and FCGR3A+)
- Dendritic cells
- Megakaryocytes

## File Formats

- **matrix.mtx**: Market Matrix format containing the gene expression counts
- **genes.tsv**: Tab-separated file with gene identifiers and symbols
- **barcodes.tsv**: Tab-separated file with cell barcodes
- **pbmc3k_processed.h5ad**: HDF5-based AnnData format containing the complete processed dataset

## Notes

- Large data files (`.h5ad`, `.mtx`, `.tar.gz`) are excluded from version control via `.gitignore`
- The processed `.h5ad` file will be generated after running the analysis notebook
- For different datasets, adjust the data path in `config.yaml`
