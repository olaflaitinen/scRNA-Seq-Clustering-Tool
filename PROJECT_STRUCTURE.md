# Project Structure

This document describes the organization and contents of the scRNA-Seq Clustering Tool project.

## Directory Tree

```
scRNA-Seq-Clustering-Tool/
├── README.md                              # Project overview and documentation
├── LICENSE                                # MIT License
├── PROJECT_STRUCTURE.md                   # This file
├── .gitignore                             # Git ignore rules
├── config.yaml                            # Configuration parameters
├── requirements.txt                       # Python dependencies
│
├── notebooks/                             # Jupyter notebooks
│   └── scrna_analysis_workflow.ipynb      # Main analysis workflow
│
├── scripts/                               # Utility scripts
│   ├── download_data.sh                   # Bash script to download data
│   ├── download_data.py                   # Python script to download data
│   └── utils.py                           # Python utility functions
│
├── data/                                  # Data directory
│   ├── .gitkeep                           # Keep directory in git
│   ├── README.md                          # Data directory documentation
│   ├── pbmc3k/                            # Raw PBMC dataset (after download)
│   │   └── filtered_gene_bc_matrices/
│   │       └── hg19/
│   │           ├── matrix.mtx
│   │           ├── genes.tsv
│   │           └── barcodes.tsv
│   └── pbmc3k_processed.h5ad              # Processed results (generated)
│
└── figures/                               # Output figures
    ├── .gitkeep                           # Keep directory in git
    ├── umap_leiden.png                    # UMAP with clusters
    ├── umap_annotated.png                 # UMAP with cell types
    ├── cell_annotations.csv               # Cell type annotations
    └── marker_genes.csv                   # Marker genes per cluster
```

## Component Descriptions

### Root Files

#### README.md
Main project documentation including:
- Project overview and goals
- Analysis workflow description
- Installation instructions
- Usage guide
- Key visualizations

#### config.yaml
Centralized configuration file containing:
- Data paths
- Quality control parameters
- Preprocessing settings
- Analysis parameters (PCA, UMAP, clustering)
- Cell type annotation mappings
- Visualization settings

#### requirements.txt
Python package dependencies:
- Core: scanpy, anndata
- Data processing: numpy, pandas, scipy
- Visualization: matplotlib, seaborn
- Clustering: leiden, umap-learn
- Environment: jupyterlab, notebook

### Notebooks

#### scrna_analysis_workflow.ipynb
Complete single-cell RNA-Seq analysis pipeline:

**Section 1-2: Setup and Data Loading**
- Import libraries and configure environment
- Load 10x Genomics formatted data

**Section 3: Quality Control**
- Calculate QC metrics (genes/cell, counts, mitochondrial %)
- Visualize QC distributions
- Filter low-quality cells and genes

**Section 4-6: Preprocessing**
- Normalize and log-transform data
- Identify highly variable genes
- Scale and regress out technical variation
- Perform PCA

**Section 7-9: Clustering and Visualization**
- Compute neighborhood graph
- Generate UMAP embedding
- Perform Leiden clustering

**Section 10-11: Biological Interpretation**
- Identify marker genes per cluster
- Annotate cell types based on known markers
- Create publication-quality visualizations

**Section 12-14: Export Results**
- Save processed data
- Export annotations and marker genes
- Generate summary statistics

### Scripts

#### download_data.sh / download_data.py
Automated data download scripts that:
- Check for existing data
- Download 3k PBMC dataset from 10x Genomics
- Extract and organize files
- Verify data integrity
- Clean up temporary files

Both scripts perform the same function - use `.py` for cross-platform compatibility.

#### utils.py
Python utility module providing:

**AnalysisConfig Class**
- Load and manage configuration parameters
- Access settings via dot notation
- Get paths, QC params, marker genes

**Helper Functions**
- `setup_scanpy()`: Configure scanpy environment
- `apply_qc_filters()`: Apply quality control filters
- `preprocess_data()`: Standard preprocessing pipeline
- `run_dimensionality_reduction()`: PCA, neighbors, UMAP
- `perform_clustering()`: Leiden/Louvain clustering
- `find_marker_genes()`: Differential expression
- `annotate_cell_types()`: Map clusters to cell types
- `export_results()`: Save analysis outputs

### Data

#### data/pbmc3k/ (after download)
Raw 10x Genomics data:
- `matrix.mtx`: Sparse expression matrix (genes × cells)
- `genes.tsv`: Gene identifiers and symbols
- `barcodes.tsv`: Cell barcodes

#### data/pbmc3k_processed.h5ad (generated)
Processed AnnData object containing:
- Normalized expression data
- QC metrics
- Dimensionality reductions (PCA, UMAP)
- Cluster assignments
- Cell type annotations
- Marker genes

### Figures

Output directory for all visualizations:
- QC violin plots
- Scatter plots (genes vs counts)
- PCA variance ratio
- Highly variable genes plot
- UMAP embeddings (clusters and cell types)
- Marker gene heatmaps and dotplots
- Cell type annotation plots

## File Formats

### Input Formats

**10x Genomics Matrix Format**
- `matrix.mtx`: Market Matrix format (sparse)
- `genes.tsv`: Tab-separated (gene_id, gene_symbol)
- `barcodes.tsv`: Tab-separated (barcode)

### Output Formats

**HDF5 (.h5ad)**
- AnnData's native format
- Compressed hierarchical data
- Contains all analysis results

**CSV/TSV**
- Cell annotations
- Marker genes
- Summary statistics

**PNG/PDF/SVG**
- Publication-quality figures
- Configurable DPI and format

## Configuration System

The project uses a YAML-based configuration system (`config.yaml`) that allows customization of:

1. **Paths**: Data and output directories
2. **QC Parameters**: Filter thresholds
3. **Preprocessing**: Normalization, scaling
4. **Analysis**: PCA, UMAP, clustering parameters
5. **Annotation**: Cell type markers and mappings
6. **Visualization**: Figure settings
7. **Export**: Output formats and options

Benefits:
- No code modification needed for parameter changes
- Reproducible analysis with version-controlled configs
- Easy experimentation with different settings
- Centralized parameter management

## Workflow Summary

```
1. Setup
   └─> Install dependencies (requirements.txt)

2. Data Acquisition
   └─> Run download script (download_data.py)

3. Analysis
   └─> Execute Jupyter notebook (scrna_analysis_workflow.ipynb)
       ├─> Load data
       ├─> Quality control
       ├─> Preprocessing
       ├─> Dimensionality reduction
       ├─> Clustering
       ├─> Cell type annotation
       └─> Export results

4. Results
   └─> Processed data (pbmc3k_processed.h5ad)
   └─> Figures (figures/*.png)
   └─> Tables (figures/*.csv)
```

## Adding Custom Datasets

To analyze your own data:

1. Place data in `data/` directory
2. Update `config.yaml`:
   - Set `paths.data_dir` to your data location
   - Adjust QC parameters if needed
   - Modify cell type markers for your expected populations
3. Run the analysis notebook
4. Update cluster-to-cell-type mapping based on results

## Best Practices

### Code Organization
- Keep analysis logic in notebooks
- Reusable functions in `scripts/utils.py`
- Configuration in `config.yaml`
- Documentation in markdown files

### Data Management
- Raw data in `data/` (gitignored if large)
- Processed data saved as `.h5ad`
- Version control small metadata files

### Reproducibility
- Set random seeds in config
- Document software versions
- Save complete analysis parameters
- Archive processed results

### Version Control
- Commit code and configs
- Exclude large data files
- Include example outputs
- Document analysis decisions

## Additional Resources

### Documentation
- [Scanpy Documentation](https://scanpy.readthedocs.io/)
- [AnnData Documentation](https://anndata.readthedocs.io/)
- [10x Genomics Datasets](https://www.10xgenomics.com/resources/datasets)

### Tutorials
- [Scanpy Tutorial](https://scanpy-tutorials.readthedocs.io/)
- [Single-cell Best Practices](https://www.sc-best-practices.org/)

### Related Tools
- Seurat (R-based scRNA-Seq analysis)
- Cell Ranger (10x preprocessing)
- CellTypist (automated annotation)
