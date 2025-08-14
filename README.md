# scRNA-Seq Clustering and Annotation Tool

[![License: MIT](https://img.shields.io/badge/License-MIT-purple.svg)](https://opensource.org/licenses/MIT)
![Language](https://img.shields.io/badge/language-Python-blue.svg)
[![Built with Scanpy](https://img.shields.io/badge/built%20with-Scanpy-brightgreen.svg)](https://scanpy.readthedocs.io/en/stable/)

An implementation of a standard single-cell RNA-Seq (scRNA-Seq) analysis workflow using [Scanpy](https://scanpy.readthedocs.io/en/stable/). This project demonstrates the key steps from raw count data to annotated cell type clusters.

![Example UMAP Plot](assets/umap_cell_types_example.png)

### Workflow Overview

The pipeline follows the best practices for scRNA-Seq analysis as established by the community.

```
┌────────────────────────┐
│     Raw Count Matrix   │
└──────────┬─────────────┘
           │
           ▼
┌────────────────────────┐
│   QC & Filtering       │ (Filter low-quality cells/genes)
└──────────┬─────────────┘
           │
           ▼
┌────────────────────────┐
│ Normalization & Scaling│
└──────────┬─────────────┘
           │
           ▼
┌────────────────────────┐
│Dimensionality Reduction│ (PCA -> UMAP)
└──────────┬─────────────┘
           │
           ▼
┌────────────────────────┐
│      Clustering        │ (Leiden Algorithm)
└──────────┬─────────────┘
           │
           ▼
┌────────────────────────┐
│ Find Marker Genes &    │
│  Cell Type Annotation  │
└──────────┬─────────────┘
           │
           ▼
┌────────────────────────┐
│   Annotated UMAP Plot  │
└────────────────────────┘
```

### Features

-   **Data Loading**: Uses `scanpy`'s built-in functions to load common scRNA-Seq data formats.
-   **Quality Control**: Implements standard QC, including filtering by gene counts, cell counts, and mitochondrial gene percentage.
-   **Dimensionality Reduction**: Utilizes Principal Component Analysis (PCA) and Uniform Manifold Approximation and Projection (UMAP).
-   **Graph-based Clustering**: Employs the highly performant Leiden algorithm.
-   **Marker Gene Identification**: Ranks genes by cluster to identify specific markers.
-   **Cell Type Annotation**: Provides a semi-automated approach to annotate clusters based on known marker genes.
-   **Visualization**: Generates key plots like UMAPs, dot plots, and violin plots to interpret results.

### Installation

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/olaflaitinen/scRNA-Seq-Clustering-Tool.git
    cd scRNA-Seq-Clustering-Tool
    ```

2.  **Create and activate a virtual environment (recommended):**
    ```bash
    python3 -m venv venv
    source venv/bin/activate
    ```

3.  **Install the required dependencies:**
    ```bash
    pip install -r requirements.txt
    ```

### Quick Start

The best way to get started is to run the guided walkthrough in the Jupyter Notebook.

```bash
jupyter-lab notebooks/analysis_walkthrough.ipynb
```

This notebook will load a public dataset (PBMC 3k from 10x Genomics) and take you through every step of the analysis, from loading data to generating the final annotated UMAP plot.

### License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
