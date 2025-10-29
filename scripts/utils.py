"""
Utility functions for scRNA-Seq analysis workflow.

This module provides helper functions for loading configuration,
performing common preprocessing tasks, and managing analysis parameters.
"""

import yaml
from pathlib import Path
from typing import Dict, Any, Optional, List
import scanpy as sc
import pandas as pd


class AnalysisConfig:
    """
    Class to manage analysis configuration parameters.

    This class loads and provides access to all configuration parameters
    defined in config.yaml.
    """

    def __init__(self, config_path: str = "config.yaml"):
        """
        Initialize configuration.

        Args:
            config_path: Path to the YAML configuration file
        """
        self.config_path = Path(config_path)
        self.config = self._load_config()

    def _load_config(self) -> Dict[str, Any]:
        """
        Load configuration from YAML file.

        Returns:
            Dictionary containing configuration parameters
        """
        if not self.config_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {self.config_path}")

        with open(self.config_path, 'r') as f:
            config = yaml.safe_load(f)

        return config

    def get(self, key: str, default: Any = None) -> Any:
        """
        Get configuration value by key.

        Args:
            key: Configuration key (supports dot notation, e.g., 'qc.min_genes')
            default: Default value if key not found

        Returns:
            Configuration value
        """
        keys = key.split('.')
        value = self.config

        for k in keys:
            if isinstance(value, dict) and k in value:
                value = value[k]
            else:
                return default

        return value

    def get_paths(self) -> Dict[str, Path]:
        """
        Get all path configurations as Path objects.

        Returns:
            Dictionary of paths
        """
        paths = self.config.get('paths', {})
        return {key: Path(value) for key, value in paths.items()}

    def get_qc_params(self) -> Dict[str, Any]:
        """Get quality control parameters."""
        return self.config.get('qc', {})

    def get_marker_genes(self) -> Dict[str, List[str]]:
        """Get marker genes for cell type annotation."""
        return self.config.get('annotation', {}).get('marker_genes', {})

    def get_cluster_mapping(self) -> Dict[str, str]:
        """Get cluster to cell type mapping."""
        return self.config.get('annotation', {}).get('cluster_mapping', {})


def setup_scanpy(config: AnalysisConfig) -> None:
    """
    Configure scanpy settings based on configuration.

    Args:
        config: AnalysisConfig object
    """
    verbosity = config.get('computation.verbosity', 3)
    sc.settings.verbosity = verbosity

    dpi = config.get('visualization.dpi', 80)
    sc.settings.set_figure_params(dpi=dpi, facecolor='white', frameon=False)

    figures_dir = config.get('paths.figures_dir', 'figures/')
    sc.settings.figdir = Path(figures_dir)

    sc.logging.print_header()


def apply_qc_filters(adata: sc.AnnData, config: AnalysisConfig) -> sc.AnnData:
    """
    Apply quality control filters to AnnData object.

    Args:
        adata: AnnData object
        config: AnalysisConfig object

    Returns:
        Filtered AnnData object
    """
    qc_params = config.get_qc_params()

    print(f"Cells before filtering: {adata.n_obs}")
    print(f"Genes before filtering: {adata.n_vars}")

    # Filter cells and genes
    sc.pp.filter_cells(adata, min_genes=qc_params.get('min_genes', 200))
    sc.pp.filter_genes(adata, min_cells=qc_params.get('min_cells', 3))

    # Filter based on QC metrics
    max_genes = qc_params.get('max_genes', 2500)
    max_pct_mt = qc_params.get('max_pct_mt', 5)

    adata = adata[adata.obs.n_genes_by_counts < max_genes, :]
    adata = adata[adata.obs.pct_counts_mt < max_pct_mt, :]

    print(f"Cells after filtering: {adata.n_obs}")
    print(f"Genes after filtering: {adata.n_vars}")

    return adata


def preprocess_data(adata: sc.AnnData, config: AnalysisConfig) -> sc.AnnData:
    """
    Apply standard preprocessing pipeline.

    Args:
        adata: AnnData object
        config: AnalysisConfig object

    Returns:
        Preprocessed AnnData object
    """
    # Save raw counts
    if config.get('normalization.save_raw', True):
        adata.layers['counts'] = adata.X.copy()

    # Normalize
    target_sum = config.get('normalization.target_sum', 1e4)
    sc.pp.normalize_total(adata, target_sum=target_sum)

    # Log transform
    sc.pp.log1p(adata)

    # Save normalized data
    adata.raw = adata

    # Identify highly variable genes
    n_top_genes = config.get('feature_selection.n_top_genes', 2000)
    flavor = config.get('feature_selection.flavor', 'seurat')
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor=flavor, subset=False)

    # Subset to HVG if specified
    if config.get('feature_selection.subset', True):
        adata = adata[:, adata.var.highly_variable]

    # Regress out unwanted variation
    regress_vars = config.get('preprocessing.regress_out', [])
    if regress_vars:
        sc.pp.regress_out(adata, regress_vars)

    # Scale
    max_value = config.get('preprocessing.scale_max_value', 10)
    sc.pp.scale(adata, max_value=max_value)

    return adata


def run_dimensionality_reduction(adata: sc.AnnData, config: AnalysisConfig) -> None:
    """
    Perform PCA, neighborhood graph, and UMAP.

    Args:
        adata: AnnData object
        config: AnalysisConfig object
    """
    # PCA
    n_comps = config.get('pca.n_comps', 50)
    svd_solver = config.get('pca.svd_solver', 'arpack')
    sc.tl.pca(adata, svd_solver=svd_solver, n_comps=n_comps)

    # Neighborhood graph
    n_neighbors = config.get('neighbors.n_neighbors', 10)
    n_pcs = config.get('neighbors.n_pcs', 40)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

    # UMAP
    min_dist = config.get('umap.min_dist', 0.5)
    random_state = config.get('umap.random_state', 42)
    sc.tl.umap(adata, min_dist=min_dist, random_state=random_state)


def perform_clustering(adata: sc.AnnData, config: AnalysisConfig) -> None:
    """
    Perform clustering using specified algorithm.

    Args:
        adata: AnnData object
        config: AnalysisConfig object
    """
    algorithm = config.get('clustering.algorithm', 'leiden')
    resolution = config.get('clustering.resolution', 0.9)
    random_state = config.get('clustering.random_state', 42)

    if algorithm == 'leiden':
        sc.tl.leiden(adata, resolution=resolution, random_state=random_state)
    elif algorithm == 'louvain':
        sc.tl.louvain(adata, resolution=resolution, random_state=random_state)
    else:
        raise ValueError(f"Unknown clustering algorithm: {algorithm}")

    print(f"Number of clusters: {len(adata.obs[algorithm].unique())}")


def find_marker_genes(adata: sc.AnnData, config: AnalysisConfig) -> None:
    """
    Identify marker genes for each cluster.

    Args:
        adata: AnnData object
        config: AnalysisConfig object
    """
    method = config.get('marker_genes.method', 'wilcoxon')
    groupby = config.get('marker_genes.groupby', 'leiden')

    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        method=method,
        key_added=f'rank_genes_{groupby}'
    )


def annotate_cell_types(adata: sc.AnnData, config: AnalysisConfig,
                        cluster_col: str = 'leiden') -> None:
    """
    Annotate cell types based on cluster mapping.

    Args:
        adata: AnnData object
        config: AnalysisConfig object
        cluster_col: Name of the cluster column
    """
    cluster_mapping = config.get_cluster_mapping()

    # Convert cluster mapping keys to strings
    cluster_mapping = {str(k): v for k, v in cluster_mapping.items()}

    adata.obs['cell_type'] = adata.obs[cluster_col].map(cluster_mapping).astype('category')

    print("\nCell type distribution:")
    print(adata.obs['cell_type'].value_counts())


def export_results(adata: sc.AnnData, config: AnalysisConfig) -> None:
    """
    Export analysis results.

    Args:
        adata: AnnData object
        config: AnalysisConfig object
    """
    export_config = config.get('export', {})
    paths = config.get_paths()
    figures_dir = paths.get('figures_dir', Path('figures'))
    results_dir = paths.get('results_dir', Path('data'))
    table_format = export_config.get('table_format', 'csv')

    # Save AnnData object
    if export_config.get('save_adata', True):
        output_file = results_dir / config.get('paths.output_file', 'processed.h5ad')
        adata.write(output_file)
        print(f"Saved processed data to: {output_file}")

    # Export annotations
    if export_config.get('export_annotations', True):
        annotations_file = figures_dir / f'cell_annotations.{table_format}'
        annotations_df = adata.obs[['leiden', 'cell_type']].copy()

        if table_format == 'csv':
            annotations_df.to_csv(annotations_file)
        elif table_format == 'tsv':
            annotations_df.to_csv(annotations_file, sep='\t')
        elif table_format == 'excel':
            annotations_df.to_excel(annotations_file.with_suffix('.xlsx'))

        print(f"Exported annotations to: {annotations_file}")

    # Export marker genes
    if export_config.get('export_markers', True):
        markers_file = figures_dir / f'marker_genes.{table_format}'
        markers_df = pd.DataFrame(adata.uns['rank_genes_leiden']['names'])

        if table_format == 'csv':
            markers_df.to_csv(markers_file, index=False)
        elif table_format == 'tsv':
            markers_df.to_csv(markers_file, sep='\t', index=False)
        elif table_format == 'excel':
            markers_df.to_excel(markers_file.with_suffix('.xlsx'), index=False)

        print(f"Exported marker genes to: {markers_file}")


def get_all_marker_genes(config: AnalysisConfig) -> List[str]:
    """
    Get flattened list of all marker genes from config.

    Args:
        config: AnalysisConfig object

    Returns:
        List of unique marker genes
    """
    marker_dict = config.get_marker_genes()
    all_markers = [gene for genes in marker_dict.values() for gene in genes]

    # Remove duplicates while preserving order
    return list(dict.fromkeys(all_markers))
