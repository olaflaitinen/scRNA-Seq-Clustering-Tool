#!/usr/bin/env python3
"""
Generate Sample UMAP Visualizations

This script creates sample UMAP plots for documentation purposes.
These represent typical outputs from the scRNA-Seq analysis pipeline.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

# Set random seed for reproducibility
np.random.seed(42)

# Configure matplotlib
plt.style.use('default')
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['legend.fontsize'] = 9

# Output directory
OUTPUT_DIR = Path('assets')
OUTPUT_DIR.mkdir(exist_ok=True)


def generate_cluster_centers():
    """Generate realistic cluster centers for PBMC data."""
    # Define 9 cluster centers in 2D space (UMAP coordinates)
    centers = np.array([
        [-5, 8],    # Cluster 0: CD4 T cells (large)
        [8, 5],     # Cluster 1: CD14+ Monocytes
        [-8, -5],   # Cluster 2: B cells
        [-3, 5],    # Cluster 3: CD4 T cells (another group)
        [5, -8],    # Cluster 4: CD8 T cells
        [10, -2],   # Cluster 5: NK cells
        [12, 8],    # Cluster 6: FCGR3A+ Monocytes
        [-10, 2],   # Cluster 7: Dendritic cells
        [0, -12],   # Cluster 8: Megakaryocytes
    ])
    return centers


def generate_cluster_data(centers, n_cells_per_cluster):
    """Generate synthetic UMAP coordinates for each cluster."""
    all_points = []
    all_labels = []

    for cluster_id, center in enumerate(centers):
        n_cells = n_cells_per_cluster[cluster_id]

        # Generate points around cluster center with varying spread
        spread = np.random.uniform(0.8, 1.5)
        points = np.random.randn(n_cells, 2) * spread + center

        all_points.append(points)
        all_labels.extend([cluster_id] * n_cells)

    return np.vstack(all_points), np.array(all_labels)


def plot_umap_clusters(coords, labels, output_path):
    """Create UMAP plot colored by cluster."""
    fig, ax = plt.subplots(figsize=(10, 8))

    # Define color palette (similar to scanpy default)
    colors = [
        '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
        '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22'
    ]

    # Plot each cluster
    for cluster_id in np.unique(labels):
        mask = labels == cluster_id
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            c=colors[cluster_id],
            label=f'{cluster_id}',
            s=5,
            alpha=0.7,
            edgecolors='none'
        )

    # Add cluster labels on the plot
    for cluster_id in np.unique(labels):
        mask = labels == cluster_id
        center = coords[mask].mean(axis=0)
        ax.text(
            center[0], center[1],
            str(cluster_id),
            fontsize=16,
            fontweight='bold',
            ha='center',
            va='center',
            color='white',
            bbox=dict(boxstyle='circle', facecolor=colors[cluster_id],
                     edgecolor='white', linewidth=2, alpha=0.9)
        )

    ax.set_xlabel('UMAP 1', fontsize=12)
    ax.set_ylabel('UMAP 2', fontsize=12)
    ax.set_title('UMAP by Leiden Cluster', fontsize=14, fontweight='bold', pad=15)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(False)

    # Legend
    handles = [mpatches.Patch(color=colors[i], label=f'Cluster {i}')
               for i in range(len(np.unique(labels)))]
    ax.legend(
        handles=handles,
        title='Leiden',
        loc='center left',
        bbox_to_anchor=(1, 0.5),
        frameon=True,
        fancybox=False,
        shadow=False
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"✓ Generated: {output_path}")


def plot_umap_annotated(coords, labels, output_path):
    """Create UMAP plot colored by cell type."""
    fig, ax = plt.subplots(figsize=(10, 8))

    # Map cluster IDs to cell types
    cluster_to_celltype = {
        0: 'CD4 T cells',
        1: 'CD14+ Monocytes',
        2: 'B cells',
        3: 'CD4 T cells',
        4: 'CD8 T cells',
        5: 'NK cells',
        6: 'FCGR3A+ Monocytes',
        7: 'Dendritic cells',
        8: 'Megakaryocytes'
    }

    # Define colors for cell types
    celltype_colors = {
        'CD4 T cells': '#1f77b4',
        'CD14+ Monocytes': '#ff7f0e',
        'B cells': '#2ca02c',
        'CD8 T cells': '#d62728',
        'NK cells': '#9467bd',
        'FCGR3A+ Monocytes': '#8c564b',
        'Dendritic cells': '#e377c2',
        'Megakaryocytes': '#bcbd22'
    }

    # Create cell type labels
    celltypes = [cluster_to_celltype[label] for label in labels]

    # Plot each cell type
    for celltype, color in celltype_colors.items():
        mask = np.array([ct == celltype for ct in celltypes])
        if mask.any():
            ax.scatter(
                coords[mask, 0],
                coords[mask, 1],
                c=color,
                label=celltype,
                s=5,
                alpha=0.7,
                edgecolors='none'
            )

    ax.set_xlabel('UMAP 1', fontsize=12)
    ax.set_ylabel('UMAP 2', fontsize=12)
    ax.set_title('UMAP by Annotated Cell Type', fontsize=14, fontweight='bold', pad=15)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(False)

    # Legend
    handles = [mpatches.Patch(color=celltype_colors[ct], label=ct)
               for ct in celltype_colors.keys()]
    ax.legend(
        handles=handles,
        title='Cell Type',
        loc='center left',
        bbox_to_anchor=(1, 0.5),
        frameon=True,
        fancybox=False,
        shadow=False
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"✓ Generated: {output_path}")


def main():
    """Main function to generate sample UMAP plots."""
    print("Generating sample UMAP visualizations...")
    print()

    # Define cluster centers
    centers = generate_cluster_centers()

    # Define number of cells per cluster (mimicking typical PBMC distribution)
    n_cells_per_cluster = [600, 400, 250, 300, 200, 150, 100, 80, 50]  # Total ~2130 cells

    # Generate synthetic UMAP coordinates
    coords, labels = generate_cluster_data(centers, n_cells_per_cluster)

    print(f"Generated {len(coords)} synthetic cells across {len(centers)} clusters")
    print()

    # Generate UMAP by cluster
    cluster_output = OUTPUT_DIR / 'umap_clusters.png'
    plot_umap_clusters(coords, labels, cluster_output)

    # Generate UMAP by cell type
    annotated_output = OUTPUT_DIR / 'umap_annotated.png'
    plot_umap_annotated(coords, labels, annotated_output)

    print()
    print("✓ All visualizations generated successfully!")
    print(f"  Output directory: {OUTPUT_DIR.absolute()}")


if __name__ == "__main__":
    main()
