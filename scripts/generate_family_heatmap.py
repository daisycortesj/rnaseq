#!/usr/bin/env python3
"""
Generate heatmaps for verified gene families (CYP, OMT)

Reads a verified gene family TSV (from extract_gene_families.py) and the
original count matrix to produce expression heatmaps showing root vs leaf
patterns across replicates.

Usage:
  python scripts/generate_family_heatmap.py \
      --genes 06_analysis/gene_families_DC/DC_CYP_verified.tsv \
      --counts 03_count_tables/00_1_DC/gene_count_matrix.tsv \
      --metadata 03_count_tables/00_1_DC/sample_metadata.tsv \
      --title "D. carota CYP Genes" \
      --output 06_analysis/gene_families_DC/cyp_heatmap.pdf

  # Combined CYP + OMT heatmap:
  python scripts/generate_family_heatmap.py \
      --genes 06_analysis/gene_families_DC/DC_CYP_OMT_combined.tsv \
      --counts 03_count_tables/00_1_DC/gene_count_matrix.tsv \
      --metadata 03_count_tables/00_1_DC/sample_metadata.tsv \
      --title "D. carota CYP + OMT Genes" \
      --output 06_analysis/gene_families_DC/cyp_omt_heatmap.pdf \
      --color-by-family
"""

import pandas as pd
import numpy as np
import argparse
import sys
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.patches import Patch
except ImportError as e:
    print(f"ERROR: Missing required package: {e}")
    print("Install with: pip install matplotlib seaborn")
    sys.exit(1)


def normalize_counts(count_matrix):
    """DESeq2-style median-of-ratios normalization."""
    geo_means = count_matrix.apply(
        lambda row: np.exp(np.log(row[row > 0]).mean()) if (row > 0).sum() > 0 else 0,
        axis=1
    )

    size_factors = []
    for col in count_matrix.columns:
        ratios = count_matrix[col] / geo_means
        ratios = ratios[np.isfinite(ratios) & (ratios > 0)]
        size_factors.append(ratios.median())

    size_factors = pd.Series(size_factors, index=count_matrix.columns)
    return count_matrix.div(size_factors, axis=1)


def generate_heatmap(genes_df, count_matrix, metadata, output_file,
                     title="Gene Family Heatmap", scale='center',
                     color_by_family=False, cluster_rows=True,
                     condition_col='condition', sig_only=False):
    """Generate a publication-quality heatmap."""

    print("=" * 60)
    print(f"Generating Heatmap: {title}")
    print("=" * 60)

    # Get gene IDs from verified list
    gene_ids = genes_df['gene_id'].tolist()

    # Optionally filter to significant DE genes
    if sig_only and 'direction' in genes_df.columns:
        sig_mask = genes_df['direction'].isin(['root_up', 'leaf_up'])
        gene_ids = genes_df.loc[sig_mask, 'gene_id'].tolist()
        print(f"  Filtering to significant DE genes: {len(gene_ids)}")

    if not gene_ids:
        print("  No genes to plot!")
        return

    # Load and normalize counts
    print(f"  Loading count matrix...")
    counts = pd.read_csv(count_matrix, sep='\t', index_col=0)

    # Drop STAR summary rows
    counts = counts[~counts.index.str.startswith('N_')]

    present = [g for g in gene_ids if g in counts.index]
    missing = [g for g in gene_ids if g not in counts.index]
    if missing:
        print(f"  WARNING: {len(missing)} genes not found in count matrix")
    if not present:
        print("  ERROR: No genes found in count matrix")
        return

    print(f"  Normalizing counts ({len(present)} genes)...")
    norm_counts = normalize_counts(counts)
    heatmap_data = np.log2(norm_counts.loc[present] + 1)

    # Scale per gene
    if scale == 'center':
        row_means = heatmap_data.mean(axis=1)
        heatmap_data = heatmap_data.subtract(row_means, axis=0)
        cbar_label = "log2(norm+1), centered"
    elif scale == 'zscore':
        row_means = heatmap_data.mean(axis=1)
        row_stds = heatmap_data.std(axis=1).replace(0, 1)
        heatmap_data = heatmap_data.subtract(row_means, axis=0).div(row_stds, axis=0)
        cbar_label = "z-score"
    else:
        cbar_label = "log2(norm+1)"

    # Order columns by condition (root first, then leaf)
    meta = pd.read_csv(metadata, sep='\t')
    if condition_col in meta.columns:
        sample_col = meta.columns[0]
        meta = meta.sort_values(condition_col)
        ordered_samples = [s for s in meta[sample_col] if s in heatmap_data.columns]
        heatmap_data = heatmap_data[ordered_samples]

    # Build row labels with gene family + description
    row_labels = []
    gene_info = genes_df.set_index('gene_id')
    family_colors = []
    family_palette = {'CYP': '#e74c3c', 'OMT': '#3498db'}

    for gid in heatmap_data.index:
        label = gid
        if gid in gene_info.index:
            info = gene_info.loc[gid]
            if isinstance(info, pd.DataFrame):
                info = info.iloc[0]
            desc = info.get('blast_description', '')
            if pd.notna(desc) and desc:
                short_desc = str(desc)[:40]
                label = f"{gid} | {short_desc}"
            fam = info.get('gene_family', '')
            family_colors.append(family_palette.get(fam, '#95a5a6'))
        else:
            family_colors.append('#95a5a6')
        row_labels.append(label)

    heatmap_data.index = row_labels

    # Column colors for condition
    col_colors = None
    if condition_col in meta.columns:
        cond_palette = {'R': '#d35400', 'L': '#27ae60',
                        'root': '#d35400', 'leaf': '#27ae60'}
        col_colors_list = []
        for sample in heatmap_data.columns:
            match = meta[meta[sample_col] == sample]
            if not match.empty:
                cond = match[condition_col].values[0]
                col_colors_list.append(cond_palette.get(cond, '#bdc3c7'))
            else:
                col_colors_list.append('#bdc3c7')
        col_colors = pd.Series(col_colors_list, index=heatmap_data.columns)

    # Figure sizing
    n_genes = len(heatmap_data)
    n_samples = len(heatmap_data.columns)
    fig_height = max(6, 0.35 * n_genes + 2)
    fig_width = max(8, 0.8 * n_samples + 4)

    # Row colors for gene family
    row_colors = None
    if color_by_family and len(set(family_colors)) > 1:
        row_colors = pd.Series(family_colors, index=heatmap_data.index)

    print(f"  Plotting {n_genes} genes x {n_samples} samples...")

    try:
        g = sns.clustermap(
            heatmap_data,
            cmap='RdBu_r',
            center=0,
            figsize=(fig_width, fig_height),
            row_cluster=cluster_rows and n_genes > 2,
            col_cluster=False,
            col_colors=col_colors,
            row_colors=row_colors,
            linewidths=0.5,
            linecolor='white',
            cbar_kws={'label': cbar_label, 'shrink': 0.5},
            yticklabels=True,
            xticklabels=True,
            dendrogram_ratio=(0.1, 0.05),
        )

        g.ax_heatmap.set_ylabel('')
        g.ax_heatmap.set_xlabel('Samples', fontsize=11)
        g.fig.suptitle(title, fontsize=14, fontweight='bold', y=1.02)

        # Add legend for conditions
        legend_elements = []
        if condition_col in meta.columns:
            legend_elements.extend([
                Patch(facecolor='#d35400', label='Root (R)'),
                Patch(facecolor='#27ae60', label='Leaf (L)'),
            ])
        if color_by_family and len(set(family_colors)) > 1:
            for fam, color in family_palette.items():
                if any(genes_df['gene_family'] == fam):
                    legend_elements.append(
                        Patch(facecolor=color, label=fam)
                    )

        if legend_elements:
            g.ax_heatmap.legend(
                handles=legend_elements,
                loc='upper left',
                bbox_to_anchor=(1.15, 1.0),
                frameon=True,
                fontsize=9
            )

        # Save
        out_path = Path(output_file)
        g.savefig(out_path, dpi=300, bbox_inches='tight')
        png_path = out_path.with_suffix('.png')
        g.savefig(png_path, dpi=150, bbox_inches='tight')
        plt.close()

        print(f"  Saved: {out_path}")
        print(f"  Saved: {png_path}")

        # Save underlying matrix
        matrix_path = out_path.with_suffix('.tsv')
        heatmap_data.to_csv(matrix_path, sep='\t')
        print(f"  Saved matrix: {matrix_path}")

    except Exception as e:
        print(f"  ERROR generating heatmap: {e}")
        import traceback
        traceback.print_exc()


def main():
    parser = argparse.ArgumentParser(
        description='Generate heatmaps for verified gene families',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('--genes', required=True,
                        help='Verified gene family TSV (from extract_gene_families.py)')
    parser.add_argument('--counts', required=True,
                        help='Gene count matrix TSV')
    parser.add_argument('--metadata', required=True,
                        help='Sample metadata TSV')
    parser.add_argument('-o', '--output', required=True,
                        help='Output PDF path')
    parser.add_argument('--title', default='Gene Family Heatmap',
                        help='Plot title')
    parser.add_argument('--scale', choices=['center', 'zscore', 'none'],
                        default='center',
                        help='Scaling method (default: center)')
    parser.add_argument('--color-by-family', action='store_true',
                        help='Color row sidebar by gene family (CYP vs OMT)')
    parser.add_argument('--sig-only', action='store_true',
                        help='Only plot significantly DE genes')
    parser.add_argument('--no-cluster', action='store_true',
                        help='Disable row clustering')
    parser.add_argument('--condition-col', default='condition',
                        help='Metadata column for condition (default: condition)')

    args = parser.parse_args()

    # Load gene list
    genes_df = pd.read_csv(args.genes, sep='\t', comment='#')
    print(f"Loaded {len(genes_df)} genes from {args.genes}")

    generate_heatmap(
        genes_df=genes_df,
        count_matrix=args.counts,
        metadata=args.metadata,
        output_file=args.output,
        title=args.title,
        scale=args.scale,
        color_by_family=args.color_by_family,
        cluster_rows=not args.no_cluster,
        condition_col=args.condition_col,
        sig_only=args.sig_only,
    )

    print("\nDone!")
    return 0


if __name__ == '__main__':
    sys.exit(main())
