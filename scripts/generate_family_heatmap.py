#!/usr/bin/env python3
"""
Generate heatmaps for verified gene families (CYP, OMT)

Reads a verified gene family TSV (from extract_gene_families.py) and the
original count matrix to produce expression heatmaps showing root vs leaf
patterns across replicates.

Single-species usage:
  python scripts/generate_family_heatmap.py \
      --genes 06_analysis/gene_families_DC/DC_CYP_verified.tsv \
      --counts 03_count_tables/00_1_DC/gene_count_matrix.tsv \
      --metadata 03_count_tables/00_1_DC/sample_metadata.tsv \
      --title "D. carota CYP Genes" \
      --output 06_analysis/gene_families_DC/cyp_heatmap.pdf

Combined two-species usage (DC + DG on one plot):
  python scripts/generate_family_heatmap.py \
      --genes 06_analysis/gene_families_DC_DG/DC_DG_CYP_verified.tsv \
      --counts 03_count_tables/00_1_DC/gene_count_matrix.tsv \
      --metadata 03_count_tables/00_1_DC/sample_metadata.tsv \
      --counts2 03_count_tables/00_2_DG/gene_count_matrix.tsv \
      --metadata2 03_count_tables/00_2_DG/sample_metadata.tsv \
      --species1 DC --species2 DG \
      --title "CYP Expression: DC vs DG" \
      --output 06_analysis/gene_families_DC_DG/cyp_heatmap_combined.pdf
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


BIOTYPE_PALETTE = {
    'protein_coding': '#e41a1c',
    'lncRNA':         '#377eb8',
    'pseudogene':     '#4daf4a',
    'miRNA':          '#984ea3',
    'snRNA':          '#ff7f00',
    'snoRNA':         '#ffff33',
    'rRNA':           '#a65628',
    'misc_RNA':       '#f781bf',
    'ncRNA':          '#377eb8',
}


def parse_biotype_from_gtf(gtf_path):
    """Extract gene_id -> gene_biotype from a GTF file."""
    if gtf_path is None or not Path(gtf_path).exists():
        return {}
    print(f"  Parsing gene biotype from GTF: {gtf_path}")
    biotype_map = {}
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] != 'gene':
                continue
            gid = bt = None
            for token in fields[8].split(';'):
                token = token.strip()
                if not token:
                    continue
                kv = token.split(' ', 1)
                if len(kv) != 2:
                    continue
                k, v = kv[0], kv[1].strip('"')
                if k == 'gene_id':
                    gid = v
                elif k in ('gene_biotype', 'gene_type'):
                    bt = v
            if gid and bt:
                biotype_map[gid] = bt
    print(f"    Found biotype for {len(biotype_map)} genes")
    return biotype_map


def _biotype_row_colors(gene_ids, biotype_map):
    """Build a list of colors from biotype_map for the given gene_ids."""
    if not biotype_map:
        return None
    return [BIOTYPE_PALETTE.get(biotype_map.get(gid, 'other'), '#999999')
            for gid in gene_ids]


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
                     condition_col='condition', sig_only=False,
                     biotype_map=None, top_n=None):
    """Generate a publication-quality heatmap.
    top_n: if set, keep only the top N genes ranked by padj then |log2FC|.
    """

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

    if top_n and len(gene_ids) > top_n:
        if 'padj' in genes_df.columns:
            ranked = genes_df.copy()
            ranked['_abs_lfc'] = ranked['log2FoldChange'].abs() if 'log2FoldChange' in ranked.columns else 0
            ranked = ranked.sort_values(['padj', '_abs_lfc'], ascending=[True, False])
            gene_ids = ranked['gene_id'].tolist()[:top_n]
        else:
            gene_ids = gene_ids[:top_n]
        print(f"  Limiting to top {top_n} genes (by padj, then |log2FC|)")

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

    gene_ids_present = list(present)
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

    # Row colors: gene family + biotype
    row_colors_parts = {}
    if color_by_family and len(set(family_colors)) > 1:
        row_colors_parts['Family'] = pd.Series(family_colors, index=heatmap_data.index)

    biotype_legend = []
    bt_colors = _biotype_row_colors(gene_ids_present, biotype_map or {})
    if bt_colors:
        row_colors_parts['Biotype'] = pd.Series(bt_colors, index=heatmap_data.index)
        seen = set()
        for gid, c in zip(gene_ids_present, bt_colors):
            bt = (biotype_map or {}).get(gid, 'other')
            if bt not in seen:
                seen.add(bt)
                biotype_legend.append(Patch(facecolor=c, label=bt))

    row_colors = pd.DataFrame(row_colors_parts) if row_colors_parts else None

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
        if biotype_legend:
            legend_elements.append(Patch(facecolor='none', edgecolor='none', label=''))
            legend_elements.extend(biotype_legend)

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


def _load_and_normalize(count_matrix_path, gene_ids):
    """Load a count matrix, normalize, log-transform, and subset to gene_ids."""
    counts = pd.read_csv(count_matrix_path, sep='\t', index_col=0)
    counts = counts[~counts.index.str.startswith('N_')]
    present = [g for g in gene_ids if g in counts.index]
    norm = normalize_counts(counts)
    log_norm = np.log2(norm.loc[present] + 1)
    return log_norm, present


def _order_samples_by_condition(meta, columns, condition_col):
    """Return sample names ordered root-first then leaf within the metadata."""
    sample_col = meta.columns[0]
    root_conds = {'R', 'root'}
    leaf_conds = {'L', 'leaf'}
    root_samples, leaf_samples, other_samples = [], [], []
    for _, row in meta.iterrows():
        s = row[sample_col]
        if s not in columns:
            continue
        cond = row.get(condition_col, '')
        if cond in root_conds:
            root_samples.append(s)
        elif cond in leaf_conds:
            leaf_samples.append(s)
        else:
            other_samples.append(s)
    return root_samples + leaf_samples + other_samples


def _build_row_labels(heatmap_data, genes_df, color_by_family):
    """Build row labels and family color list from gene info."""
    gene_info = genes_df.set_index('gene_id')
    family_palette = {'CYP': '#e74c3c', 'OMT': '#3498db'}
    row_labels = []
    family_colors = []
    for gid in heatmap_data.index:
        label = gid
        if gid in gene_info.index:
            info = gene_info.loc[gid]
            if isinstance(info, pd.DataFrame):
                info = info.iloc[0]
            desc = info.get('blast_description', '')
            if pd.notna(desc) and desc:
                label = f"{gid} | {str(desc)[:40]}"
            fam = info.get('gene_family', '')
            family_colors.append(family_palette.get(fam, '#95a5a6'))
        else:
            family_colors.append('#95a5a6')
        row_labels.append(label)
    return row_labels, family_colors, family_palette


def generate_combined_heatmap(genes_df, count_matrix1, metadata1,
                              count_matrix2, metadata2,
                              species1, species2, output_file,
                              title="Gene Family Heatmap", scale='center',
                              color_by_family=False, cluster_rows=True,
                              condition_col='condition', sig_only=False,
                              biotype_map=None, top_n=None):
    """Generate a single heatmap with two species side by side.

    Columns are ordered: SP1-Root | SP1-Leaf | SP2-Root | SP2-Leaf.
    Two stacked column color bars show species and tissue type.
    Each species is normalized independently before merging.
    top_n: if set, keep only the top N genes ranked by padj then |log2FC|.
    """

    print("=" * 60)
    print(f"Generating Combined Heatmap: {title}")
    print(f"  Species: {species1} + {species2}")
    print("=" * 60)

    gene_ids = genes_df['gene_id'].tolist()
    if sig_only and 'direction' in genes_df.columns:
        sig_mask = genes_df['direction'].isin(['root_up', 'leaf_up'])
        gene_ids = genes_df.loc[sig_mask, 'gene_id'].tolist()
        print(f"  Filtering to significant DE genes: {len(gene_ids)}")

    if top_n and len(gene_ids) > top_n:
        if 'padj' in genes_df.columns:
            ranked = genes_df.copy()
            ranked['_abs_lfc'] = ranked['log2FoldChange'].abs() if 'log2FoldChange' in ranked.columns else 0
            ranked = ranked.sort_values(['padj', '_abs_lfc'], ascending=[True, False])
            gene_ids = ranked['gene_id'].tolist()[:top_n]
        else:
            gene_ids = gene_ids[:top_n]
        print(f"  Limiting to top {top_n} genes (by padj, then |log2FC|)")

    if not gene_ids:
        print("  No genes to plot!")
        return

    # Normalize each species independently
    print(f"  Loading and normalizing {species1} counts...")
    log1, present1 = _load_and_normalize(count_matrix1, gene_ids)
    print(f"    {species1}: {len(present1)} genes found")

    print(f"  Loading and normalizing {species2} counts...")
    log2, present2 = _load_and_normalize(count_matrix2, gene_ids)
    print(f"    {species2}: {len(present2)} genes found")

    # Use genes present in both species
    shared_genes = [g for g in gene_ids if g in log1.index and g in log2.index]
    if not shared_genes:
        print("  ERROR: No genes found in both count matrices")
        return
    print(f"  Shared genes in both species: {len(shared_genes)}")

    log1 = log1.loc[shared_genes]
    log2 = log2.loc[shared_genes]

    # Order samples within each species: root then leaf
    meta1 = pd.read_csv(metadata1, sep='\t')
    meta2 = pd.read_csv(metadata2, sep='\t')

    sp1_ordered = _order_samples_by_condition(meta1, set(log1.columns), condition_col)
    sp2_ordered = _order_samples_by_condition(meta2, set(log2.columns), condition_col)

    # Concatenate: SP1 samples | SP2 samples
    heatmap_data = pd.concat([log1[sp1_ordered], log2[sp2_ordered]], axis=1)

    # Center per gene across all samples
    if scale == 'center':
        row_means = heatmap_data.mean(axis=1)
        heatmap_data = heatmap_data.subtract(row_means, axis=0)
        cbar_label = "Log2 Expression (centered)"
    elif scale == 'zscore':
        row_means = heatmap_data.mean(axis=1)
        row_stds = heatmap_data.std(axis=1).replace(0, 1)
        heatmap_data = heatmap_data.subtract(row_means, axis=0).div(row_stds, axis=0)
        cbar_label = "z-score"
    else:
        cbar_label = "log2(norm+1)"

    # Build row labels
    row_labels, family_colors, family_palette = _build_row_labels(
        heatmap_data, genes_df, color_by_family
    )
    heatmap_data.index = row_labels

    # Build two-row column color bars: Species (top) and Tissue Type (bottom)
    sp1_set = set(sp1_ordered)
    sp2_set = set(sp2_ordered)

    species_palette = {species1: '#8e44ad', species2: '#2980b9'}
    tissue_palette = {'Root': '#d35400', 'Leaf': '#27ae60'}

    sp1_sample_col = meta1.columns[0]
    sp2_sample_col = meta2.columns[0]

    root_conds = {'R', 'root'}
    leaf_conds = {'L', 'leaf'}

    species_colors = []
    tissue_colors = []
    for sample in heatmap_data.columns:
        if sample in sp1_set:
            species_colors.append(species_palette[species1])
            match = meta1[meta1[sp1_sample_col] == sample]
            if not match.empty:
                cond = match[condition_col].values[0]
            else:
                cond = ''
        else:
            species_colors.append(species_palette[species2])
            match = meta2[meta2[sp2_sample_col] == sample]
            if not match.empty:
                cond = match[condition_col].values[0]
            else:
                cond = ''

        if cond in root_conds:
            tissue_colors.append(tissue_palette['Root'])
        elif cond in leaf_conds:
            tissue_colors.append(tissue_palette['Leaf'])
        else:
            tissue_colors.append('#bdc3c7')

    col_colors_df = pd.DataFrame({
        'Species': species_colors,
        'Tissue Type': tissue_colors,
    }, index=heatmap_data.columns)

    # Figure sizing
    n_genes = len(heatmap_data)
    n_samples = len(heatmap_data.columns)
    fig_height = max(6, 0.35 * n_genes + 2)
    fig_width = max(10, 0.8 * n_samples + 4)

    row_colors_parts = {}
    if color_by_family and len(set(family_colors)) > 1:
        row_colors_parts['Family'] = pd.Series(family_colors, index=heatmap_data.index)

    biotype_legend = []
    bt_colors = _biotype_row_colors(shared_genes, biotype_map or {})
    if bt_colors:
        row_colors_parts['Biotype'] = pd.Series(bt_colors, index=heatmap_data.index)
        seen = set()
        for gid, c in zip(shared_genes, bt_colors):
            bt = (biotype_map or {}).get(gid, 'other')
            if bt not in seen:
                seen.add(bt)
                biotype_legend.append(Patch(facecolor=c, label=bt))

    row_colors = pd.DataFrame(row_colors_parts) if row_colors_parts else None

    print(f"  Plotting {n_genes} genes x {n_samples} samples...")

    try:
        g = sns.clustermap(
            heatmap_data,
            cmap='RdBu_r',
            center=0,
            figsize=(fig_width, fig_height),
            row_cluster=cluster_rows and n_genes > 2,
            col_cluster=False,
            col_colors=col_colors_df,
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

        legend_elements = [
            Patch(facecolor=species_palette[species1], label=species1),
            Patch(facecolor=species_palette[species2], label=species2),
            Patch(facecolor='none', edgecolor='none', label=''),
            Patch(facecolor=tissue_palette['Root'], label='Root'),
            Patch(facecolor=tissue_palette['Leaf'], label='Leaf'),
        ]
        if color_by_family and len(set(family_colors)) > 1:
            legend_elements.append(Patch(facecolor='none', edgecolor='none', label=''))
            for fam, color in family_palette.items():
                if 'gene_family' in genes_df.columns and any(genes_df['gene_family'] == fam):
                    legend_elements.append(Patch(facecolor=color, label=fam))
        if biotype_legend:
            legend_elements.append(Patch(facecolor='none', edgecolor='none', label=''))
            legend_elements.extend(biotype_legend)

        g.ax_heatmap.legend(
            handles=legend_elements,
            loc='upper left',
            bbox_to_anchor=(1.15, 1.0),
            frameon=True,
            fontsize=9,
        )

        out_path = Path(output_file)
        g.savefig(out_path, dpi=300, bbox_inches='tight')
        png_path = out_path.with_suffix('.png')
        g.savefig(png_path, dpi=150, bbox_inches='tight')
        plt.close()

        print(f"  Saved: {out_path}")
        print(f"  Saved: {png_path}")

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
                        help='Gene count matrix TSV (species 1)')
    parser.add_argument('--metadata', required=True,
                        help='Sample metadata TSV (species 1)')
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
    parser.add_argument('--gtf', default=None,
                        help='GTF annotation file (optional, adds gene biotype sidebar)')

    # Combined two-species mode
    parser.add_argument('--counts2',
                        help='Second species count matrix (enables combined mode)')
    parser.add_argument('--metadata2',
                        help='Second species metadata TSV')
    parser.add_argument('--species1', default='SP1',
                        help='Label for species 1 (default: SP1)')
    parser.add_argument('--species2', default='SP2',
                        help='Label for species 2 (default: SP2)')
    parser.add_argument('--top-n', type=int, default=None,
                        help='Only plot the top N genes per family, ranked by padj then |log2FC| (default: all)')

    args = parser.parse_args()

    # Load gene list
    genes_df = pd.read_csv(args.genes, sep='\t', comment='#')
    print(f"Loaded {len(genes_df)} genes from {args.genes}")

    biotype_map = parse_biotype_from_gtf(args.gtf)

    if args.counts2:
        if not args.metadata2:
            print("ERROR: --metadata2 is required when using --counts2")
            return 1
        generate_combined_heatmap(
            genes_df=genes_df,
            count_matrix1=args.counts,
            metadata1=args.metadata,
            count_matrix2=args.counts2,
            metadata2=args.metadata2,
            species1=args.species1,
            species2=args.species2,
            output_file=args.output,
            title=args.title,
            scale=args.scale,
            color_by_family=args.color_by_family,
            cluster_rows=not args.no_cluster,
            condition_col=args.condition_col,
            sig_only=args.sig_only,
            biotype_map=biotype_map,
            top_n=args.top_n,
        )
    else:
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
            biotype_map=biotype_map,
            top_n=args.top_n,
        )

    print("\nDone!")
    return 0


if __name__ == '__main__':
    sys.exit(main())
