#!/usr/bin/env python3
"""
PyDESeq2 Step 3: Generate Plots
================================

Generates publication-quality plots from DESeq2 results:
  - MA plot (baseMean vs log2FC, blue significant points)
  - Enhanced Volcano plot (4-color: NS / Log2FC / padj / both)
  - PCA plot (top 50 variable genes, PC1 vs PC2)
  - Sample correlation heatmap (Euclidean distance, blue gradient)
  - CYP heatmap (auto-detected from BLAST annotations or gene_family column)
  - OMT heatmap (auto-detected from BLAST annotations or gene_family column)

Accepts any of these as input:
  - Unfiltered DESeq2 results (from Step 1)
  - Filtered DESeq2 results (from Step 2)
  - Combined BLAST+DESeq2 annotated file (from run_combine_filter.sbatch)
  - Gene-family-verified file (from extract_gene_families.py)

CYP and OMT genes are identified automatically:
  - If the input has a 'gene_family' column, that is used directly
  - If the input has a 'blast_description' column, regex patterns identify
    CYP and OMT genes (same patterns as extract_gene_families.py)
  - If neither column exists, gene-family heatmaps are skipped

Usage:
  # From combined annotated file (auto-detects CYP/OMT):
  python pydeseq2_generate_plots.py combined_annotated.tsv \\
      --count-matrix count_matrix.tsv --metadata metadata.tsv -o plots/

  # With HMMER domain annotations on heatmap rows:
  python pydeseq2_generate_plots.py combined_annotated.tsv \\
      --count-matrix count_matrix.tsv --metadata metadata.tsv \\
      --hmmer pfam_domains.txt -o plots/

  # Combined two-species heatmap (DC + DG on one plot):
  python pydeseq2_generate_plots.py DC_annotated.tsv \\
      --count-matrix DC_counts.tsv --metadata DC_metadata.tsv \\
      --count-matrix2 DG_counts.tsv --metadata2 DG_metadata.tsv \\
      --species1 DC --species2 DG -o plots/

  # From raw DESeq2 results (MA + volcano only, no heatmaps):
  python pydeseq2_generate_plots.py pydeseq2_results.tsv -o plots/
"""

import os
import re
import sys
import argparse
import pandas as pd
import numpy as np
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
    print("  pip install matplotlib seaborn")
    sys.exit(1)


# ============================================================================
# GENE FAMILY PATTERNS (same as extract_gene_families.py)
# ============================================================================

GENE_FAMILIES = {
    'CYP': {
        'full_name': 'Cytochrome P450',
        'blast_patterns': [
            r'cytochrome\s+P[- ]?450',
            r'\bCYP\d',
            r'\bP450\b',
            r'monooxygenase',
        ],
        'blast_exclude': [
            r'NADPH.cytochrome\s+P450\s+reductase',
            r'cytochrome\s+P450\s+reductase',
            r'cytochrome\s+b5',
        ],
        'hmmer_pfam': ['PF00067'],
        'hmmer_names': ['p450'],
    },
    'OMT': {
        'full_name': 'O-Methyltransferase',
        'blast_patterns': [
            r'O-methyltransferase',
            r'\bOMT\b',
            r'\bCOMT\b',
            r'caffeic\s+acid.*methyltransferase',
            r'caffeoyl.CoA\s+O-methyltransferase',
            r'\bCCoAOMT\b',
            r'catechol\s+O-methyltransferase',
            r'flavonoid.*O-methyltransferase',
            r'isoflavone.*O-methyltransferase',
            r'myricetin.*O-methyltransferase',
            r'trans-resveratrol.*O-methyltransferase',
        ],
        'blast_exclude': [
            r'DNA\s+methyltransferase',
            r'histone.*methyltransferase',
            r'rRNA.*methyltransferase',
            r'tRNA.*methyltransferase',
            r'N-methyltransferase',
            r'protein.*methyltransferase',
        ],
        'hmmer_pfam': ['PF00891', 'PF01596', 'PF08100'],
        'hmmer_names': ['Methyltransf_2', 'Methyltransf_3', 'Dimerisation'],
    },
}

FAMILY_COLORS = {'CYP': '#e74c3c', 'OMT': '#3498db'}
COND_COLORS = {'R': '#d35400', 'L': '#27ae60',
               'root': '#d35400', 'leaf': '#27ae60'}


# ============================================================================
# GENE FAMILY DETECTION
# ============================================================================

def detect_gene_families(results_df):
    """
    Identify CYP and OMT genes from either gene_family column or
    blast_description patterns. Returns a dict mapping gene_id -> family name.
    """
    family_map = {}

    if 'gene_family' in results_df.columns:
        print("  Detecting gene families from 'gene_family' column...")
        for idx, row in results_df.iterrows():
            gid = row.get('gene_id', idx)
            fam = row['gene_family']
            if pd.notna(fam) and fam in GENE_FAMILIES:
                family_map[gid] = fam
        print(f"    Found {len(family_map)} annotated genes")
        return family_map

    if 'blast_description' not in results_df.columns:
        print("  No blast_description or gene_family column -- skipping family detection")
        return family_map

    print("  Detecting gene families from blast_description patterns...")
    desc_col = results_df['blast_description'].fillna('')

    for fam_name, fam_def in GENE_FAMILIES.items():
        include_mask = pd.Series(False, index=results_df.index)
        for pattern in fam_def['blast_patterns']:
            include_mask |= desc_col.str.contains(pattern, case=False, regex=True)

        exclude_mask = pd.Series(False, index=results_df.index)
        for pattern in fam_def['blast_exclude']:
            exclude_mask |= desc_col.str.contains(pattern, case=False, regex=True)

        hits = results_df[include_mask & ~exclude_mask]
        for idx in hits.index:
            gid = hits.at[idx, 'gene_id'] if 'gene_id' in hits.columns else idx
            if gid not in family_map:
                family_map[gid] = fam_name

        n = (include_mask & ~exclude_mask).sum()
        print(f"    {fam_name}: {n} genes matched")

    return family_map


def parse_hmmer_domains(hmmer_file):
    """Parse HMMER domtblout file, return dict: gene_id -> domain string."""
    if hmmer_file is None or not Path(hmmer_file).exists():
        return {}

    print(f"  Reading HMMER domains: {hmmer_file}")
    domain_map = {}
    with open(hmmer_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 23:
                continue
            gid = parts[3].split('|')[0]
            domain = f"{parts[0]}({parts[1]})"
            if gid not in domain_map:
                domain_map[gid] = set()
            domain_map[gid].add(domain)

    for gid in domain_map:
        domain_map[gid] = '; '.join(sorted(domain_map[gid]))

    print(f"    Loaded domains for {len(domain_map)} genes")
    return domain_map


# ============================================================================
# GTF BIOTYPE PARSING
# ============================================================================

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
    """Extract gene_id -> gene_biotype mapping from a GTF file.

    Scans gene-level rows for the 'gene_biotype' or 'gene_type' attribute.
    Returns an empty dict if the file doesn't exist or contains no biotypes.
    """
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
            attrs = fields[8]
            gid = bt = None
            for token in attrs.split(';'):
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

    print(f"    Found biotype for {len(biotype_map)} genes "
          f"({len(set(biotype_map.values()))} categories)")
    return biotype_map


def _biotype_row_colors(gene_ids, biotype_map):
    """Build a Series of colors from biotype_map for the given gene_ids."""
    if not biotype_map:
        return None
    colors = []
    for gid in gene_ids:
        bt = biotype_map.get(gid, 'other')
        colors.append(BIOTYPE_PALETTE.get(bt, '#999999'))
    return colors


# ============================================================================
# NORMALIZATION
# ============================================================================

def calculate_normalized_counts(count_matrix):
    """DESeq2-style median-of-ratios normalization."""
    print("\n  Calculating normalized counts...")

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
    print("    Size factors: " + ", ".join(f"{s}={sf:.3f}" for s, sf in size_factors.items()))

    return count_matrix.div(size_factors, axis=1)


# ============================================================================
# MA + VOLCANO PLOTS
# ============================================================================

def generate_ma_plot(results_df, output_dir, alpha=0.05):
    """Generate MA plot with blue significant points."""
    print("  Creating MA plot...")
    fig, ax = plt.subplots(figsize=(8, 6))

    valid = (results_df['log2FoldChange'].notna() & results_df['baseMean'].notna()
             & np.isfinite(results_df['log2FoldChange']) & np.isfinite(results_df['baseMean']))

    ax.scatter(results_df.loc[valid, 'baseMean'],
               results_df.loc[valid, 'log2FoldChange'],
               alpha=0.4, s=1.5, c='#AAAAAA', rasterized=True)

    sig = valid & results_df['padj'].notna() & (results_df['padj'] < alpha)
    if sig.sum() > 0:
        ax.scatter(results_df.loc[sig, 'baseMean'],
                   results_df.loc[sig, 'log2FoldChange'],
                   alpha=0.6, s=2, c='#3366CC', rasterized=True,
                   label=f'padj < {alpha} (n={sig.sum():,})')

    ax.axhline(y=0, color='black', linestyle='-', linewidth=0.8)
    ax.set_xlabel('Mean of normalized counts', fontsize=12)
    ax.set_ylabel('Log$_2$ fold change', fontsize=12)
    ax.set_title(f'MA plot with alpha = {alpha}', fontsize=14, fontweight='bold')
    ax.set_xscale('log')
    ax.legend(loc='best', fontsize=9, framealpha=0.9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_dir / "ma_plot.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / "ma_plot.png", dpi=150, bbox_inches='tight')
    plt.close()
    print(f"    Saved: {output_dir / 'ma_plot.pdf'}")


def generate_volcano_plot(results_df, output_dir, padj_cutoff=0.05,
                          lfc_cutoff=2.0, top_n_labels=10):
    """Generate Enhanced Volcano plot with 4-color categories."""
    print("  Creating enhanced volcano plot...")
    fig, ax = plt.subplots(figsize=(10, 8))

    valid = (results_df['log2FoldChange'].notna() & results_df['padj'].notna()
             & np.isfinite(results_df['log2FoldChange']) & np.isfinite(results_df['padj']))
    df = results_df.loc[valid].copy()
    df['neg_log10_padj'] = -np.log10(df['padj'] + 1e-300)

    pass_padj = df['padj'] < padj_cutoff
    pass_lfc = df['log2FoldChange'].abs() > lfc_cutoff

    cat_ns = ~pass_padj & ~pass_lfc
    cat_lfc = pass_lfc & ~pass_padj
    cat_padj = pass_padj & ~pass_lfc
    cat_both = pass_padj & pass_lfc

    for mask, color, label in [
        (cat_ns,   '#AAAAAA', 'NS'),
        (cat_lfc,  '#2ca02c', f'Log$_2$ FC'),
        (cat_padj, '#1f77b4', f'adj. p-value'),
        (cat_both, '#d62728', f'adj. p-value and Log$_2$ FC'),
    ]:
        subset = df.loc[mask]
        if len(subset) == 0:
            continue
        ax.scatter(subset['log2FoldChange'], subset['neg_log10_padj'],
                   alpha=0.5, s=2, c=color, label=label, rasterized=True)

    ax.axhline(y=-np.log10(padj_cutoff), color='black', linestyle='--',
               linewidth=0.6, alpha=0.5)
    ax.axvline(x=lfc_cutoff, color='black', linestyle='--',
               linewidth=0.6, alpha=0.5)
    ax.axvline(x=-lfc_cutoff, color='black', linestyle='--',
               linewidth=0.6, alpha=0.5)

    if top_n_labels > 0 and cat_both.sum() > 0:
        sig_df = df.loc[cat_both].nlargest(top_n_labels, 'neg_log10_padj')
        gene_id_col = 'gene_id' if 'gene_id' in sig_df.columns else None
        for idx, row in sig_df.iterrows():
            label_text = row[gene_id_col] if gene_id_col else str(idx)
            ax.annotate(label_text,
                        xy=(row['log2FoldChange'], row['neg_log10_padj']),
                        fontsize=5.5, alpha=0.8,
                        xytext=(5, 3), textcoords='offset points',
                        arrowprops=dict(arrowstyle='-', color='gray',
                                        lw=0.4, alpha=0.5))

    ax.set_xlabel('Log$_2$ fold change', fontsize=12)
    ax.set_ylabel('-Log$_{10}$ p', fontsize=12)
    ax.set_title('Volcano plot', fontsize=14, fontweight='bold')
    ax.text(0.5, 1.01, 'EnhancedVolcano', transform=ax.transAxes,
            ha='center', va='bottom', fontsize=9, fontstyle='italic',
            color='#666666')

    n_total = len(df)
    ax.text(0.5, -0.08, f'Total = {n_total:,} variables',
            transform=ax.transAxes, ha='center', va='top', fontsize=9)

    ax.legend(loc='upper right', fontsize=8, framealpha=0.9, markerscale=3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_dir / "volcano_plot.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / "volcano_plot.png", dpi=150, bbox_inches='tight')
    plt.close()
    print(f"    Saved: {output_dir / 'volcano_plot.pdf'}")


def generate_ma_volcano_plots(results_df, output_dir):
    """Generate MA and Enhanced Volcano plots."""
    print("\nGenerating MA and Volcano plots...")
    generate_ma_plot(results_df, output_dir)
    generate_volcano_plot(results_df, output_dir)


# ============================================================================
# PCA PLOT
# ============================================================================

def generate_pca_plot(count_matrix, metadata, output_dir, n_top=50,
                      condition_col='condition'):
    """PCA of top-N most variable genes, colored by condition."""
    print("\nGenerating PCA plot...")

    counts = count_matrix[~count_matrix.index.str.startswith('N_')]
    norm = calculate_normalized_counts(counts)
    log_norm = np.log2(norm + 1)

    variances = log_norm.var(axis=1)
    top_genes = variances.nlargest(n_top).index
    subset = log_norm.loc[top_genes]

    X = subset.values.T  # samples x genes
    X_centered = X - X.mean(axis=0)

    U, S, Vt = np.linalg.svd(X_centered, full_matrices=False)
    explained = (S ** 2) / (S ** 2).sum() * 100
    pc1 = U[:, 0] * S[0]
    pc2 = U[:, 1] * S[1]

    meta = metadata.copy()
    sample_col = meta.columns[0]
    if sample_col != condition_col:
        meta = meta.set_index(sample_col)

    fig, ax = plt.subplots(figsize=(7, 5))
    cond_palette = {'R': '#d35400', 'L': '#27ae60',
                    'root': '#d35400', 'leaf': '#27ae60'}
    cond_labels = {'R': 'Root', 'L': 'Leaf',
                   'root': 'Root', 'leaf': 'Leaf'}

    for i, sample in enumerate(subset.columns):
        cond = meta.loc[sample, condition_col] if sample in meta.index else 'unknown'
        color = cond_palette.get(cond, '#bdc3c7')
        label = cond_labels.get(cond, str(cond))
        ax.scatter(pc1[i], pc2[i], c=color, s=80, edgecolors='white',
                   linewidths=0.5, zorder=3)

    handles = []
    for cond_code, color in cond_palette.items():
        label = cond_labels.get(cond_code, cond_code)
        if label not in [h.get_label() for h in handles]:
            handles.append(plt.scatter([], [], c=color, s=60,
                                       edgecolors='white', label=label))

    ax.set_xlabel(f'PC1: {explained[0]:.0f}% variance', fontsize=12)
    ax.set_ylabel(f'PC2: {explained[1]:.0f}% variance', fontsize=12)
    ax.set_title(f'PC1 vs PC2: top {n_top} variable genes',
                 fontsize=13, fontweight='bold')
    ax.legend(handles=handles, title='condition', fontsize=9,
              title_fontsize=10, loc='best', framealpha=0.9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.axhline(0, color='gray', linewidth=0.3, alpha=0.5)
    ax.axvline(0, color='gray', linewidth=0.3, alpha=0.5)

    plt.tight_layout()
    plt.savefig(output_dir / "pca_plot.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / "pca_plot.png", dpi=150, bbox_inches='tight')
    plt.close()
    print(f"    Saved: {output_dir / 'pca_plot.pdf'}")


# ============================================================================
# SAMPLE CORRELATION HEATMAP
# ============================================================================

def generate_sample_correlation_heatmap(count_matrix, metadata, output_dir,
                                        condition_col='condition'):
    """Sample-to-sample distance heatmap with blue gradient."""
    print("\nGenerating sample correlation heatmap...")

    counts = count_matrix[~count_matrix.index.str.startswith('N_')]
    norm = calculate_normalized_counts(counts)
    log_norm = np.log2(norm + 1)

    from scipy.spatial.distance import pdist, squareform
    dist_vec = pdist(log_norm.T.values, metric='euclidean')
    dist_mat = squareform(dist_vec)
    dist_df = pd.DataFrame(dist_mat, index=log_norm.columns,
                           columns=log_norm.columns)

    meta = metadata.copy()
    sample_col = meta.columns[0]
    if sample_col != condition_col:
        meta = meta.set_index(sample_col)

    cond_palette = {'R': '#d35400', 'L': '#27ae60',
                    'root': '#d35400', 'leaf': '#27ae60'}
    col_colors = []
    for sample in dist_df.columns:
        cond = meta.loc[sample, condition_col] if sample in meta.index else ''
        col_colors.append(cond_palette.get(cond, '#bdc3c7'))
    col_colors_series = pd.Series(col_colors, index=dist_df.columns)

    n = len(dist_df)
    fig_size = max(6, 0.6 * n + 2)

    g = sns.clustermap(
        dist_df,
        cmap='Blues_r',
        figsize=(fig_size, fig_size),
        linewidths=0.5,
        linecolor='white',
        col_colors=col_colors_series,
        row_colors=col_colors_series,
        xticklabels=True,
        yticklabels=True,
        cbar_kws={'label': 'Euclidean distance', 'shrink': 0.5},
        dendrogram_ratio=0.12,
    )

    g.ax_heatmap.set_xlabel('')
    g.ax_heatmap.set_ylabel('')
    g.fig.suptitle('Sample-to-Sample Distance', fontsize=13,
                   fontweight='bold', y=1.02)

    legend_elements = [
        Patch(facecolor='#d35400', label='Root'),
        Patch(facecolor='#27ae60', label='Leaf'),
    ]
    g.ax_heatmap.legend(handles=legend_elements, loc='upper left',
                        bbox_to_anchor=(1.12, 1.0), frameon=True, fontsize=9)

    plt.savefig(output_dir / "sample_correlation_heatmap.pdf",
                dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / "sample_correlation_heatmap.png",
                dpi=150, bbox_inches='tight')
    plt.close()
    print(f"    Saved: {output_dir / 'sample_correlation_heatmap.pdf'}")


# ============================================================================
# GENE FAMILY HEATMAP
# ============================================================================

def generate_family_heatmap(results_df, gene_ids, family_name, full_name,
                            count_matrix, metadata, domain_map, output_dir,
                            scale='center', cluster_rows=True,
                            condition_col='condition', biotype_map=None):
    """
    Generate a heatmap for a single gene family (CYP or OMT).
    gene_ids: list of gene IDs belonging to this family.
    """
    print(f"\n{'='*60}")
    print(f"Generating {family_name} Heatmap ({full_name})")
    print(f"{'='*60}")

    # Drop STAR summary rows from count matrix
    counts = count_matrix[~count_matrix.index.str.startswith('N_')]

    present = [g for g in gene_ids if g in counts.index]
    missing = [g for g in gene_ids if g not in counts.index]
    if missing:
        print(f"  WARNING: {len(missing)} {family_name} genes not in count matrix")
    if not present:
        print(f"  No {family_name} genes in count matrix -- skipping heatmap")
        return

    print(f"  Normalizing counts for {len(present)} {family_name} genes...")
    norm_counts = calculate_normalized_counts(counts)
    heatmap_data = np.log2(norm_counts.loc[present] + 1)

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

    # Order columns by condition
    meta = metadata.copy()
    if condition_col in meta.columns:
        sample_col = meta.columns[0] if meta.index.name is None else None
        if sample_col and sample_col != condition_col:
            meta = meta.set_index(sample_col)
        ordered = [s for s in meta.sort_values(condition_col).index if s in heatmap_data.columns]
        heatmap_data = heatmap_data[ordered]

    # Build row labels: gene_id | short description [domain]
    row_labels = []
    id_col = 'gene_id'
    has_gene_id = id_col in results_df.columns

    for gid in heatmap_data.index:
        label = gid
        if has_gene_id:
            match = results_df[results_df['gene_id'] == gid]
        else:
            match = results_df.loc[[gid]] if gid in results_df.index else pd.DataFrame()

        if not match.empty:
            row = match.iloc[0]
            desc = row.get('blast_description', '')
            if pd.notna(desc) and desc:
                label = f"{gid} | {str(desc)[:40]}"

        if gid in domain_map:
            label += f"  [{domain_map[gid]}]"

        row_labels.append(label)

    heatmap_data.index = row_labels

    # Column colors for condition
    col_colors = None
    if condition_col in meta.columns:
        col_colors_list = []
        for sample in heatmap_data.columns:
            cond = meta.loc[sample, condition_col] if sample in meta.index else None
            col_colors_list.append(COND_COLORS.get(cond, '#bdc3c7'))
        col_colors = pd.Series(col_colors_list, index=heatmap_data.columns)

    # Figure sizing
    n_genes = len(heatmap_data)
    n_samples = len(heatmap_data.columns)
    fig_height = max(6, 0.35 * n_genes + 2)
    fig_width = max(8, 0.8 * n_samples + 4)

    # Biotype row colors
    row_colors = None
    biotype_legend = []
    bt_colors = _biotype_row_colors(present, biotype_map or {})
    if bt_colors:
        row_colors = pd.Series(bt_colors, index=heatmap_data.index)
        seen = set()
        for gid, c in zip(present, bt_colors):
            bt = (biotype_map or {}).get(gid, 'other')
            if bt not in seen:
                seen.add(bt)
                biotype_legend.append(Patch(facecolor=c, label=bt))

    print(f"  Plotting {n_genes} genes x {n_samples} samples...")

    prefix = family_name.lower()
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
        g.fig.suptitle(f'{full_name} ({family_name}) Expression ({n_genes} genes)',
                       fontsize=14, fontweight='bold', y=1.02)

        legend_elements = [
            Patch(facecolor='#d35400', label='Root (R)'),
            Patch(facecolor='#27ae60', label='Leaf (L)'),
        ]
        if biotype_legend:
            legend_elements.append(Patch(facecolor='none', edgecolor='none', label=''))
            legend_elements.extend(biotype_legend)
        g.ax_heatmap.legend(handles=legend_elements, loc='upper left',
                            bbox_to_anchor=(1.15, 1.0), frameon=True, fontsize=9)

        pdf_path = output_dir / f"{prefix}_heatmap.pdf"
        g.savefig(pdf_path, dpi=300, bbox_inches='tight')
        g.savefig(output_dir / f"{prefix}_heatmap.png", dpi=150, bbox_inches='tight')
        plt.close()
        print(f"    Saved: {pdf_path}")

        matrix_path = output_dir / f"{prefix}_heatmap_matrix.tsv"
        heatmap_data.to_csv(matrix_path, sep='\t')
        print(f"    Saved: {matrix_path}")

    except Exception as e:
        print(f"  ERROR generating {family_name} heatmap: {e}")
        import traceback
        traceback.print_exc()

    # Save gene list for this family
    gene_list_path = output_dir / f"{prefix}_gene_list.tsv"
    if has_gene_id:
        family_df = results_df[results_df['gene_id'].isin(gene_ids)]
    else:
        family_df = results_df.loc[results_df.index.intersection(gene_ids)]
    family_df.to_csv(gene_list_path, sep='\t', index=not has_gene_id)
    print(f"    Saved gene list: {gene_list_path} ({len(family_df)} genes)")


# ============================================================================
# COMBINED TWO-SPECIES HEATMAP
# ============================================================================

def _order_samples(meta, available_cols, condition_col):
    """Return sample names ordered root-first then leaf."""
    sample_col = meta.columns[0]
    root_conds = {'R', 'root'}
    leaf_conds = {'L', 'leaf'}
    root_samples, leaf_samples, other_samples = [], [], []
    for _, row in meta.iterrows():
        s = row[sample_col]
        if s not in available_cols:
            continue
        cond = row.get(condition_col, '')
        if cond in root_conds:
            root_samples.append(s)
        elif cond in leaf_conds:
            leaf_samples.append(s)
        else:
            other_samples.append(s)
    return root_samples + leaf_samples + other_samples


def _normalize_and_log(count_path, gene_ids):
    """Load count matrix, normalize, log2-transform, subset to gene_ids."""
    counts = pd.read_csv(count_path, sep='\t', index_col=0)
    counts = counts[~counts.index.str.startswith('N_')]
    present = [g for g in gene_ids if g in counts.index]
    norm = calculate_normalized_counts(counts)
    return np.log2(norm.loc[present] + 1), present


def generate_combined_family_heatmap(
        results_df, gene_ids, family_name, full_name,
        count_matrix_path1, metadata_path1,
        count_matrix_path2, metadata_path2,
        species1, species2, domain_map, output_dir,
        scale='center', cluster_rows=True, condition_col='condition',
        biotype_map=None):
    """Generate a single heatmap with two species side by side.

    Columns: SP1-Root | SP1-Leaf | SP2-Root | SP2-Leaf.
    Two stacked column color bars for species and tissue type.
    """

    print(f"\n{'='*60}")
    print(f"Generating Combined {family_name} Heatmap: {species1} + {species2}")
    print(f"{'='*60}")

    if not gene_ids:
        print(f"  No {family_name} genes -- skipping")
        return

    print(f"  Loading and normalizing {species1} counts...")
    log1, present1 = _normalize_and_log(count_matrix_path1, gene_ids)
    print(f"    {species1}: {len(present1)} genes found")

    print(f"  Loading and normalizing {species2} counts...")
    log2, present2 = _normalize_and_log(count_matrix_path2, gene_ids)
    print(f"    {species2}: {len(present2)} genes found")

    shared = [g for g in gene_ids if g in log1.index and g in log2.index]
    if not shared:
        print(f"  No {family_name} genes found in both count matrices -- skipping")
        return
    print(f"  Shared genes: {len(shared)}")

    log1 = log1.loc[shared]
    log2 = log2.loc[shared]

    meta1 = pd.read_csv(metadata_path1, sep='\t')
    meta2 = pd.read_csv(metadata_path2, sep='\t')

    sp1_ordered = _order_samples(meta1, set(log1.columns), condition_col)
    sp2_ordered = _order_samples(meta2, set(log2.columns), condition_col)

    heatmap_data = pd.concat([log1[sp1_ordered], log2[sp2_ordered]], axis=1)

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

    # Row labels with description + HMMER domain
    has_gene_id = 'gene_id' in results_df.columns
    row_labels = []
    for gid in heatmap_data.index:
        label = gid
        if has_gene_id:
            match = results_df[results_df['gene_id'] == gid]
        else:
            match = results_df.loc[[gid]] if gid in results_df.index else pd.DataFrame()
        if not match.empty:
            desc = match.iloc[0].get('blast_description', '')
            if pd.notna(desc) and desc:
                label = f"{gid} | {str(desc)[:40]}"
        if gid in domain_map:
            label += f"  [{domain_map[gid]}]"
        row_labels.append(label)
    heatmap_data.index = row_labels

    # Two-row column color bars
    sp1_set = set(sp1_ordered)
    species_palette = {species1: '#8e44ad', species2: '#2980b9'}
    tissue_palette = {'Root': '#d35400', 'Leaf': '#27ae60'}
    root_conds = {'R', 'root'}
    leaf_conds = {'L', 'leaf'}

    sp1_sample_col = meta1.columns[0]
    sp2_sample_col = meta2.columns[0]

    species_colors = []
    tissue_colors = []
    for sample in heatmap_data.columns:
        if sample in sp1_set:
            species_colors.append(species_palette[species1])
            match = meta1[meta1[sp1_sample_col] == sample]
            cond = match[condition_col].values[0] if not match.empty else ''
        else:
            species_colors.append(species_palette[species2])
            match = meta2[meta2[sp2_sample_col] == sample]
            cond = match[condition_col].values[0] if not match.empty else ''

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

    n_genes = len(heatmap_data)
    n_samples = len(heatmap_data.columns)
    fig_height = max(6, 0.35 * n_genes + 2)
    fig_width = max(10, 0.8 * n_samples + 4)

    # Biotype row colors
    row_colors = None
    biotype_legend = []
    bt_colors = _biotype_row_colors(shared, biotype_map or {})
    if bt_colors:
        row_colors = pd.Series(bt_colors, index=heatmap_data.index)
        seen = set()
        for gid, c in zip(shared, bt_colors):
            bt = (biotype_map or {}).get(gid, 'other')
            if bt not in seen:
                seen.add(bt)
                biotype_legend.append(Patch(facecolor=c, label=bt))

    print(f"  Plotting {n_genes} genes x {n_samples} samples...")
    prefix = family_name.lower()

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
        g.fig.suptitle(
            f'{full_name} ({family_name}) Expression: {species1} vs {species2} ({n_genes} genes)',
            fontsize=14, fontweight='bold', y=1.02,
        )

        legend_elements = [
            Patch(facecolor=species_palette[species1], label=species1),
            Patch(facecolor=species_palette[species2], label=species2),
            Patch(facecolor='none', edgecolor='none', label=''),
            Patch(facecolor=tissue_palette['Root'], label='Root'),
            Patch(facecolor=tissue_palette['Leaf'], label='Leaf'),
        ]
        if biotype_legend:
            legend_elements.append(Patch(facecolor='none', edgecolor='none', label=''))
            legend_elements.extend(biotype_legend)
        g.ax_heatmap.legend(handles=legend_elements, loc='upper left',
                            bbox_to_anchor=(1.15, 1.0), frameon=True, fontsize=9)

        pdf_path = output_dir / f"{prefix}_heatmap_combined.pdf"
        g.savefig(pdf_path, dpi=300, bbox_inches='tight')
        g.savefig(output_dir / f"{prefix}_heatmap_combined.png", dpi=150, bbox_inches='tight')
        plt.close()
        print(f"    Saved: {pdf_path}")

        matrix_path = output_dir / f"{prefix}_heatmap_combined_matrix.tsv"
        heatmap_data.to_csv(matrix_path, sep='\t')
        print(f"    Saved: {matrix_path}")

    except Exception as e:
        print(f"  ERROR generating combined {family_name} heatmap: {e}")
        import traceback
        traceback.print_exc()


# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="PyDESeq2 Step 3: Generate plots (MA, volcano, CYP/OMT heatmaps)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # From combined annotated file (auto-detects CYP/OMT):
  python pydeseq2_generate_plots.py combined_annotated.tsv \\
      --count-matrix count_matrix.tsv --metadata metadata.tsv -o plots/

  # With HMMER domain labels:
  python pydeseq2_generate_plots.py combined_annotated.tsv \\
      --count-matrix count_matrix.tsv --metadata metadata.tsv \\
      --hmmer pfam_domains.txt -o plots/

  # Raw DESeq2 results (MA + volcano only):
  python pydeseq2_generate_plots.py pydeseq2_results.tsv -o plots/
        """
    )

    parser.add_argument("results_file",
                        help="DESeq2 results TSV (unfiltered, filtered, combined annotated, or gene-family verified)")

    parser.add_argument("-o", "--output", default="plots",
                        help="Output directory (default: plots)")

    parser.add_argument("--count-matrix", default=None,
                        help="Count matrix TSV (required for heatmaps)")
    parser.add_argument("--metadata", default=None,
                        help="Metadata TSV (required for heatmaps)")
    parser.add_argument("--hmmer", default=None,
                        help="HMMER domtblout file (optional, adds domain labels to heatmap rows)")

    parser.add_argument("--contrast-factor", default="condition",
                        help="Metadata column for contrast (default: condition)")
    parser.add_argument("--contrast-A", default="R",
                        help="Numerator condition code (default: R)")
    parser.add_argument("--contrast-B", default="L",
                        help="Denominator condition code (default: L)")

    parser.add_argument("--scale", choices=['center', 'zscore'], default='center',
                        help="Heatmap scaling method (default: center)")
    parser.add_argument("--no-row-cluster", action="store_true",
                        help="Disable row clustering in heatmaps")
    parser.add_argument("--gtf", default=None,
                        help="GTF annotation file (optional, adds gene biotype sidebar to heatmaps)")

    # Combined two-species heatmap mode
    parser.add_argument("--count-matrix2", default=None,
                        help="Second species count matrix (enables combined heatmap)")
    parser.add_argument("--metadata2", default=None,
                        help="Second species metadata TSV")
    parser.add_argument("--species1", default="SP1",
                        help="Label for species 1 (default: SP1)")
    parser.add_argument("--species2", default="SP2",
                        help="Label for species 2 (default: SP2)")

    args = parser.parse_args()

    if not os.path.exists(args.results_file):
        print(f"ERROR: Results file not found: {args.results_file}")
        return 1

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("PyDESeq2 Plot Generation (Step 3)")
    print("=" * 60)
    print(f"Input:  {args.results_file}")
    print(f"Output: {output_dir}")
    print()

    try:
        results_df = pd.read_csv(args.results_file, sep='\t', comment='#', index_col=0)
        print(f"Loaded {len(results_df)} genes")

        # Detect input type
        has_blast = 'blast_description' in results_df.columns
        has_family = 'gene_family' in results_df.columns
        print(f"  Input has blast_description: {has_blast}")
        print(f"  Input has gene_family:       {has_family}")

        # --- MA + Volcano (always) ---
        if 'log2FoldChange' in results_df.columns and 'baseMean' in results_df.columns:
            generate_ma_volcano_plots(results_df, output_dir)
        else:
            print("\n  Skipping MA/Volcano (no log2FoldChange or baseMean columns)")

        # --- Plots requiring count matrix + metadata ---
        if not args.count_matrix or not args.metadata:
            print("\n  Count matrix and/or metadata not provided -- skipping heatmaps, PCA, correlation.")
            print("  For heatmaps, provide --count-matrix and --metadata")
        elif not os.path.exists(args.count_matrix):
            print(f"\n  WARNING: Count matrix not found: {args.count_matrix}")
        elif not os.path.exists(args.metadata):
            print(f"\n  WARNING: Metadata not found: {args.metadata}")
        else:
            count_matrix = pd.read_csv(args.count_matrix, sep='\t', index_col=0)
            metadata = pd.read_csv(args.metadata, sep='\t')

            # --- PCA plot ---
            try:
                generate_pca_plot(count_matrix, metadata, output_dir,
                                  condition_col=args.contrast_factor)
            except Exception as e:
                print(f"  WARNING: PCA plot failed: {e}")

            # --- Sample correlation heatmap ---
            try:
                generate_sample_correlation_heatmap(
                    count_matrix, metadata, output_dir,
                    condition_col=args.contrast_factor)
            except Exception as e:
                print(f"  WARNING: Sample correlation heatmap failed: {e}")

            # --- Gene family heatmaps ---
            domain_map = parse_hmmer_domains(args.hmmer)
            biotype_map = parse_biotype_from_gtf(args.gtf)
            family_map = detect_gene_families(results_df)

            if not family_map:
                print("\n  No CYP or OMT genes detected -- skipping family heatmaps.")
                print("  (Input may be raw DESeq2 results without BLAST annotation)")
            else:
                for fam_name, fam_def in GENE_FAMILIES.items():
                    gene_ids = [gid for gid, fam in family_map.items() if fam == fam_name]
                    if not gene_ids:
                        print(f"\n  No {fam_name} genes found -- skipping {fam_name} heatmap")
                        continue
                    generate_family_heatmap(
                        results_df, gene_ids, fam_name, fam_def['full_name'],
                        count_matrix, metadata, domain_map, output_dir,
                        scale=args.scale,
                        cluster_rows=not args.no_row_cluster,
                        condition_col=args.contrast_factor,
                        biotype_map=biotype_map,
                    )

                # Combined two-species heatmaps
                if (args.count_matrix2 and args.metadata2
                        and os.path.exists(args.count_matrix2)
                        and os.path.exists(args.metadata2)):
                    print(f"\n  Generating combined heatmaps ({args.species1} + {args.species2})...")
                    for fam_name, fam_def in GENE_FAMILIES.items():
                        fam_gene_ids = [gid for gid, fam in family_map.items() if fam == fam_name]
                        if not fam_gene_ids:
                            continue
                        generate_combined_family_heatmap(
                            results_df, fam_gene_ids, fam_name, fam_def['full_name'],
                            args.count_matrix, args.metadata,
                            args.count_matrix2, args.metadata2,
                            args.species1, args.species2,
                            domain_map, output_dir,
                            scale=args.scale,
                            cluster_rows=not args.no_row_cluster,
                            condition_col=args.contrast_factor,
                            biotype_map=biotype_map,
                        )

        print("\n" + "=" * 60)
        print("Plot generation complete!")
        print("=" * 60)
        print(f"\nOutput: {output_dir}")
        return 0

    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
