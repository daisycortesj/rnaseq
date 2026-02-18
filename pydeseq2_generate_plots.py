#!/usr/bin/env python3
"""
PyDESeq2 Step 3: Generate Plots
================================

Generates publication-quality plots from DESeq2 results:
  - MA plot (baseMean vs log2FC)
  - Volcano plot (log2FC vs -log10 padj)
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

def generate_ma_volcano_plots(results_df, output_dir):
    """Generate MA and Volcano plots."""
    print("\nGenerating MA and Volcano plots...")

    # -- MA plot --
    print("  Creating MA plot...")
    fig, ax = plt.subplots(figsize=(8, 6))

    valid = (results_df['log2FoldChange'].notna() & results_df['baseMean'].notna()
             & np.isfinite(results_df['log2FoldChange']) & np.isfinite(results_df['baseMean']))

    ax.scatter(results_df.loc[valid, 'baseMean'],
               results_df.loc[valid, 'log2FoldChange'],
               alpha=0.5, s=2, c='gray', label='Not significant')

    sig = valid & results_df['padj'].notna() & (results_df['padj'] < 0.05)
    if sig.sum() > 0:
        ax.scatter(results_df.loc[sig, 'baseMean'],
                   results_df.loc[sig, 'log2FoldChange'],
                   alpha=0.7, s=3, c='red', label=f'padj < 0.05 (n={sig.sum()})')

    ax.axhline(y=0, color='black', linestyle='--', linewidth=0.5)
    ax.set_xlabel('Mean Normalized Counts (baseMean)', fontsize=12)
    ax.set_ylabel('Log2 Fold Change', fontsize=12)
    ax.set_title('MA Plot', fontsize=14, fontweight='bold')
    ax.set_xscale('log')
    ax.legend()
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_dir / "ma_plot.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / "ma_plot.png", dpi=150, bbox_inches='tight')
    plt.close()
    print(f"    Saved: {output_dir / 'ma_plot.pdf'}")

    # -- Volcano plot --
    print("  Creating volcano plot...")
    fig, ax = plt.subplots(figsize=(10, 8))

    valid = (results_df['log2FoldChange'].notna() & results_df['padj'].notna()
             & np.isfinite(results_df['log2FoldChange']) & np.isfinite(results_df['padj']))

    neg_log10 = -np.log10(results_df.loc[valid, 'padj'] + 1e-300)

    ax.scatter(results_df.loc[valid, 'log2FoldChange'], neg_log10,
               alpha=0.5, s=2, c='gray', label='Not significant')

    sig = valid & (results_df['padj'] < 0.05)
    if sig.sum() > 0:
        sig_neg = -np.log10(results_df.loc[sig, 'padj'] + 1e-300)
        ax.scatter(results_df.loc[sig, 'log2FoldChange'], sig_neg,
                   alpha=0.7, s=3, c='red', label=f'padj < 0.05 (n={sig.sum()})')

    ax.axhline(y=-np.log10(0.05), color='black', linestyle='--', linewidth=0.5, label='padj = 0.05')
    ax.axvline(x=2, color='blue', linestyle='--', linewidth=0.5, alpha=0.5)
    ax.axvline(x=-2, color='blue', linestyle='--', linewidth=0.5, alpha=0.5, label='|log2FC| = 2')

    ax.set_xlabel('Log2 Fold Change', fontsize=12)
    ax.set_ylabel('-Log10 Adjusted P-value', fontsize=12)
    ax.set_title('Volcano Plot', fontsize=14, fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_dir / "volcano_plot.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / "volcano_plot.png", dpi=150, bbox_inches='tight')
    plt.close()
    print(f"    Saved: {output_dir / 'volcano_plot.pdf'}")


# ============================================================================
# GENE FAMILY HEATMAP
# ============================================================================

def generate_family_heatmap(results_df, gene_ids, family_name, full_name,
                            count_matrix, metadata, domain_map, output_dir,
                            scale='center', cluster_rows=True, condition_col='condition'):
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

        # --- Gene family heatmaps ---
        if not args.count_matrix or not args.metadata:
            print("\n  Count matrix and/or metadata not provided -- skipping heatmaps.")
            print("  For heatmaps, provide --count-matrix and --metadata")
        elif not os.path.exists(args.count_matrix):
            print(f"\n  WARNING: Count matrix not found: {args.count_matrix}")
        elif not os.path.exists(args.metadata):
            print(f"\n  WARNING: Metadata not found: {args.metadata}")
        else:
            count_matrix = pd.read_csv(args.count_matrix, sep='\t', index_col=0)
            metadata = pd.read_csv(args.metadata, sep='\t')

            domain_map = parse_hmmer_domains(args.hmmer)
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
