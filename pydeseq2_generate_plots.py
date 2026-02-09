#!/usr/bin/env python3
"""
PyDESeq2 Step 3: Generate Plots from Filtered Results
======================================================

This script generates publication-quality plots from filtered PyDESeq2 results:
  - CYP gene heatmap (Figure 6A style)
  - General DEG heatmap (all significant genes)
  - MA plot
  - Volcano plot

Input: 
  - Filtered results TSV (from Step 2)
  - Original count matrix (for expression data)
  - Metadata (for sample grouping)
  - CYP family mapping (optional, for CYP heatmap)

Usage:
  python pydeseq2_generate_plots.py filtered_results.tsv \\
      --count-matrix count_matrix.tsv \\
      --metadata metadata.tsv \\
      --cyp-family-map cyp_families.tsv \\
      -o plots/
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

try:
    import matplotlib.pyplot as plt
    import seaborn as sns
    from scipy.cluster.hierarchy import linkage, leaves_list
    from matplotlib.gridspec import GridSpec
except ImportError as e:
    print(f"ERROR: Missing required package: {e}")
    print("Please install required packages:")
    print("  pip install matplotlib seaborn scipy")
    sys.exit(1)


def load_cyp_family_map(filepath):
    """Load CYP gene family mapping."""
    print(f"Loading CYP family map: {filepath}")
    df = pd.read_csv(filepath, sep='\t')
    
    if 'gene_id' not in df.columns:
        raise ValueError("CYP family map must have 'gene_id' column")
    if 'cyp_family' not in df.columns:
        raise ValueError("CYP family map must have 'cyp_family' column")
    
    df = df.set_index('gene_id')
    print(f"  Loaded {len(df)} CYP gene mappings")
    return df


def calculate_normalized_counts(count_matrix):
    """
    Calculate normalized counts using median-of-ratios (DESeq2-style).
    This approximates what PyDESeq2 does internally.
    """
    print("\nCalculating normalized counts...")
    
    # Calculate geometric mean for each gene
    gene_geometric_means = count_matrix.apply(
        lambda row: np.exp(np.log(row[row > 0]).mean()) if (row > 0).sum() > 0 else 0,
        axis=1
    )
    
    # Calculate size factors for each sample
    size_factors = []
    for col in count_matrix.columns:
        ratios = count_matrix[col] / gene_geometric_means
        ratios = ratios[np.isfinite(ratios) & (ratios > 0)]
        size_factor = ratios.median()
        size_factors.append(size_factor)
    
    size_factors = pd.Series(size_factors, index=count_matrix.columns)
    
    print("  Size factors:")
    for sample, sf in size_factors.items():
        print(f"    {sample}: {sf:.3f}")
    
    # Normalize counts
    normalized_counts = count_matrix.div(size_factors, axis=1)
    
    return normalized_counts


def generate_ma_volcano_plots(results_df, output_dir):
    """Generate MA and Volcano plots."""
    print("\nGenerating MA and Volcano plots...")
    
    # MA plot
    print("  Creating MA plot...")
    fig, ax = plt.subplots(figsize=(8, 6))
    
    valid = (results_df['log2FoldChange'].notna()) & (results_df['baseMean'].notna())
    valid = valid & (np.isfinite(results_df['log2FoldChange'])) & (np.isfinite(results_df['baseMean']))
    
    # All genes in gray
    ax.scatter(
        results_df.loc[valid, 'baseMean'],
        results_df.loc[valid, 'log2FoldChange'],
        alpha=0.5, s=2, c='gray', label='Not significant'
    )
    
    # Significant genes in red
    sig_valid = valid & (results_df['padj'].notna()) & (results_df['padj'] < 0.05)
    if sig_valid.sum() > 0:
        ax.scatter(
            results_df.loc[sig_valid, 'baseMean'],
            results_df.loc[sig_valid, 'log2FoldChange'],
            alpha=0.7, s=3, c='red', label=f'padj < 0.05 (n={sig_valid.sum()})'
        )
    
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
    
    # Volcano plot
    print("  Creating volcano plot...")
    fig, ax = plt.subplots(figsize=(10, 8))
    
    valid = (results_df['log2FoldChange'].notna()) & (results_df['padj'].notna())
    valid = valid & (np.isfinite(results_df['log2FoldChange'])) & (np.isfinite(results_df['padj']))
    
    neg_log10_padj = -np.log10(results_df.loc[valid, 'padj'] + 1e-300)
    
    # All genes in gray
    ax.scatter(
        results_df.loc[valid, 'log2FoldChange'],
        neg_log10_padj,
        alpha=0.5, s=2, c='gray', label='Not significant'
    )
    
    # Significant genes in red
    sig_valid = valid & (results_df['padj'] < 0.05)
    if sig_valid.sum() > 0:
        sig_neg_log10 = -np.log10(results_df.loc[sig_valid, 'padj'] + 1e-300)
        ax.scatter(
            results_df.loc[sig_valid, 'log2FoldChange'],
            sig_neg_log10,
            alpha=0.7, s=3, c='red', label=f'padj < 0.05 (n={sig_valid.sum()})'
        )
    
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


def generate_cyp_heatmap(filtered_results, count_matrix, metadata, cyp_family_map, output_dir,
                         contrast_A='root', contrast_B='leaf', contrast_factor='condition',
                         scale_method='center', row_cluster=True):
    """Generate CYP gene expression heatmap (Figure 6A style)."""
    print("\n" + "=" * 60)
    print("Generating CYP Gene Heatmap")
    print("=" * 60)
    
    # Step 1: Filter to CYP genes
    print("\nStep 1: Filtering to CYP genes...")
    
    deg_gene_ids = set(filtered_results.index)
    cyp_gene_ids = set(cyp_family_map.index)
    cyp_deg_ids = deg_gene_ids & cyp_gene_ids
    
    if len(cyp_deg_ids) == 0:
        print("  ERROR: No CYP genes found in filtered results!")
        print(f"    Filtered DEGs: {len(deg_gene_ids)}")
        print(f"    CYP genes in map: {len(cyp_gene_ids)}")
        return
    
    print(f"  CYP DEGs found: {len(cyp_deg_ids)}")
    
    cyp_deg_df = filtered_results.loc[list(cyp_deg_ids)].copy()
    cyp_deg_df['cyp_family'] = cyp_family_map.loc[cyp_deg_df.index, 'cyp_family']
    
    # Save CYP DEG list
    cyp_deg_file = output_dir / "cyp_deg_list.tsv"
    cyp_deg_df.to_csv(cyp_deg_file, sep='\t')
    print(f"  Saved CYP DEG list: {cyp_deg_file}")
    
    # Step 2: Get normalized counts
    print("\nStep 2: Getting normalized expression data...")
    normalized_counts = calculate_normalized_counts(count_matrix)
    
    cyp_genes_in_matrix = [g for g in cyp_deg_df.index if g in normalized_counts.index]
    if len(cyp_genes_in_matrix) == 0:
        print("  ERROR: No CYP genes found in count matrix!")
        return
    
    cyp_norm_counts = normalized_counts.loc[cyp_genes_in_matrix]
    print(f"  Expression matrix: {cyp_norm_counts.shape[0]} genes × {cyp_norm_counts.shape[1]} samples")
    
    # Step 3: Transform: log2(normalized + 1)
    print(f"\nStep 3: Transforming data (log2 + {scale_method})...")
    heatmap_data = np.log2(cyp_norm_counts + 1)
    
    # Step 4: Scale per gene
    if scale_method == 'center':
        row_means = heatmap_data.mean(axis=1)
        heatmap_data = heatmap_data.subtract(row_means, axis=0)
        scale_label = "log2(norm+1), centered"
    elif scale_method == 'zscore':
        row_means = heatmap_data.mean(axis=1)
        row_stds = heatmap_data.std(axis=1)
        row_stds = row_stds.replace(0, 1)
        heatmap_data = heatmap_data.subtract(row_means, axis=0).div(row_stds, axis=0)
        scale_label = "z-score"
    
    # Step 5: Order samples by condition
    print("\nStep 4: Ordering samples...")
    if contrast_factor in metadata.columns:
        samples_B = sorted([s for s in heatmap_data.columns 
                          if s in metadata.index and metadata.loc[s, contrast_factor] == contrast_B])
        samples_A = sorted([s for s in heatmap_data.columns 
                          if s in metadata.index and metadata.loc[s, contrast_factor] == contrast_A])
        column_order = samples_B + samples_A
    else:
        column_order = list(heatmap_data.columns)
    
    heatmap_data = heatmap_data[column_order]
    
    # Step 6: Order genes by CYP family
    print("\nStep 5: Ordering genes by CYP family...")
    gene_families = cyp_deg_df.loc[heatmap_data.index, 'cyp_family']
    unique_families = sorted(gene_families.unique())
    print(f"  CYP families: {', '.join(unique_families)}")
    
    # Group genes by family and optionally cluster within families
    row_order = []
    for family in unique_families:
        family_genes = gene_families[gene_families == family].index.tolist()
        
        if len(family_genes) > 1 and row_cluster:
            family_data = heatmap_data.loc[family_genes]
            try:
                linkage_matrix = linkage(family_data, method='ward', metric='euclidean')
                ordered_idx = leaves_list(linkage_matrix)
                family_genes = [family_genes[i] for i in ordered_idx]
            except:
                pass
        
        row_order.extend(family_genes)
    
    heatmap_data = heatmap_data.loc[row_order]
    gene_families = gene_families.loc[row_order]
    
    # Step 7: Save heatmap matrix
    matrix_file = output_dir / "cyp_heatmap_matrix.tsv"
    heatmap_data.to_csv(matrix_file, sep='\t')
    print(f"\n  Saved heatmap matrix: {matrix_file}")
    
    # Step 8: Create plot
    print("\nStep 6: Creating heatmap plot...")
    
    n_genes = len(heatmap_data)
    n_samples = len(heatmap_data.columns)
    
    fig_width = max(10, 2.5 + n_samples * 0.5)
    fig_height = max(8, 1 + n_genes * 0.35)
    
    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = GridSpec(1, 4, figure=fig, width_ratios=[0.08, 0.18, 0.66, 0.08], wspace=0.02)
    
    ax_family = fig.add_subplot(gs[0, 0])
    ax_genes = fig.add_subplot(gs[0, 1])
    ax_heatmap = fig.add_subplot(gs[0, 2])
    ax_cbar = fig.add_subplot(gs[0, 3])
    
    # Draw heatmap
    vmin = -3 if scale_method == 'zscore' else None
    vmax = 3 if scale_method == 'zscore' else None
    
    im = ax_heatmap.imshow(heatmap_data.values, aspect='auto', cmap='RdBu_r',
                           vmin=vmin, vmax=vmax, interpolation='nearest')
    
    # X-axis: sample labels
    ax_heatmap.set_xticks(range(n_samples))
    ax_heatmap.set_xticklabels(column_order, rotation=45, ha='right', fontsize=10, fontweight='bold')
    ax_heatmap.xaxis.set_ticks_position('bottom')
    ax_heatmap.set_yticks([])
    ax_heatmap.set_xlim(-0.5, n_samples - 0.5)
    ax_heatmap.set_ylim(n_genes - 0.5, -0.5)
    
    # Gene IDs on left
    ax_genes.set_xlim(0, 1)
    ax_genes.set_ylim(n_genes - 0.5, -0.5)
    ax_genes.axis('off')
    
    for i, gene_id in enumerate(heatmap_data.index):
        ax_genes.text(0.95, i, gene_id, fontsize=7, va='center', ha='right',
                     fontfamily='monospace')
    
    # Family labels with brackets
    ax_family.set_xlim(0, 1)
    ax_family.set_ylim(n_genes - 0.5, -0.5)
    ax_family.axis('off')
    
    # Find family blocks
    current_family = None
    block_start = 0
    family_blocks = []
    
    for i, (gene, family) in enumerate(gene_families.items()):
        if family != current_family:
            if current_family is not None:
                family_blocks.append((current_family, block_start, i - 1))
            current_family = family
            block_start = i
    if current_family is not None:
        family_blocks.append((current_family, block_start, len(gene_families) - 1))
    
    # Draw family brackets
    for family, start, end in family_blocks:
        mid = (start + end) / 2
        height = end - start + 1
        
        bracket_x = 0.85
        ax_family.plot([bracket_x, bracket_x], [start - 0.3, end + 0.3],
                      color='black', linewidth=1.5, clip_on=False)
        ax_family.plot([bracket_x, bracket_x + 0.12], [start - 0.3, start - 0.3],
                      color='black', linewidth=1.5, clip_on=False)
        ax_family.plot([bracket_x, bracket_x + 0.12], [end + 0.3, end + 0.3],
                      color='black', linewidth=1.5, clip_on=False)
        
        ax_family.text(0.4, mid, family, fontsize=9, fontweight='bold',
                      ha='center', va='center', rotation=90 if height > 3 else 0)
    
    # Colorbar
    cbar = fig.colorbar(im, cax=ax_cbar)
    cbar.set_label(scale_label, fontsize=10)
    
    # Title
    fig.suptitle(f'CYP Gene Expression Heatmap ({n_genes} genes)',
                fontsize=14, fontweight='bold', y=0.98)
    
    # Save
    pdf_file = output_dir / "cyp_heatmap.pdf"
    plt.savefig(pdf_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {pdf_file}")
    
    png_file = output_dir / "cyp_heatmap.png"
    plt.savefig(png_file, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {png_file}")
    
    plt.close()
    
    print("\n✓ CYP heatmap complete!")


def main():
    parser = argparse.ArgumentParser(
        description="PyDESeq2 Step 3: Generate plots from filtered results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  # Generate all plots including CYP heatmap
  python pydeseq2_generate_plots.py filtered_results.tsv \\
      --count-matrix count_matrix.tsv \\
      --metadata metadata.tsv \\
      --cyp-family-map cyp_families.tsv \\
      --contrast-A root --contrast-B leaf \\
      -o plots/

  # MA and Volcano plots only (no heatmap)
  python pydeseq2_generate_plots.py filtered_results.tsv -o plots/
        """
    )
    
    parser.add_argument("filtered_results",
                       help="Path to filtered PyDESeq2 results TSV")
    
    parser.add_argument("-o", "--output", default="plots",
                       help="Output directory (default: plots)")
    
    parser.add_argument("--count-matrix", default=None,
                       help="Count matrix TSV (required for heatmaps)")
    
    parser.add_argument("--metadata", default=None,
                       help="Metadata TSV (required for heatmaps)")
    
    parser.add_argument("--cyp-family-map", default=None,
                       help="CYP family mapping TSV (for CYP heatmap)")
    
    parser.add_argument("--contrast-factor", default="condition",
                       help="Metadata column for contrast (default: condition)")
    parser.add_argument("--contrast-A", default="root",
                       help="Numerator condition (default: root)")
    parser.add_argument("--contrast-B", default="leaf",
                       help="Denominator condition (default: leaf)")
    
    parser.add_argument("--scale", choices=['center', 'zscore'], default='center',
                       help="Heatmap scaling method (default: center)")
    parser.add_argument("--no-row-cluster", action="store_true",
                       help="Disable row clustering in heatmap")
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.filtered_results):
        print(f"ERROR: Filtered results file not found: {args.filtered_results}")
        return 1
    
    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("PyDESeq2 Plot Generation (Step 3)")
    print("=" * 60)
    print(f"Filtered results: {args.filtered_results}")
    print(f"Output directory: {output_dir}")
    print()
    
    try:
        # Read filtered results
        filtered_results = pd.read_csv(args.filtered_results, sep='\t', index_col=0)
        print(f"Filtered genes: {len(filtered_results)}")
        
        # Generate MA and Volcano plots (always)
        generate_ma_volcano_plots(filtered_results, output_dir)
        
        # Generate heatmap if count matrix provided
        if args.count_matrix and args.metadata:
            if not os.path.exists(args.count_matrix):
                print(f"\nWARNING: Count matrix not found: {args.count_matrix}")
                print("Skipping heatmap generation.")
            elif not os.path.exists(args.metadata):
                print(f"\nWARNING: Metadata not found: {args.metadata}")
                print("Skipping heatmap generation.")
            else:
                print(f"\nReading count matrix: {args.count_matrix}")
                count_matrix = pd.read_csv(args.count_matrix, sep='\t', index_col=0)
                
                print(f"Reading metadata: {args.metadata}")
                metadata = pd.read_csv(args.metadata, sep='\t')
                if 'sample' in metadata.columns:
                    metadata = metadata.set_index('sample')
                
                # Generate CYP heatmap if mapping provided
                if args.cyp_family_map:
                    if os.path.exists(args.cyp_family_map):
                        cyp_family_map = load_cyp_family_map(args.cyp_family_map)
                        generate_cyp_heatmap(
                            filtered_results, count_matrix, metadata, cyp_family_map, output_dir,
                            contrast_A=args.contrast_A,
                            contrast_B=args.contrast_B,
                            contrast_factor=args.contrast_factor,
                            scale_method=args.scale,
                            row_cluster=not args.no_row_cluster
                        )
                    else:
                        print(f"\nWARNING: CYP family map not found: {args.cyp_family_map}")
                else:
                    print("\nNOTE: No CYP family map provided. Skipping CYP heatmap.")
                    print("      Use --cyp-family-map to generate CYP heatmap.")
        else:
            print("\nNOTE: Count matrix and/or metadata not provided.")
            print("      Only MA and Volcano plots generated.")
            print("      For heatmaps, provide --count-matrix and --metadata")
        
        print("\n" + "=" * 60)
        print("✓ Plot generation complete!")
        print("=" * 60)
        print(f"\nOutput directory: {output_dir}")
        
        return 0
        
    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
