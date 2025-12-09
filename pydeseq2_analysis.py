#!/usr/bin/env python3
"""
PyDESeq2 Differential Expression Analysis for RNA-seq Data

This script performs differential expression analysis using PyDESeq2,
a Python implementation of the DESeq2 method for bulk RNA-seq data.

Usage: python pydeseq2_analysis.py <count_matrix> <metadata> [output_dir] [--design DESIGN]
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
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    import matplotlib.pyplot as plt
    import seaborn as sns
except ImportError as e:
    print(f"ERROR: Missing required package: {e}")
    print("Please install PyDESeq2 and dependencies:")
    print("  pip install pydeseq2 matplotlib seaborn")
    sys.exit(1)


def read_data(count_file, metadata_file):
    """Read count matrix and metadata files."""
    print(f"Reading count matrix: {count_file}")
    count_matrix = pd.read_csv(count_file, sep='\t', index_col=0)
    
    print(f"Reading metadata: {metadata_file}")
    metadata = pd.read_csv(metadata_file, sep='\t')
    
    # Ensure sample names match
    common_samples = list(set(count_matrix.columns) & set(metadata['sample']))
    if len(common_samples) == 0:
        raise ValueError("No common samples found between count matrix and metadata")
    
    count_matrix = count_matrix[common_samples]
    metadata = metadata[metadata['sample'].isin(common_samples)].copy()
    metadata = metadata.set_index('sample')
    metadata = metadata.loc[common_samples]  # Ensure same order
    
    print(f"Samples: {len(common_samples)}")
    print(f"Genes: {len(count_matrix)}")
    
    return count_matrix, metadata


def generate_qc_plots(count_matrix, metadata, output_dir):
    """Generate quality control plots."""
    print("Generating quality control plots...")
    
    # Total read counts per sample
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    total_counts = count_matrix.sum()
    axes[0].bar(range(len(total_counts)), total_counts.values)
    axes[0].set_xticks(range(len(total_counts)))
    axes[0].set_xticklabels(total_counts.index, rotation=45, ha='right')
    axes[0].set_ylabel('Total Read Counts')
    axes[0].set_title('Total Read Counts per Sample')
    axes[0].grid(axis='y', alpha=0.3)
    
    axes[1].boxplot(total_counts.values)
    axes[1].set_ylabel('Total Read Counts')
    axes[1].set_title('Distribution of Total Counts')
    axes[1].grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / "qc_total_counts.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Gene count distribution
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    gene_totals = count_matrix.sum(axis=1)
    axes[0].hist(gene_totals, bins=50, edgecolor='black', alpha=0.7)
    axes[0].set_xlabel('Total Counts')
    axes[0].set_ylabel('Number of Genes')
    axes[0].set_title('Total Counts per Gene')
    axes[0].set_yscale('log')
    axes[0].grid(axis='y', alpha=0.3)
    
    axes[1].hist(np.log10(gene_totals + 1), bins=50, edgecolor='black', alpha=0.7)
    axes[1].set_xlabel('Log10(Total Counts + 1)')
    axes[1].set_ylabel('Number of Genes')
    axes[1].set_title('Log10 Total Counts per Gene')
    axes[1].grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / "qc_gene_counts.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    
    print("QC plots saved")


def run_pydeseq2_analysis(count_matrix, metadata, design_formula, output_dir):
    """Run PyDESeq2 differential expression analysis."""
    print(f"Running PyDESeq2 analysis with design: {design_formula}")
    
    # Create DeseqDataSet
    print("Creating DeseqDataSet...")
    # PyDESeq2 expects count matrix as (samples × genes), but we have (genes × samples)
    # Transpose the count matrix
    print(f"Count matrix shape (before transpose): {count_matrix.shape} (genes × samples)")
    count_matrix_t = count_matrix.T  # Transpose to (samples × genes)
    print(f"Count matrix shape (after transpose): {count_matrix_t.shape} (samples × genes)")
    
    # Ensure metadata index matches count matrix index (sample names)
    metadata_indexed = metadata.copy()
    if metadata_indexed.index.name != 'sample' and 'sample' in metadata_indexed.columns:
        metadata_indexed = metadata_indexed.set_index('sample')
    
    # Ensure same order as count matrix rows (which are now samples after transpose)
    metadata_indexed = metadata_indexed.loc[count_matrix_t.index]
    
    print(f"Metadata shape: {metadata_indexed.shape} (samples × variables)")
    print(f"Count matrix samples: {list(count_matrix_t.index[:5])}...")
    print(f"Metadata samples: {list(metadata_indexed.index[:5])}...")
    
    try:
        # Newer PyDESeq2 API uses 'counts' parameter
        # Count matrix should be (samples × genes)
        dds = DeseqDataSet(
            counts=count_matrix_t,
            metadata=metadata_indexed,
            design_factors=design_formula.split('+'),
            refit_cooks=True,
            n_cpus=1
        )
    except TypeError as e:
        # Try alternative parameter names
        error_msg = str(e)
        if 'count_matrix' in error_msg or 'counts' in error_msg:
            # Try with 'count_matrix' (older API)
            try:
                dds = DeseqDataSet(
                    count_matrix=count_matrix,
                    metadata=metadata_indexed,
                    design_factors=design_formula.split('+'),
                    refit_cooks=True,
                    n_cpus=1
                )
            except Exception as e2:
                raise ValueError(f"Failed to create DeseqDataSet. Error: {e2}\n"
                               f"Tried both 'counts' and 'count_matrix' parameters.")
        else:
            raise
    
    # Filter low count genes (similar to DESeq2 default)
    print("Filtering low count genes...")
    min_counts = 10
    # count_matrix_t is (samples × genes), so sum along axis=0 (columns) to get gene totals
    genes_to_keep = (count_matrix_t.sum(axis=0) >= min_counts)
    dds = dds[:, genes_to_keep]
    print(f"Genes after filtering (>= {min_counts} counts): {len(dds)}")
    
    # Fit dispersions and LFCs
    print("Fitting dispersions and LFCs...")
    dds.deseq2()
    
    # Get statistical test results
    print("Computing statistical tests...")
    
    # Determine contrast based on design formula
    if 'treatment' in design_formula:
        factor = 'treatment'
    elif 'group' in design_formula:
        factor = 'group'
    elif 'condition' in design_formula:
        factor = 'condition'
    else:
        # Use first factor in design
        factor = design_formula.split('+')[0].strip()
    
    # Get unique levels
    if factor in metadata.columns:
        levels = sorted(metadata[factor].unique())
        if len(levels) >= 2:
            contrast = (factor, levels[1], levels[0])
            print(f"Comparing {factor}: {levels[1]} vs {levels[0]}")
        else:
            contrast = None
            print(f"Only one level found for {factor}, using default contrast")
    else:
        contrast = None
        print(f"Factor {factor} not found in metadata, using default contrast")
    
    # Create DeseqStats object
    stat_res = DeseqStats(dds, contrast=contrast, n_cpus=1)
    stat_res.summary()
    
    # Get results
    results_df = stat_res.results_df
    
    # Save results
    results_file = output_dir / "pydeseq2_results.tsv"
    results_df.to_csv(results_file, sep='\t')
    print(f"Results saved: {results_file}")
    
    # Summary statistics
    sig_genes = results_df[
        (results_df['padj'].notna()) & 
        (results_df['padj'] < 0.05)
    ]
    
    print(f"\nSignificant genes (padj < 0.05): {len(sig_genes)}")
    if len(sig_genes) > 0:
        upreg = sig_genes[sig_genes['log2FoldChange'] > 0]
        downreg = sig_genes[sig_genes['log2FoldChange'] < 0]
        print(f"  Upregulated: {len(upreg)}")
        print(f"  Downregulated: {len(downreg)}")
    
    return dds, stat_res, results_df


def generate_plots(dds, stat_res, results_df, output_dir):
    """Generate visualization plots."""
    print("Generating visualization plots...")
    
    # MA plot
    print("  Creating MA plot...")
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Filter out infinite and NaN values
    valid = (results_df['log2FoldChange'].notna()) & (results_df['baseMean'].notna())
    valid = valid & (np.isfinite(results_df['log2FoldChange'])) & (np.isfinite(results_df['baseMean']))
    
    ax.scatter(
        results_df.loc[valid, 'baseMean'],
        results_df.loc[valid, 'log2FoldChange'],
        alpha=0.5,
        s=1,
        c='gray'
    )
    
    # Highlight significant genes
    sig_valid = valid & (results_df['padj'].notna()) & (results_df['padj'] < 0.05)
    if sig_valid.sum() > 0:
        ax.scatter(
            results_df.loc[sig_valid, 'baseMean'],
            results_df.loc[sig_valid, 'log2FoldChange'],
            alpha=0.7,
            s=2,
            c='red'
        )
    
    ax.axhline(y=0, color='black', linestyle='--', linewidth=0.5)
    ax.set_xlabel('Mean Normalized Counts')
    ax.set_ylabel('Log2 Fold Change')
    ax.set_title('PyDESeq2 MA Plot')
    ax.set_xscale('log')
    ax.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / "pydeseq2_ma_plot.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Volcano plot
    print("  Creating volcano plot...")
    fig, ax = plt.subplots(figsize=(10, 8))
    
    valid = (results_df['log2FoldChange'].notna()) & (results_df['padj'].notna())
    valid = valid & (np.isfinite(results_df['log2FoldChange'])) & (np.isfinite(results_df['padj']))
    
    # Convert p-values to -log10
    neg_log10_padj = -np.log10(results_df.loc[valid, 'padj'] + 1e-300)
    
    ax.scatter(
        results_df.loc[valid, 'log2FoldChange'],
        neg_log10_padj,
        alpha=0.5,
        s=1,
        c='gray'
    )
    
    # Highlight significant genes
    sig_valid = valid & (results_df['padj'] < 0.05)
    if sig_valid.sum() > 0:
        sig_neg_log10 = -np.log10(results_df.loc[sig_valid, 'padj'] + 1e-300)
        ax.scatter(
            results_df.loc[sig_valid, 'log2FoldChange'],
            sig_neg_log10,
            alpha=0.7,
            s=2,
            c='red'
        )
    
    ax.axhline(y=-np.log10(0.05), color='black', linestyle='--', linewidth=0.5, label='padj = 0.05')
    ax.axvline(x=1, color='black', linestyle='--', linewidth=0.5)
    ax.axvline(x=-1, color='black', linestyle='--', linewidth=0.5)
    ax.set_xlabel('Log2 Fold Change')
    ax.set_ylabel('-Log10 Adjusted P-value')
    ax.set_title('PyDESeq2 Volcano Plot')
    ax.legend()
    ax.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / "pydeseq2_volcano_plot.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Heatmap of top variable genes
    print("  Creating heatmap of top variable genes...")
    try:
        # Get normalized counts
        # dds.count_matrix is (samples × genes) after transpose
        normalized_counts = dds.count_matrix.div(dds.size_factors, axis=0)
        
        # Calculate variance across samples (axis=0) to find variable genes
        # normalized_counts is (samples × genes), so var(axis=0) gives variance per gene
        gene_vars = normalized_counts.var(axis=0)
        top_var_genes = gene_vars.nlargest(50).index
        
        # Prepare data for heatmap - select columns (genes) not rows
        heatmap_data = normalized_counts[top_var_genes]
        # Center by subtracting mean across samples (axis=0)
        heatmap_data = heatmap_data.subtract(heatmap_data.mean(axis=0), axis=1)
        
        # Create annotation
        if 'treatment' in dds.metadata.columns:
            annotation_col = dds.metadata[['treatment']]
        elif 'group' in dds.metadata.columns:
            annotation_col = dds.metadata[['group']]
        else:
            annotation_col = None
        
        fig, ax = plt.subplots(figsize=(10, 8))
        sns.heatmap(
            heatmap_data,
            cmap='RdBu_r',
            center=0,
            annot=False,
            fmt='.2f',
            cbar_kws={'label': 'Centered Log Counts'},
            yticklabels=False,
            xticklabels=True,
            ax=ax
        )
        ax.set_title('Top 50 Most Variable Genes')
        ax.set_xlabel('Samples')
        ax.set_ylabel('Genes')
        
        plt.tight_layout()
        plt.savefig(output_dir / "heatmap_top_variable_genes.pdf", dpi=300, bbox_inches='tight')
        plt.close()
    except Exception as e:
        print(f"  Warning: Could not create heatmap: {e}")
    
    print("Visualization plots saved")


def generate_summary(count_file, metadata_file, results_df, dds, output_dir):
    """Generate summary report."""
    print("Generating summary report...")
    
    summary_file = output_dir / "analysis_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("PyDESeq2 Differential Expression Analysis Summary\n")
        f.write("=" * 60 + "\n\n")
        f.write("Input files:\n")
        f.write(f"  Count matrix: {count_file}\n")
        f.write(f"  Metadata: {metadata_file}\n\n")
        
        f.write("Data summary:\n")
        f.write(f"  Total samples: {len(dds.count_matrix)}\n")
        f.write(f"  Total genes: {len(dds.count_matrix.columns)}\n")
        f.write(f"  Genes after filtering: {len(dds)}\n\n")
        
        f.write("PyDESeq2 results:\n")
        sig_genes = results_df[
            (results_df['padj'].notna()) & 
            (results_df['padj'] < 0.05)
        ]
        f.write(f"  Significant genes (padj < 0.05): {len(sig_genes)}\n")
        if len(sig_genes) > 0:
            upreg = sig_genes[sig_genes['log2FoldChange'] > 0]
            downreg = sig_genes[sig_genes['log2FoldChange'] < 0]
            f.write(f"  Upregulated genes: {len(upreg)}\n")
            f.write(f"  Downregulated genes: {len(downreg)}\n")
        
        f.write("\nOutput files:\n")
        f.write("  PyDESeq2 results: pydeseq2_results.tsv\n")
        f.write("  Quality control plots: qc_*.pdf\n")
        f.write("  Analysis plots: pydeseq2_*.pdf\n")
        f.write("  Heatmap: heatmap_top_variable_genes.pdf\n")
    
    print(f"Summary saved: {summary_file}")


def main():
    parser = argparse.ArgumentParser(
        description="PyDESeq2 differential expression analysis for RNA-seq data"
    )
    parser.add_argument("count_matrix", help="Path to gene count matrix TSV file")
    parser.add_argument("metadata", help="Path to sample metadata TSV file")
    parser.add_argument("-o", "--output", default="pydeseq2_results",
                       help="Output directory (default: pydeseq2_results)")
    parser.add_argument("--design", default=None,
                       help="Design formula (e.g., 'treatment' or 'group+condition'). "
                            "If not specified, will auto-detect from metadata")
    
    args = parser.parse_args()
    
    # Validate input files
    if not os.path.exists(args.count_matrix):
        print(f"ERROR: Count matrix file not found: {args.count_matrix}")
        return 1
    
    if not os.path.exists(args.metadata):
        print(f"ERROR: Metadata file not found: {args.metadata}")
        return 1
    
    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("PyDESeq2 Differential Expression Analysis")
    print("=" * 60)
    print(f"Count matrix: {args.count_matrix}")
    print(f"Metadata: {args.metadata}")
    print(f"Output directory: {output_dir}\n")
    
    try:
        # Read data
        count_matrix, metadata = read_data(args.count_matrix, args.metadata)
        
        # Determine design formula
        if args.design:
            design_formula = args.design
        else:
            # Auto-detect: prefer 'treatment', then 'group', then 'condition'
            if 'treatment' in metadata.columns:
                design_formula = 'treatment'
            elif 'group' in metadata.columns:
                design_formula = 'group'
            elif 'condition' in metadata.columns:
                design_formula = 'condition'
            else:
                # Use first categorical column
                categorical_cols = metadata.select_dtypes(include=['object', 'category']).columns
                if len(categorical_cols) > 0:
                    design_formula = categorical_cols[0]
                    print(f"Auto-detected design factor: {design_formula}")
                else:
                    raise ValueError("No suitable grouping variable found in metadata. "
                                   "Please specify --design or add 'treatment', 'group', or 'condition' column.")
        
        print(f"Using design formula: {design_formula}\n")
        
        # Generate QC plots
        generate_qc_plots(count_matrix, metadata, output_dir)
        
        # Run PyDESeq2 analysis
        dds, stat_res, results_df = run_pydeseq2_analysis(
            count_matrix, metadata, design_formula, output_dir
        )
        
        # Generate plots
        generate_plots(dds, stat_res, results_df, output_dir)
        
        # Generate summary
        generate_summary(args.count_matrix, args.metadata, results_df, dds, output_dir)
        
        print("\n" + "=" * 60)
        print("Analysis complete! Results saved to:", output_dir)
        print("=" * 60)
        
        return 0
        
    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())

