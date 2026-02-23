#!/usr/bin/env python3
"""
PyDESeq2 Differential Expression Analysis for RNA-seq Data

This script performs differential expression analysis using PyDESeq2,
a Python implementation of the DESeq2 method for bulk RNA-seq data.

Usage: python pydeseq2_countmatrix.py <count_matrix> <metadata> [output_dir] [--design DESIGN]
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


def run_pydeseq2_analysis(count_matrix, metadata, design_formula, output_dir,
                          contrast_A="R", contrast_B="L"):
    """Run PyDESeq2 differential expression analysis."""
    print(f"Running PyDESeq2 analysis with design: {design_formula}")
    
    print(f"Count matrix shape (before transpose): {count_matrix.shape} (genes x samples)")
    count_matrix_t = count_matrix.T
    print(f"Count matrix shape (after transpose): {count_matrix_t.shape} (samples x genes)")
    
    metadata_indexed = metadata.copy()
    if metadata_indexed.index.name != 'sample' and 'sample' in metadata_indexed.columns:
        metadata_indexed = metadata_indexed.set_index('sample')
    
    metadata_indexed = metadata_indexed.loc[count_matrix_t.index]
    
    print(f"Metadata shape: {metadata_indexed.shape} (samples x variables)")
    print(f"Count matrix samples: {list(count_matrix_t.index[:5])}...")
    print(f"Metadata samples: {list(metadata_indexed.index[:5])}...")
    
    try:
        dds = DeseqDataSet(
            counts=count_matrix_t,
            metadata=metadata_indexed,
            design_factors=design_formula.split('+'),
            refit_cooks=True,
            n_cpus=1
        )
    except TypeError as e:
        error_msg = str(e)
        if 'count_matrix' in error_msg or 'counts' in error_msg:
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
    
    print("Filtering low count genes...")
    min_counts = 10
    genes_to_keep = (count_matrix_t.sum(axis=0) >= min_counts)
    print(f"Genes before filtering: {len(count_matrix_t.columns)}")
    print(f"Genes with >= {min_counts} counts: {genes_to_keep.sum()}")
    
    count_matrix_filtered = count_matrix_t.loc[:, genes_to_keep]
    
    dds = DeseqDataSet(
        counts=count_matrix_filtered,
        metadata=metadata_indexed,
        design_factors=design_formula.split('+'),
        refit_cooks=True,
        n_cpus=1
    )
    
    count_matrix_for_plots = count_matrix_filtered.copy()
    
    print(f"Genes after filtering: {len(dds)}")
    
    if not isinstance(dds, DeseqDataSet):
        raise TypeError(f"Expected DeseqDataSet, but got {type(dds)}.")
    
    print("Fitting dispersions and LFCs...")
    
    try:
        if hasattr(dds, 'deseq2'):
            print("Using dds.deseq2() method...")
            dds.deseq2()
        elif hasattr(dds, 'fit'):
            print("Using dds.fit() method...")
            dds.fit()
        elif hasattr(dds, 'estimate_size_factors'):
            print("Using step-by-step workflow...")
            dds.estimate_size_factors()
            dds.estimate_dispersions()
            dds.fit()
        elif hasattr(dds, 'fit_size_factors'):
            print("Using alternative step-by-step workflow...")
            dds.fit_size_factors()
            dds.fit_genewise_dispersions()
            dds.fit_dispersion_trend()
            dds.fit_dispersion_prior()
            dds.fit_MAP_dispersions()
            dds.fit_LFC()
        else:
            available_methods = [m for m in dir(dds) if not m.startswith('_') and callable(getattr(dds, m, None))]
            fitting_methods = [m for m in available_methods if any(x in m.lower() for x in ['fit', 'deseq', 'estimate', 'size', 'dispersion'])]
            print(f"Available fitting-related methods: {fitting_methods}")
            raise AttributeError("No recognized fitting method found.")
    except AttributeError as e:
        print(f"ERROR: {e}")
        raise
    
    print("Computing statistical tests...")
    
    factor = design_formula.split('+')[0].strip()
    
    # contrast_A = numerator (root), contrast_B = denominator/baseline (leaf)
    #   log2FC = log2( A / B ) = log2( root / leaf )
    #   Positive log2FC → gene is HIGHER in root
    #   Negative log2FC → gene is HIGHER in leaf
    contrast = [factor, contrast_A, contrast_B]
    print(f"Contrast: {factor} — {contrast_A} vs {contrast_B} (baseline)")
    print(f"  Positive log2FC = higher in {contrast_A} (root)")
    print(f"  Negative log2FC = higher in {contrast_B} (leaf)")
    
    stat_res = DeseqStats(dds, contrast=contrast, n_cpus=1)
    stat_res.summary()
    
    results_df = stat_res.results_df
    
    results_file = output_dir / "pydeseq2_results.tsv"
    results_df.to_csv(results_file, sep='\t')
    print(f"Results saved: {results_file}")
    
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
    
    return dds, stat_res, results_df, count_matrix_for_plots


def generate_plots(dds, stat_res, results_df, output_dir, count_matrix_for_plots=None):
    """Generate visualization plots."""
    print("Generating visualization plots...")
    
    # MA plot
    print("  Creating MA plot...")
    fig, ax = plt.subplots(figsize=(8, 6))
    
    valid = (results_df['log2FoldChange'].notna()) & (results_df['baseMean'].notna())
    valid = valid & (np.isfinite(results_df['log2FoldChange'])) & (np.isfinite(results_df['baseMean']))
    
    ax.scatter(
        results_df.loc[valid, 'baseMean'],
        results_df.loc[valid, 'log2FoldChange'],
        alpha=0.5, s=1, c='gray'
    )
    
    sig_valid = valid & (results_df['padj'].notna()) & (results_df['padj'] < 0.05)
    if sig_valid.sum() > 0:
        ax.scatter(
            results_df.loc[sig_valid, 'baseMean'],
            results_df.loc[sig_valid, 'log2FoldChange'],
            alpha=0.7, s=2, c='red'
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
    
    neg_log10_padj = -np.log10(results_df.loc[valid, 'padj'] + 1e-300)
    
    ax.scatter(
        results_df.loc[valid, 'log2FoldChange'],
        neg_log10_padj,
        alpha=0.5, s=1, c='gray'
    )
    
    sig_valid = valid & (results_df['padj'] < 0.05)
    if sig_valid.sum() > 0:
        sig_neg_log10 = -np.log10(results_df.loc[sig_valid, 'padj'] + 1e-300)
        ax.scatter(
            results_df.loc[sig_valid, 'log2FoldChange'],
            sig_neg_log10,
            alpha=0.7, s=2, c='red'
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
        if count_matrix_for_plots is None or not isinstance(count_matrix_for_plots, pd.DataFrame):
            raise ValueError("Count matrix for plots is not available or not a DataFrame")
        
        count_data = count_matrix_for_plots.copy()
        
        if hasattr(dds, 'size_factors') and dds.size_factors is not None:
            normalized_counts = count_data.div(dds.size_factors, axis=0)
        elif hasattr(dds, 'obs') and 'size_factors' in dds.obs.columns:
            normalized_counts = count_data.div(dds.obs['size_factors'], axis=0)
        else:
            normalized_counts = count_data
        
        gene_vars = normalized_counts.var(axis=0)
        top_var_genes = gene_vars.nlargest(50).index.tolist()
        
        heatmap_data = normalized_counts[top_var_genes]
        heatmap_data = heatmap_data.T
        heatmap_data = heatmap_data.subtract(heatmap_data.mean(axis=1), axis=0)
        
        metadata_df = None
        if hasattr(dds, 'obs'):
            metadata_df = dds.obs
        elif hasattr(dds, 'metadata'):
            metadata_df = dds.metadata
        
        annotation_col = None
        if metadata_df is not None and len(metadata_df) > 0:
            try:
                metadata_aligned = metadata_df.reindex(heatmap_data.columns)
                
                if 'treatment' in metadata_aligned.columns:
                    annotation_col = metadata_aligned[['treatment']]
                elif 'group' in metadata_aligned.columns:
                    annotation_col = metadata_aligned[['group']]
                elif 'condition' in metadata_aligned.columns:
                    annotation_col = metadata_aligned[['condition']]
            except Exception as e:
                print(f"  Warning: Could not align metadata: {e}")
                annotation_col = None
        
        fig, ax = plt.subplots(figsize=(12, 10))
        
        heatmap_kwargs = {
            'data': heatmap_data,
            'cmap': 'RdBu_r',
            'center': 0,
            'annot': False,
            'fmt': '.2f',
            'cbar_kws': {'label': 'Centered Normalized Counts'},
            'yticklabels': False,
            'xticklabels': True,
            'ax': ax
        }
        
        if annotation_col is not None:
            try:
                heatmap_kwargs['col_colors'] = annotation_col
            except Exception as e:
                print(f"  Warning: Could not add annotations: {e}")
        
        sns.heatmap(**heatmap_kwargs)
        
        ax.set_title('Top 50 Most Variable Genes', fontsize=14, fontweight='bold')
        ax.set_xlabel('Samples', fontsize=12)
        ax.set_ylabel('Genes', fontsize=12)
        
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
        if hasattr(dds, 'count_matrix'):
            n_samples = len(dds.count_matrix)
            n_genes = len(dds.count_matrix.columns)
        elif hasattr(dds, 'X'):
            n_samples = dds.X.shape[0]
            n_genes = dds.X.shape[1]
        elif hasattr(dds, 'counts'):
            if hasattr(dds.counts, 'shape'):
                n_samples = dds.counts.shape[0]
                n_genes = dds.counts.shape[1]
            else:
                n_samples = len(dds.counts)
                n_genes = len(dds.counts.columns) if hasattr(dds.counts, 'columns') else 0
        else:
            n_samples = len(dds) if hasattr(dds, '__len__') else 0
            n_genes = len(dds.var) if hasattr(dds, 'var') else 0
        
        f.write(f"  Total samples: {n_samples}\n")
        f.write(f"  Total genes: {n_genes}\n")
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
    parser.add_argument("--contrast-A", default="R",
                       help="Numerator condition - positive log2FC means higher in A (default: R = root)")
    parser.add_argument("--contrast-B", default="L",
                       help="Denominator/baseline condition (default: L = leaf)")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.count_matrix):
        print(f"ERROR: Count matrix file not found: {args.count_matrix}")
        return 1
    
    if not os.path.exists(args.metadata):
        print(f"ERROR: Metadata file not found: {args.metadata}")
        return 1
    
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("PyDESeq2 Differential Expression Analysis")
    print("=" * 60)
    print(f"Count matrix: {args.count_matrix}")
    print(f"Metadata: {args.metadata}")
    print(f"Output directory: {output_dir}\n")
    
    try:
        count_matrix, metadata = read_data(args.count_matrix, args.metadata)
        
        if args.design:
            design_formula = args.design
        else:
            if 'treatment' in metadata.columns:
                design_formula = 'treatment'
            elif 'group' in metadata.columns:
                design_formula = 'group'
            elif 'condition' in metadata.columns:
                design_formula = 'condition'
            else:
                categorical_cols = metadata.select_dtypes(include=['object', 'category']).columns
                if len(categorical_cols) > 0:
                    design_formula = categorical_cols[0]
                    print(f"Auto-detected design factor: {design_formula}")
                else:
                    raise ValueError("No suitable grouping variable found in metadata. "
                                   "Please specify --design or add 'treatment', 'group', or 'condition' column.")
        
        print(f"Using design formula: {design_formula}\n")
        
        generate_qc_plots(count_matrix, metadata, output_dir)
        
        dds, stat_res, results_df, count_matrix_for_plots = run_pydeseq2_analysis(
            count_matrix, metadata, design_formula, output_dir,
            contrast_A=args.contrast_A, contrast_B=args.contrast_B
        )
        
        generate_plots(dds, stat_res, results_df, output_dir, count_matrix_for_plots)
        
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
