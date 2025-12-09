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
    print(f"Genes before filtering: {len(count_matrix_t.columns)}")
    print(f"Genes with >= {min_counts} counts: {genes_to_keep.sum()}")
    
    # Filter before creating DeseqDataSet
    count_matrix_filtered = count_matrix_t.loc[:, genes_to_keep]
    print(f"Creating DeseqDataSet with filtered counts...")
    
    # Recreate DeseqDataSet with filtered counts
    print("Creating DeseqDataSet object...")
    dds = DeseqDataSet(
        counts=count_matrix_filtered,
        metadata=metadata_indexed,
        design_factors=design_formula.split('+'),
        refit_cooks=True,
        n_cpus=1
    )
    
    # Store reference to count matrix for later use (since dds.count_matrix may not exist)
    # In PyDESeq2, counts are stored in dds.X (AnnData format) or we can use the original
    count_matrix_for_plots = count_matrix_filtered.copy()
    
    print(f"Genes after filtering: {len(dds)}")
    print(f"DeseqDataSet type: {type(dds)}")
    print(f"Is instance of DeseqDataSet: {isinstance(dds, DeseqDataSet)}")
    
    # Verify the object has the expected methods
    if not isinstance(dds, DeseqDataSet):
        raise TypeError(f"Expected DeseqDataSet, but got {type(dds)}. "
                       f"Please check PyDESeq2 installation and version.")
    
    # Fit dispersions and LFCs
    print("Fitting dispersions and LFCs...")
    
    # PyDESeq2 workflow: Try different API versions
    # Some versions have dds.deseq2(), others use step-by-step methods
    try:
        if hasattr(dds, 'deseq2'):
            print("Using dds.deseq2() method...")
            dds.deseq2()
        elif hasattr(dds, 'fit'):
            print("Using dds.fit() method...")
            dds.fit()
        elif hasattr(dds, 'estimate_size_factors'):
            # Step-by-step workflow (newer API)
            print("Using step-by-step workflow (estimate_size_factors, estimate_dispersions, fit)...")
            dds.estimate_size_factors()
            dds.estimate_dispersions()
            dds.fit()
        elif hasattr(dds, 'fit_size_factors'):
            # Alternative step-by-step workflow
            print("Using alternative step-by-step workflow...")
            dds.fit_size_factors()
            dds.fit_genewise_dispersions()
            dds.fit_dispersion_trend()
            dds.fit_dispersion_prior()
            dds.fit_MAP_dispersions()
            dds.fit_LFC()
        else:
            # Check available methods for debugging
            available_methods = [m for m in dir(dds) if not m.startswith('_') and callable(getattr(dds, m, None))]
            fitting_methods = [m for m in available_methods if any(x in m.lower() for x in ['fit', 'deseq', 'estimate', 'size', 'dispersion'])]
            print(f"Available fitting-related methods: {fitting_methods}")
            raise AttributeError("No recognized fitting method found. Please check PyDESeq2 version.")
    except AttributeError as e:
        print(f"ERROR: {e}")
        print(f"DeseqDataSet type: {type(dds)}")
        print(f"Is instance of DeseqDataSet: {isinstance(dds, DeseqDataSet)}")
        raise
    
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
    
    return dds, stat_res, results_df, count_matrix_for_plots


def generate_plots(dds, stat_res, results_df, output_dir, count_matrix_for_plots=None):
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
        # Try different ways to access count matrix from DeseqDataSet
        if hasattr(dds, 'count_matrix'):
            count_data = dds.count_matrix
        elif hasattr(dds, 'X'):
            # AnnData stores data in .X
            count_data = pd.DataFrame(dds.X, index=dds.obs.index, columns=dds.var.index)
        elif hasattr(dds, 'counts'):
            count_data = dds.counts
        elif count_matrix_for_plots is not None:
            # Use the original count matrix we stored
            count_data = count_matrix_for_plots
        else:
            # Last resort: try to reconstruct from dds
            raise ValueError("Could not access count matrix from DeseqDataSet. "
                          "Please ensure PyDESeq2 is properly installed.")
        
        # Normalize by size factors
        # count_data is (samples × genes)
        if hasattr(dds, 'size_factors') and dds.size_factors is not None:
            # size_factors should be a Series or array with same length as samples
            if isinstance(count_data, pd.DataFrame):
                normalized_counts = count_data.div(dds.size_factors, axis=0)
            else:
                # If count_data is numpy array, convert to DataFrame first
                normalized_counts = pd.DataFrame(count_data, index=dds.obs.index, columns=dds.var.index)
                normalized_counts = normalized_counts.div(dds.size_factors, axis=0)
        elif hasattr(dds, 'obs') and 'size_factors' in dds.obs.columns:
            # Size factors might be stored in obs
            normalized_counts = count_data.div(dds.obs['size_factors'], axis=0)
        else:
            # If no size factors, use raw counts (log transform for visualization)
            print("  Warning: No size factors found, using raw counts")
            normalized_counts = count_data
        
        # Calculate variance across samples (axis=0) to find variable genes
        # normalized_counts is (samples × genes), so var(axis=0) gives variance per gene
        gene_vars = normalized_counts.var(axis=0)
        top_var_genes = gene_vars.nlargest(50).index
        
        # Prepare data for heatmap - select columns (genes) not rows
        # normalized_counts is (samples × genes), so we select genes (columns)
        heatmap_data = normalized_counts[top_var_genes]
        
        # Transpose for visualization: genes as rows, samples as columns
        heatmap_data = heatmap_data.T
        
        # Center by subtracting mean across samples (axis=1 now since transposed)
        heatmap_data = heatmap_data.subtract(heatmap_data.mean(axis=1), axis=0)
        
        # Create annotation for samples (columns)
        # In AnnData/PyDESeq2, metadata is stored in dds.obs, not dds.metadata
        metadata_df = dds.obs if hasattr(dds, 'obs') else (dds.metadata if hasattr(dds, 'metadata') else None)
        
        annotation_col = None
        if metadata_df is not None:
            # Ensure metadata index matches heatmap column names (sample names)
            if metadata_df.index.equals(heatmap_data.columns):
                if 'treatment' in metadata_df.columns:
                    annotation_col = metadata_df[['treatment']]
                elif 'group' in metadata_df.columns:
                    annotation_col = metadata_df[['group']]
                elif 'condition' in metadata_df.columns:
                    annotation_col = metadata_df[['condition']]
            else:
                # Reindex to match heatmap columns
                metadata_aligned = metadata_df.reindex(heatmap_data.columns)
                if 'treatment' in metadata_aligned.columns:
                    annotation_col = metadata_aligned[['treatment']]
                elif 'group' in metadata_aligned.columns:
                    annotation_col = metadata_aligned[['group']]
                elif 'condition' in metadata_aligned.columns:
                    annotation_col = metadata_aligned[['condition']]
        
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # Create heatmap with optional annotation
        heatmap_kwargs = {
            'data': heatmap_data,
            'cmap': 'RdBu_r',
            'center': 0,
            'annot': False,
            'fmt': '.2f',
            'cbar_kws': {'label': 'Centered Normalized Counts'},
            'yticklabels': True,
            'xticklabels': True,
            'ax': ax
        }
        
        if annotation_col is not None:
            heatmap_kwargs['col_colors'] = annotation_col
        
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
        # Get count matrix dimensions - try multiple ways to access
        if hasattr(dds, 'count_matrix'):
            n_samples = len(dds.count_matrix)
            n_genes = len(dds.count_matrix.columns)
        elif hasattr(dds, 'X'):
            # AnnData stores data in .X (samples × genes)
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
            # Fallback: use shape from dds itself
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
        dds, stat_res, results_df, count_matrix_for_plots = run_pydeseq2_analysis(
            count_matrix, metadata, design_formula, output_dir
        )
        
        # Generate plots
        generate_plots(dds, stat_res, results_df, output_dir, count_matrix_for_plots)
        
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

