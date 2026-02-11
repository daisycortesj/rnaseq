#!/usr/bin/env python3
"""
PyDESeq2 Differential Expression Analysis for RNA-seq Data
===========================================================

This script performs differential expression analysis using PyDESeq2,
a Python implementation of the DESeq2 method for bulk RNA-seq data.

Additionally, it generates a CYP gene expression heatmap matching
Figure 6A from the Project Narrative (CYP genes upregulated in root vs leaf).

=== HOW TO MATCH FIGURE 6A (CYP Heatmap) ===

Required inputs:
  1. count_matrix.tsv: Gene count matrix (genes x samples)
  2. metadata.tsv:     Sample metadata with 'sample' and 'condition' columns
  3. cyp_family_map.tsv: CYP gene family mapping (gene_id, cyp_family, cyp_clan)

Minimal command for Figure 6A:
  python pydeseq2_analysis.py \\
      count_matrix.tsv \\
      metadata.tsv \\
      -o results \\
      --contrast-factor condition \\
      --contrast-A root \\
      --contrast-B leaf \\
      --cyp-family-map cyp_family_map.tsv \\
      --root-up-only \\
      --lfc 2.0 \\
      --scale center

Output files:
  - pydeseq2_results.tsv:           Full DE results
  - deg_filtered.tsv:               DEGs passing padj/lfc filters
  - cyp_deg_filtered.tsv:           CYP DEGs only
  - cyp_heatmap_matrix.tsv:         Transformed values used for heatmap
  - cyp_heatmap.pdf / cyp_heatmap.png: Figure 6A-style heatmap

=== CYP Family Map Format ===

Tab-separated file with at minimum:
  gene_id    cyp_family    [cyp_clan]
  LOC1234    CYP71         CYP71_clan
  LOC5678    CYP72         CYP72_clan
  ...

The cyp_family column is used for grouping and color-coding rows.
The cyp_clan column (if present) can provide additional grouping.

Usage: python pydeseq2_analysis.py <count_matrix> <metadata> [options]
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import re
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

try:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.patches import Patch
    from scipy.cluster.hierarchy import linkage, leaves_list
    from scipy.spatial.distance import pdist
except ImportError as e:
    print(f"ERROR: Missing required package: {e}")
    print("Please install PyDESeq2 and dependencies:")
    print("  pip install pydeseq2 matplotlib seaborn scipy")
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


def load_cyp_family_map(filepath):
    """
    Load CYP gene family mapping from TSV file.
    
    Expected columns: gene_id, cyp_family (required), cyp_clan (optional)
    Returns a DataFrame indexed by gene_id.
    """
    print(f"Loading CYP family map: {filepath}")
    df = pd.read_csv(filepath, sep='\t')
    
    # Check required columns
    if 'gene_id' not in df.columns:
        raise ValueError("CYP family map must have 'gene_id' column")
    if 'cyp_family' not in df.columns:
        raise ValueError("CYP family map must have 'cyp_family' column")
    
    df = df.set_index('gene_id')
    print(f"  Loaded {len(df)} CYP gene mappings")
    print(f"  Families: {df['cyp_family'].nunique()} unique")
    
    return df


def load_cyp_gene_list(filepath):
    """
    Load simple CYP gene list (one gene ID per line).
    Returns a set of gene IDs.
    """
    print(f"Loading CYP gene list: {filepath}")
    with open(filepath, 'r') as f:
        genes = set(line.strip() for line in f if line.strip())
    print(f"  Loaded {len(genes)} CYP gene IDs")
    return genes


def load_sample_order(filepath):
    """
    Load explicit sample ordering from file (one sample name per line).
    Returns a list of sample names.
    """
    print(f"Loading sample order: {filepath}")
    with open(filepath, 'r') as f:
        samples = [line.strip() for line in f if line.strip()]
    print(f"  Loaded order for {len(samples)} samples")
    return samples


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
                          contrast_factor=None, contrast_A=None, contrast_B=None):
    """
    Run PyDESeq2 differential expression analysis.
    
    Args:
        contrast_factor: Factor column for contrast (e.g., 'condition')
        contrast_A: Numerator level (e.g., 'root') - positive log2FC means A > B
        contrast_B: Denominator level (e.g., 'leaf')
    
    Returns:
        dds, stat_res, results_df, count_matrix_for_plots
    """
    print(f"Running PyDESeq2 analysis with design: {design_formula}")
    
    # Create DeseqDataSet
    print("Creating DeseqDataSet...")
    print(f"Count matrix shape (before transpose): {count_matrix.shape} (genes x samples)")
    count_matrix_t = count_matrix.T  # Transpose to (samples x genes)
    print(f"Count matrix shape (after transpose): {count_matrix_t.shape} (samples x genes)")
    
    # Ensure metadata index matches count matrix index (sample names)
    metadata_indexed = metadata.copy()
    if metadata_indexed.index.name != 'sample' and 'sample' in metadata_indexed.columns:
        metadata_indexed = metadata_indexed.set_index('sample')
    
    # Ensure same order as count matrix rows (which are now samples after transpose)
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
                raise ValueError(f"Failed to create DeseqDataSet. Error: {e2}")
        else:
            raise
    
    # Filter low count genes
    print("Filtering low count genes...")
    min_counts = 10
    genes_to_keep = (count_matrix_t.sum(axis=0) >= min_counts)
    print(f"Genes before filtering: {len(count_matrix_t.columns)}")
    print(f"Genes with >= {min_counts} counts: {genes_to_keep.sum()}")
    
    count_matrix_filtered = count_matrix_t.loc[:, genes_to_keep]
    print("Creating DeseqDataSet with filtered counts...")
    
    ## This flag tells PyDESeq2 to check for outliers  - refit_cooks=True
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
    
    # Fit dispersions and LFCs
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
    
    # Get statistical test results
    print("Computing statistical tests...")
    
    # Determine contrast
    if contrast_factor and contrast_A and contrast_B:
        factor = contrast_factor
        if factor not in dds.obs.columns:
            raise ValueError(f"Contrast factor '{factor}' not found in metadata. "
                           f"Available columns: {list(dds.obs.columns)}")
        
        available_levels = dds.obs[factor].unique()
        if contrast_A not in available_levels:
            raise ValueError(f"Contrast level A '{contrast_A}' not found for factor '{factor}'. "
                           f"Available levels: {list(available_levels)}")
        if contrast_B not in available_levels:
            raise ValueError(f"Contrast level B '{contrast_B}' not found for factor '{factor}'. "
                           f"Available levels: {list(available_levels)}")
        
        contrast = (factor, contrast_A, contrast_B)
        print(f"Using specified contrast: {factor}: {contrast_A} vs {contrast_B}")
        print(f"  (Positive log2FC means higher in {contrast_A})")
    else:
        # Auto-detect contrast
        if 'treatment' in design_formula:
            factor = 'treatment'
        elif 'group' in design_formula:
            factor = 'group'
        elif 'condition' in design_formula:
            factor = 'condition'
        else:
            factor = design_formula.split('+')[0].strip()
        
        if factor in dds.obs.columns:
            levels = sorted(dds.obs[factor].unique())
            if len(levels) >= 2:
                contrast = (factor, levels[1], levels[0])
                print(f"Auto-detected contrast: {factor}: {levels[1]} vs {levels[0]}")
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


def generate_ma_volcano_plots(results_df, output_dir):
    """Generate MA and Volcano plots."""
    print("Generating MA and Volcano plots...")
    
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
    
    print("  MA and Volcano plots saved")


def get_normalized_counts(dds, count_matrix_for_plots):
    """
    Extract normalized counts from PyDESeq2 using size factors.
    Returns DataFrame with genes as rows, samples as columns.
    """
    count_data = count_matrix_for_plots.copy()
    
    if hasattr(dds, 'size_factors') and dds.size_factors is not None:
        normalized_counts = count_data.div(dds.size_factors, axis=0)
        print("  Normalized counts using size factors")
    elif hasattr(dds, 'obs') and 'size_factors' in dds.obs.columns:
        normalized_counts = count_data.div(dds.obs['size_factors'], axis=0)
        print("  Normalized counts using size factors from obs")
    else:
        print("  Warning: No size factors found, using raw counts")
        normalized_counts = count_data
    
    return normalized_counts.T


def generate_cyp_heatmap(dds, results_df, output_dir, count_matrix_for_plots,
                         cyp_family_map=None, cyp_genes=None,
                         padj_cutoff=0.05, lfc_cutoff=2.0, root_up_only=True,
                         sample_order=None, scale_method='center',
                         row_cluster=True, col_cluster=False,
                         contrast_A='root', contrast_B='leaf', contrast_factor='condition'):
    """
    Generate CYP gene expression heatmap (Figure 6A style).
    """
    print("\n" + "=" * 60)
    print("Generating CYP Gene Expression Heatmap (Figure 6A style)")
    print("=" * 60)
    
    # Step 1: Filter DEGs by padj and log2FC
    print("\nStep 1: Filtering DEGs...")
    
    deg_df = results_df.copy()
    n_total = len(deg_df)
    #   # Get number of genes with padj < padj_cutoff
    # OLD CODE: deg_df = deg_df[deg_df['padj'].notna() & (deg_df['padj'] < padj_cutoff)]
    #   n_padj = len(deg_df)
    #   print(f"  Genes with padj < {padj_cutoff}: {n_padj} (from {n_total})")


    # Remove genes with missing padj values
    ## This is QC, not statistical filtering 
    deg_df = deg_df[deg_df['padj'].notna()]
    
    # Apply padj filter (unless disabled with padj >= 1.0)
    if padj_cutoff < 1.0:
        deg_df = deg_df[deg_df['padj'] < padj_cutoff]
        n_padj = len(deg_df)
        print(f"  Genes with padj < {padj_cutoff}: {n_padj} (from {n_total})")
    else:
        print(f"  padj filter DISABLED (padj_cutoff={padj_cutoff}): keeping all {len(deg_df)} genes")
    
    # Apply log2FC filter (unless disabled with lfc_cutoff = 0)
    if lfc_cutoff > 0:
        deg_df = deg_df[np.abs(deg_df['log2FoldChange']) > lfc_cutoff]
        n_lfc = len(deg_df)
        print(f"  Genes with |log2FC| > {lfc_cutoff}: {n_lfc}")
    else:
        print(f"  log2FC filter DISABLED (lfc_cutoff=0): keeping all {len(deg_df)} genes")
    
    if root_up_only:
        deg_df = deg_df[deg_df['log2FoldChange'] > 0]
        n_up = len(deg_df)
        print(f"  Genes upregulated in {contrast_A}: {n_up}")
    
    if len(deg_df) == 0:
        print("  WARNING: No DEGs passed filters. Cannot create CYP heatmap.")
        return
    
    deg_file = output_dir / "deg_filtered.tsv"
    deg_df.to_csv(deg_file, sep='\t')
    print(f"  Saved filtered DEGs: {deg_file}")
    
    # Step 2: Filter to CYP genes
    print("\nStep 2: Filtering to CYP genes...")
    
    deg_gene_ids = set(deg_df.index)
    
    if cyp_family_map is not None:
        cyp_gene_ids = set(cyp_family_map.index)
        cyp_deg_ids = deg_gene_ids & cyp_gene_ids
        
        if len(cyp_deg_ids) == 0:
            print("  WARNING: No exact matches between DEG IDs and CYP map IDs.")
            print(f"    Sample DEG IDs: {list(deg_gene_ids)[:5]}")
            print(f"    Sample CYP IDs: {list(cyp_gene_ids)[:5]}")
            print("  Cannot create CYP heatmap without matching genes.")
            return
        
        print(f"  CYP DEGs found: {len(cyp_deg_ids)} (from {len(deg_gene_ids)} DEGs, {len(cyp_gene_ids)} CYP genes)")
        
        cyp_deg_df = deg_df.loc[list(cyp_deg_ids)].copy()
        cyp_deg_df['cyp_family'] = cyp_family_map.loc[cyp_deg_df.index, 'cyp_family']
        if 'cyp_clan' in cyp_family_map.columns:
            cyp_deg_df['cyp_clan'] = cyp_family_map.loc[cyp_deg_df.index, 'cyp_clan']
        
    elif cyp_genes is not None:
        cyp_deg_ids = deg_gene_ids & cyp_genes
        
        if len(cyp_deg_ids) == 0:
            print("  WARNING: No CYP genes found in DEG list.")
            print(f"    Sample DEG IDs: {list(deg_gene_ids)[:5]}")
            print(f"    Sample CYP IDs: {list(cyp_genes)[:5]}")
            print("  Cannot create CYP heatmap without matching genes.")
            return
        
        print(f"  CYP DEGs found: {len(cyp_deg_ids)}")
        cyp_deg_df = deg_df.loc[list(cyp_deg_ids)].copy()
        cyp_deg_df['cyp_family'] = 'CYP'
        
    else:
        print("  ERROR: No CYP gene mapping provided (--cyp-family-map or --cyp-genes)")
        print("  Cannot create CYP heatmap.")
        return
    
    cyp_deg_file = output_dir / "cyp_deg_filtered.tsv"
    cyp_deg_df.to_csv(cyp_deg_file, sep='\t')
    print(f"  Saved CYP DEGs: {cyp_deg_file}")
    
    # Step 3: Get normalized counts for CYP genes
    print("\nStep 3: Extracting normalized counts...")
    
    norm_counts = get_normalized_counts(dds, count_matrix_for_plots)
    
    cyp_genes_in_matrix = [g for g in cyp_deg_df.index if g in norm_counts.index]
    if len(cyp_genes_in_matrix) < len(cyp_deg_df):
        print(f"  WARNING: {len(cyp_deg_df) - len(cyp_genes_in_matrix)} CYP genes not in count matrix")
    
    cyp_norm_counts = norm_counts.loc[cyp_genes_in_matrix]
    print(f"  Normalized counts shape: {cyp_norm_counts.shape} (genes x samples)")
    
    # Step 4: Transform counts: log2(normalized + 1)
    print("\nStep 4: Log2 transforming counts...")
    
    heatmap_data = np.log2(cyp_norm_counts + 1)
    print("  Transformed to log2(normalized + 1)")
    
    # Step 5: Scale per gene (row)
    print(f"\nStep 5: Scaling per gene (method: {scale_method})...")
    
    if scale_method == 'center':
        row_means = heatmap_data.mean(axis=1)
        heatmap_data = heatmap_data.subtract(row_means, axis=0)
        scale_label = "log2(norm+1), centered"
    elif scale_method == 'zscore':
        row_means = heatmap_data.mean(axis=1)
        row_stds = heatmap_data.std(axis=1)
        row_stds = row_stds.replace(0, 1)
        heatmap_data = heatmap_data.subtract(row_means, axis=0).div(row_stds, axis=0)
        scale_label = "log2(norm+1), z-scored"
    else:
        raise ValueError(f"Unknown scale method: {scale_method}")
    
    print("  Scaling complete")
    
    # Step 6: Order columns (samples)
    print("\nStep 6: Ordering columns (samples)...")
    
    if hasattr(dds, 'obs'):
        metadata_df = dds.obs
    elif hasattr(dds, 'metadata'):
        metadata_df = dds.metadata
    else:
        metadata_df = None
        print("  WARNING: No metadata found for sample annotation")
    
    if sample_order is not None:
        valid_samples = [s for s in sample_order if s in heatmap_data.columns]
        missing_samples = [s for s in sample_order if s not in heatmap_data.columns]
        extra_samples = [s for s in heatmap_data.columns if s not in sample_order]
        
        if missing_samples:
            print(f"  WARNING: Samples in order file but not in data: {missing_samples}")
        if extra_samples:
            print(f"  WARNING: Samples in data but not in order file: {extra_samples}")
        
        column_order = valid_samples + extra_samples
        print(f"  Using explicit sample order ({len(column_order)} samples)")
    elif metadata_df is not None and contrast_factor in metadata_df.columns:
        condition_col = metadata_df[contrast_factor]
        
        samples_B = [s for s in heatmap_data.columns if s in condition_col.index and condition_col[s] == contrast_B]
        samples_A = [s for s in heatmap_data.columns if s in condition_col.index and condition_col[s] == contrast_A]
        other_samples = [s for s in heatmap_data.columns if s not in samples_B + samples_A]
        
        column_order = samples_B + samples_A + other_samples
        print(f"  Ordering by condition: {contrast_B} ({len(samples_B)}) | {contrast_A} ({len(samples_A)})")
    else:
        column_order = list(heatmap_data.columns)
        print("  Using default column order")
    
    heatmap_data = heatmap_data[column_order]
    
    # Step 7: Order rows by CYP family
    print("\nStep 7: Ordering rows by CYP family...")
    
    gene_families = cyp_deg_df.loc[heatmap_data.index, 'cyp_family']
    unique_families = sorted(gene_families.unique())
    
    print(f"  CYP families present: {unique_families}")
    
    row_order = []
    family_blocks = {}
    
    for family in unique_families:
        family_genes = gene_families[gene_families == family].index.tolist()
        
        if len(family_genes) > 1 and row_cluster:
            family_data = heatmap_data.loc[family_genes]
            try:
                linkage_matrix = linkage(family_data, method='ward', metric='euclidean')
                ordered_idx = leaves_list(linkage_matrix)
                family_genes = [family_genes[i] for i in ordered_idx]
            except Exception as e:
                print(f"    Warning: Could not cluster {family}: {e}")
        
        family_blocks[family] = (len(row_order), len(row_order) + len(family_genes))
        row_order.extend(family_genes)
    
    heatmap_data = heatmap_data.loc[row_order]
    gene_families = gene_families.loc[row_order]
    
    print(f"  Row order determined ({len(row_order)} genes in {len(unique_families)} families)")
    
    # Step 8: Save heatmap matrix
    print("\nStep 8: Saving heatmap matrix...")
    
    matrix_file = output_dir / "cyp_heatmap_matrix.tsv"
    heatmap_data.to_csv(matrix_file, sep='\t')
    print(f"  Saved: {matrix_file}")
    
    # Step 9: Create heatmap plot (Figure 6A style)
    print("\nStep 9: Creating heatmap plot (publication style)...")
    
    n_genes = len(heatmap_data)
    n_samples = len(heatmap_data.columns)
    
    # Sort samples within each condition group (Leaf first, then Root)
    if metadata_df is not None and contrast_factor in metadata_df.columns:
        condition_col = metadata_df[contrast_factor]
        samples_B = sorted([s for s in heatmap_data.columns if s in condition_col.index and condition_col[s] == contrast_B])
        samples_A = sorted([s for s in heatmap_data.columns if s in condition_col.index and condition_col[s] == contrast_A])
        column_order = samples_B + samples_A
        heatmap_data = heatmap_data[column_order]
        display_labels = column_order  # Keep original sample names
    else:
        display_labels = list(heatmap_data.columns)
    
    # Figure dimensions - match publication style
    fig_width = max(10, 2.5 + n_samples * 0.5)
    fig_height = max(8, 1 + n_genes * 0.35)
    
    # Create figure with GridSpec for precise layout
    # Layout: [family brackets | gene IDs (locus) | heatmap | colorbar]
    # Gene IDs on LEFT side, family labels grouped above them
    from matplotlib.gridspec import GridSpec
    
    fig = plt.figure(figsize=(fig_width, fig_height))
    
    # Define grid: family bracket space (8%), gene labels (18%), heatmap (66%), colorbar (8%)
    gs = GridSpec(1, 4, figure=fig, width_ratios=[0.08, 0.18, 0.66, 0.08], wspace=0.02)
    ax_family = fig.add_subplot(gs[0, 0])
    ax_genes = fig.add_subplot(gs[0, 1])
    ax_heatmap = fig.add_subplot(gs[0, 2])
    ax_cbar = fig.add_subplot(gs[0, 3])
    
    # Perform hierarchical clustering on rows if requested (before plotting)
    if row_cluster and n_genes > 1:
        from scipy.cluster.hierarchy import dendrogram, linkage
        try:
            row_linkage = linkage(heatmap_data.values, method='ward', metric='euclidean')
            # Get the order from clustering without plotting dendrogram
            from scipy.cluster.hierarchy import leaves_list
            row_order_idx = leaves_list(row_linkage)
            heatmap_data = heatmap_data.iloc[row_order_idx]
            gene_families = gene_families.iloc[row_order_idx]
            print("  Rows clustered hierarchically")
        except Exception as e:
            print(f"  Warning: Could not cluster rows: {e}")
    
    # Draw the heatmap
    vmin = -3 if scale_method == 'zscore' else None
    vmax = 3 if scale_method == 'zscore' else None
    
    im = ax_heatmap.imshow(heatmap_data.values, aspect='auto', cmap='RdBu_r', 
                           vmin=vmin, vmax=vmax, interpolation='nearest')
    
    # X-axis labels (samples) at bottom
    ax_heatmap.set_xticks(range(n_samples))
    ax_heatmap.set_xticklabels(display_labels, rotation=45, ha='right', fontsize=10, fontweight='bold')
    ax_heatmap.xaxis.set_ticks_position('bottom')
    
    # No Y-axis labels on heatmap (they go in the gene ID axis on the left)
    ax_heatmap.set_yticks([])
    ax_heatmap.set_ylabel('')
    
    # Set heatmap limits
    ax_heatmap.set_xlim(-0.5, n_samples - 0.5)
    ax_heatmap.set_ylim(n_genes - 0.5, -0.5)
    
    # === Gene IDs (locus) on the LEFT side ===
    ax_genes.set_xlim(0, 1)
    ax_genes.set_ylim(n_genes - 0.5, -0.5)
    ax_genes.set_xticks([])
    ax_genes.set_yticks([])
    ax_genes.spines['top'].set_visible(False)
    ax_genes.spines['right'].set_visible(False)
    ax_genes.spines['bottom'].set_visible(False)
    ax_genes.spines['left'].set_visible(False)
    
    # Draw gene IDs aligned to the right (next to heatmap)
    for i, gene_id in enumerate(heatmap_data.index):
        ax_genes.text(0.95, i, gene_id, fontsize=7, va='center', ha='right', 
                      fontfamily='monospace')
    
    # === Family labels with brackets on the FAR LEFT ===
    ax_family.set_xlim(0, 1)
    ax_family.set_ylim(n_genes - 0.5, -0.5)
    ax_family.set_xticks([])
    ax_family.set_yticks([])
    ax_family.spines['top'].set_visible(False)
    ax_family.spines['right'].set_visible(False)
    ax_family.spines['bottom'].set_visible(False)
    ax_family.spines['left'].set_visible(False)
    
    # Find family blocks and add bracket annotations
    current_family = None
    block_start = 0
    family_blocks_plot = []
    
    for i, (gene, family) in enumerate(gene_families.items()):
        if family != current_family:
            if current_family is not None:
                family_blocks_plot.append((current_family, block_start, i - 1))
            current_family = family
            block_start = i
    # Add last block
    if current_family is not None:
        family_blocks_plot.append((current_family, block_start, len(gene_families) - 1))
    
    # Draw family brackets and labels
    for family, start, end in family_blocks_plot:
        mid = (start + end) / 2
        height = end - start + 1
        
        # Draw bracket: vertical line + horizontal ticks (pointing right toward genes)
        bracket_x = 0.85
        ax_family.plot([bracket_x, bracket_x], [start - 0.3, end + 0.3], 
                       color='black', linewidth=1.5, clip_on=False)
        ax_family.plot([bracket_x, bracket_x + 0.12], [start - 0.3, start - 0.3], 
                       color='black', linewidth=1.5, clip_on=False)
        ax_family.plot([bracket_x, bracket_x + 0.12], [end + 0.3, end + 0.3], 
                       color='black', linewidth=1.5, clip_on=False)
        
        # Family label (rotated if tall block, horizontal if short)
        ax_family.text(0.4, mid, family, fontsize=9, fontweight='bold',
                       ha='center', va='center', rotation=90 if height > 3 else 0)
    
    # Colorbar
    cbar = fig.colorbar(im, cax=ax_cbar)
    cbar.set_label(scale_label, fontsize=10)
    cbar.ax.tick_params(labelsize=8)
    
    # Title
    fig.suptitle(f'CYP Gene Expression Heatmap ({n_genes} genes)', 
                 fontsize=14, fontweight='bold', y=0.98)
    
    # Save
    pdf_file = output_dir / "cyp_heatmap.pdf"
    plt.savefig(pdf_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {pdf_file}")
    
    png_file = output_dir / "cyp_heatmap.png"
    plt.savefig(png_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {png_file}")
    
    plt.close()
    
    print("\nCYP heatmap generation complete!")
    print(f"  Total CYP DEGs plotted: {n_genes}")
    print(f"  CYP families: {', '.join(unique_families)}")


def generate_deg_heatmap(dds, results_df, output_dir, count_matrix_for_plots,
                         padj_cutoff=0.05, lfc_cutoff=2.0, 
                         top_n=50, scale_method='center',
                         contrast_A='root', contrast_B='leaf', contrast_factor='condition'):
    """
    Generate a general DEG heatmap showing ALL differentially expressed genes.
    
    This heatmap shows statistically significant DEGs (not just CYP genes)
    with log2FC annotations to clearly indicate effect sizes.
    
    Statistical criteria:
    - padj < padj_cutoff (adjusted p-value, FDR corrected)
    - |log2FoldChange| > lfc_cutoff (effect size)
    
    Parameters:
        top_n: Maximum number of DEGs to show (sorted by |log2FC|)
    """
    print("\n" + "=" * 60)
    print("Generating DEG Heatmap (Statistically Significant Genes)")
    print("=" * 60)
    
    # Step 1: Filter for statistically significant DEGs
    print("\nStep 1: Applying statistical filters...")
    
    deg_df = results_df.copy()
    n_total = len(deg_df)
    
    # Remove genes with missing padj values #new
    # Step 1: Remove NA (QUALITY CONTROL - you have to do this)
    deg_df = deg_df[deg_df['padj'].notna()]
    
    # Filter by adjusted p-value (statistical significance)
#   # OLD CODE: deg_df = deg_df[deg_df['padj'].notna()]
   # deg_df = deg_df[deg_df['padj'] < padj_cutoff]
   # n_padj = len(deg_df)
   # print(f"  ✓ Genes with padj < {padj_cutoff}: {n_padj} / {n_total}")
   ## Step 2: NO statistical filtering (set padj_cutoff to 1.0 to disable statistical filtering)
    if padj_cutoff < 1.0:
        deg_df = deg_df[deg_df['padj'] < padj_cutoff]
        n_padj = len(deg_df)
        print(f"  ✓ Genes with padj < {padj_cutoff}: {n_padj} / {n_total}")
    else:
        print(f"  ✓ padj filter DISABLED (cutoff={padj_cutoff}): keeping all {len(deg_df)} genes")
    
    # Filter by log2 fold change (biological significance)
    if lfc_cutoff > 0:
        deg_df = deg_df[np.abs(deg_df['log2FoldChange']) > lfc_cutoff]
        n_lfc = len(deg_df)
        print(f"  ✓ Genes with |log2FC| > {lfc_cutoff}: {n_lfc}")
    else:
        print(f"  ✓ log2FC filter DISABLED (cutoff=0): keeping all {len(deg_df)} genes")
    
    if len(deg_df) == 0:
        print("  ⚠ WARNING: No genes passed DEG filters. Cannot create DEG heatmap.")
        print("  Consider lowering --padj or --lfc thresholds.")
        return
    
    # Step 2: Sort by absolute log2FC and take top N
    deg_df['abs_log2FC'] = np.abs(deg_df['log2FoldChange'])
    deg_df = deg_df.sort_values('abs_log2FC', ascending=False)
    
    if len(deg_df) > top_n:
        print(f"\nStep 2: Selecting top {top_n} DEGs by |log2FC|...")
        deg_df = deg_df.head(top_n)
    else:
        print(f"\nStep 2: Using all {len(deg_df)} DEGs...")
    
    # Save DEG list with statistics
    deg_stats_file = output_dir / "deg_heatmap_genes.tsv"
    deg_export = deg_df[['baseMean', 'log2FoldChange', 'pvalue', 'padj']].copy()
    deg_export['direction'] = np.where(deg_df['log2FoldChange'] > 0, 
                                        f'UP_in_{contrast_A}', f'UP_in_{contrast_B}')
    deg_export = deg_export.sort_values('log2FoldChange', ascending=False)
    deg_export.to_csv(deg_stats_file, sep='\t')
    print(f"  Saved DEG statistics: {deg_stats_file}")
    
    # Step 3: Get normalized expression data
    print("\nStep 3: Getting normalized expression data...")
    
    normalized_counts = get_normalized_counts(dds, count_matrix_for_plots)
    deg_gene_ids = deg_df.index.tolist()
    
    # Filter to DEGs that exist in count matrix
    available_genes = [g for g in deg_gene_ids if g in normalized_counts.index]
    if len(available_genes) < len(deg_gene_ids):
        print(f"  Note: {len(deg_gene_ids) - len(available_genes)} genes not in count matrix")
    
    heatmap_data = normalized_counts.loc[available_genes].copy()
    print(f"  Expression matrix: {heatmap_data.shape[0]} genes × {heatmap_data.shape[1]} samples")
    
    # Step 4: Transform data (log2 + center/zscore)
    print(f"\nStep 4: Transforming data (log2 + {scale_method})...")
    
    heatmap_data = np.log2(heatmap_data + 1)
    
    if scale_method == 'center':
        row_means = heatmap_data.mean(axis=1)
        heatmap_data = heatmap_data.sub(row_means, axis=0)
        scale_label = "log2(normalized + 1)\ncentered"
    elif scale_method == 'zscore':
        row_means = heatmap_data.mean(axis=1)
        row_stds = heatmap_data.std(axis=1)
        row_stds = row_stds.replace(0, 1)
        heatmap_data = heatmap_data.sub(row_means, axis=0).div(row_stds, axis=0)
        scale_label = "z-score"
    
    # Step 5: Get log2FC values for annotation
    log2fc_values = deg_df.loc[available_genes, 'log2FoldChange']
    padj_values = deg_df.loc[available_genes, 'padj']
    
    # Step 6: Order samples by condition
    print("\nStep 5: Ordering samples...")
    
    metadata_df = dds.obs if hasattr(dds, 'obs') else None
    if metadata_df is not None and contrast_factor in metadata_df.columns:
        condition_col = metadata_df[contrast_factor]
        samples_B = sorted([s for s in heatmap_data.columns if s in condition_col.index and condition_col[s] == contrast_B])
        samples_A = sorted([s for s in heatmap_data.columns if s in condition_col.index and condition_col[s] == contrast_A])
        column_order = samples_B + samples_A
        heatmap_data = heatmap_data[column_order]
    else:
        column_order = list(heatmap_data.columns)
    
    # Sort genes by log2FC (up-regulated at top, down-regulated at bottom)
    gene_order = log2fc_values.sort_values(ascending=False).index.tolist()
    heatmap_data = heatmap_data.loc[gene_order]
    log2fc_values = log2fc_values.loc[gene_order]
    padj_values = padj_values.loc[gene_order]
    
    # Step 7: Create heatmap plot
    print("\nStep 6: Creating DEG heatmap...")
    
    n_genes = len(heatmap_data)
    n_samples = len(heatmap_data.columns)
    
    fig_width = max(12, 4 + n_samples * 0.6)
    fig_height = max(10, 2 + n_genes * 0.3)
    
    fig, axes = plt.subplots(1, 3, figsize=(fig_width, fig_height),
                              gridspec_kw={'width_ratios': [0.08, 0.82, 0.10], 'wspace': 0.02})
    ax_lfc, ax_heatmap, ax_cbar = axes
    
    # Draw the heatmap
    vmin = -3 if scale_method == 'zscore' else None
    vmax = 3 if scale_method == 'zscore' else None
    
    im = ax_heatmap.imshow(heatmap_data.values, aspect='auto', cmap='RdBu_r',
                            vmin=vmin, vmax=vmax, interpolation='nearest')
    
    # X-axis: sample labels
    ax_heatmap.set_xticks(range(n_samples))
    ax_heatmap.set_xticklabels(column_order, rotation=45, ha='right', fontsize=9)
    
    # Y-axis: gene labels with log2FC and significance
    gene_labels = []
    for gene in heatmap_data.index:
        lfc = log2fc_values[gene]
        padj = padj_values[gene]
        direction = "↑" if lfc > 0 else "↓"
        sig_stars = "***" if padj < 0.001 else "**" if padj < 0.01 else "*"
        gene_labels.append(f"{gene}")
    
    ax_heatmap.set_yticks(range(n_genes))
    ax_heatmap.set_yticklabels(gene_labels, fontsize=7, fontfamily='monospace')
    
    # Add condition labels at top
    if metadata_df is not None and contrast_factor in metadata_df.columns:
        n_B = len(samples_B)
        n_A = len(samples_A)
        ax_heatmap.text(n_B/2 - 0.5, -1.5, contrast_B, ha='center', fontsize=11, fontweight='bold')
        ax_heatmap.text(n_B + n_A/2 - 0.5, -1.5, contrast_A, ha='center', fontsize=11, fontweight='bold')
        # Add separator line
        ax_heatmap.axvline(x=n_B - 0.5, color='black', linewidth=2)
    
    # Log2FC annotation bar on the left
    ax_lfc.set_xlim(0, 1)
    ax_lfc.set_ylim(n_genes - 0.5, -0.5)
    ax_lfc.set_xticks([])
    ax_lfc.set_yticks([])
    ax_lfc.set_title('log2FC', fontsize=9, fontweight='bold')
    
    for i, gene in enumerate(heatmap_data.index):
        lfc = log2fc_values[gene]
        padj = padj_values[gene]
        color = '#d62728' if lfc > 0 else '#1f77b4'  # Red=up, Blue=down
        sig = "***" if padj < 0.001 else "**" if padj < 0.01 else "*"
        ax_lfc.text(0.5, i, f"{lfc:+.1f}{sig}", fontsize=6, ha='center', va='center',
                    color=color, fontweight='bold', fontfamily='monospace')
    
    ax_lfc.spines['top'].set_visible(False)
    ax_lfc.spines['right'].set_visible(False)
    ax_lfc.spines['bottom'].set_visible(False)
    ax_lfc.spines['left'].set_visible(False)
    
    # Colorbar
    cbar = fig.colorbar(im, cax=ax_cbar)
    cbar.set_label(scale_label, fontsize=10)
    
    # Title with statistical info
    title = f"DEG Heatmap: {n_genes} Differentially Expressed Genes\n"
    title += f"(padj < {padj_cutoff}, |log2FC| > {lfc_cutoff})"
    fig.suptitle(title, fontsize=12, fontweight='bold', y=0.98)
    
    # Add legend for significance
    legend_text = "Significance: * p<0.05, ** p<0.01, *** p<0.001"
    fig.text(0.5, 0.01, legend_text, ha='center', fontsize=8, style='italic')
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    # Save
    pdf_file = output_dir / "deg_heatmap.pdf"
    plt.savefig(pdf_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {pdf_file}")
    
    png_file = output_dir / "deg_heatmap.png"
    plt.savefig(png_file, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {png_file}")
    
    plt.close()
    
    # Summary statistics
    n_up = sum(log2fc_values > 0)
    n_down = sum(log2fc_values < 0)
    
    print(f"\n  DEG Heatmap Summary:")
    print(f"  ├── Total DEGs shown: {n_genes}")
    print(f"  ├── Upregulated in {contrast_A}: {n_up} (log2FC > 0)")
    print(f"  ├── Upregulated in {contrast_B}: {n_down} (log2FC < 0)")
    print(f"  ├── Mean |log2FC|: {np.abs(log2fc_values).mean():.2f}")
    print(f"  └── Median padj: {padj_values.median():.2e}")


def generate_summary(count_file, metadata_file, results_df, dds, output_dir, args):
    """Generate summary report."""
    print("Generating summary report...")
    
    summary_file = output_dir / "analysis_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("PyDESeq2 Differential Expression Analysis Summary\n")
        f.write("=" * 60 + "\n\n")
        f.write("Input files:\n")
        f.write(f"  Count matrix: {count_file}\n")
        f.write(f"  Metadata: {metadata_file}\n")
        if hasattr(args, 'cyp_family_map') and args.cyp_family_map:
            f.write(f"  CYP family map: {args.cyp_family_map}\n")
        f.write("\n")
        
        f.write("Analysis parameters:\n")
        f.write(f"  Design formula: {args.design or 'auto-detected'}\n")
        f.write(f"  Contrast: {args.contrast_A} vs {args.contrast_B} (factor: {args.contrast_factor})\n")
        f.write(f"  padj cutoff: {args.padj}\n")
        f.write(f"  log2FC cutoff: {args.lfc}\n")
        f.write(f"  Root-up only: {args.root_up_only}\n")
        f.write(f"  Scale method: {args.scale}\n")
        f.write("\n")
        
        f.write("Data summary:\n")
        if hasattr(dds, 'X'):
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
        
        cyp_deg_file = output_dir / "cyp_deg_filtered.tsv"
        if cyp_deg_file.exists():
            cyp_deg = pd.read_csv(cyp_deg_file, sep='\t', index_col=0)
            f.write("\nCYP gene analysis:\n")
            f.write(f"  CYP DEGs (filtered): {len(cyp_deg)}\n")
            if 'cyp_family' in cyp_deg.columns:
                family_counts = cyp_deg['cyp_family'].value_counts()
                f.write("  CYP families:\n")
                for fam, count in family_counts.items():
                    f.write(f"    - {fam}: {count}\n")
        
        f.write("\nOutput files:\n")
        f.write("  PyDESeq2 results: pydeseq2_results.tsv\n")
        f.write("  Filtered DEGs: deg_filtered.tsv\n")
        f.write("  CYP DEGs: cyp_deg_filtered.tsv\n")
        f.write("  CYP heatmap matrix: cyp_heatmap_matrix.tsv\n")
        f.write("  CYP heatmap: cyp_heatmap.pdf, cyp_heatmap.png\n")
        f.write("  Quality control plots: qc_*.pdf\n")
        f.write("  Analysis plots: pydeseq2_*.pdf\n")
    
    print(f"Summary saved: {summary_file}")


def main():
    parser = argparse.ArgumentParser(
        description="PyDESeq2 differential expression analysis with CYP heatmap generation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage for Figure 6A-style CYP heatmap:

  python pydeseq2_analysis.py count_matrix.tsv metadata.tsv \\
      -o results \\
      --contrast-factor condition \\
      --contrast-A root \\
      --contrast-B leaf \\
      --cyp-family-map cyp_families.tsv \\
      --root-up-only \\
      --lfc 2.0 \\
      --scale center

CYP family map format (TSV):
  gene_id    cyp_family    cyp_clan
  LOC1234    CYP71         CYP71_clan
  LOC5678    CYP72         CYP72_clan
        """
    )
    
    parser.add_argument("count_matrix", help="Path to gene count matrix TSV file")
    parser.add_argument("metadata", help="Path to sample metadata TSV file")
    
    parser.add_argument("-o", "--output", default="pydeseq2_results",
                       help="Output directory (default: pydeseq2_results)")
    
    parser.add_argument("--design", default=None,
                       help="Design formula (e.g., 'treatment' or 'group+condition'). "
                            "If not specified, will auto-detect from metadata")
    
    parser.add_argument("--contrast-factor", default="condition",
                       help="Metadata column for contrast (default: condition)")
    parser.add_argument("--contrast-A", default="root",
                       help="Numerator condition - positive log2FC (default: root)")
    parser.add_argument("--contrast-B", default="leaf",
                       help="Denominator condition (default: leaf)")
    
    parser.add_argument("--padj", type=float, default=0.05,
                       help="Adjusted p-value cutoff (default: 0.05, set 1.0 to disable)")
    parser.add_argument("--lfc", type=float, default=2.0,
                       help="Absolute log2 fold change cutoff (default: 2.0, set 0 to disable)")
    parser.add_argument("--root-up-only", action="store_true", default=False,
                       help="Keep only genes upregulated in contrast-A (default: False)")
    parser.add_argument("--no-root-up-only", action="store_true", default=False,
                       help="Include both up and down-regulated genes")
    parser.add_argument("--top-deg", type=int, default=50,
                       help="Number of top DEGs to show in DEG heatmap (default: 50)")
    
    parser.add_argument("--cyp-family-map", default=None,
                       help="TSV file mapping gene_id to cyp_family (preferred)")
    parser.add_argument("--cyp-genes", default=None,
                       help="Text file with CYP gene IDs (one per line)")
    
    parser.add_argument("--sample-order", default=None,
                       help="Text file with sample names in desired order")
    
    parser.add_argument("--scale", choices=['center', 'zscore'], default='center',
                       help="Scaling method for heatmap (default: center)")
    parser.add_argument("--row-cluster", action="store_true", default=True,
                       help="Cluster rows (genes) hierarchically (default: True)")
    parser.add_argument("--no-row-cluster", action="store_true", default=False,
                       help="Disable row clustering")
    parser.add_argument("--col-cluster", action="store_true", default=False,
                       help="Cluster columns (samples) - NOT recommended for Figure 6A")
    parser.add_argument("--no-col-cluster", action="store_true", default=True,
                       help="Do not cluster columns (default)")
    
    args = parser.parse_args()
    
    if args.no_root_up_only:
        args.root_up_only = False
    if args.no_row_cluster:
        args.row_cluster = False
    
    if not os.path.exists(args.count_matrix):
        print(f"ERROR: Count matrix file not found: {args.count_matrix}")
        return 1
    
    if not os.path.exists(args.metadata):
        print(f"ERROR: Metadata file not found: {args.metadata}")
        return 1
    
    cyp_family_map = None
    cyp_genes = None
    
    if args.cyp_family_map:
        if not os.path.exists(args.cyp_family_map):
            print(f"ERROR: CYP family map file not found: {args.cyp_family_map}")
            return 1
        cyp_family_map = load_cyp_family_map(args.cyp_family_map)
    
    if args.cyp_genes:
        if not os.path.exists(args.cyp_genes):
            print(f"ERROR: CYP genes file not found: {args.cyp_genes}")
            return 1
        cyp_genes = load_cyp_gene_list(args.cyp_genes)
    
    sample_order = None
    if args.sample_order:
        if not os.path.exists(args.sample_order):
            print(f"ERROR: Sample order file not found: {args.sample_order}")
            return 1
        sample_order = load_sample_order(args.sample_order)
    
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("PyDESeq2 Differential Expression Analysis")
    print("=" * 60)
    print(f"Count matrix: {args.count_matrix}")
    print(f"Metadata: {args.metadata}")
    print(f"Output directory: {output_dir}")
    print(f"Contrast: {args.contrast_A} vs {args.contrast_B} (factor: {args.contrast_factor})")
    print(f"Filters: padj < {args.padj}, |log2FC| > {args.lfc}")
    print(f"Root-up only: {args.root_up_only}")
    print(f"Scale method: {args.scale}")
    print()
    
    try:
        count_matrix, metadata = read_data(args.count_matrix, args.metadata)
        
        if args.design:
            design_formula = args.design
        else:
            if args.contrast_factor in metadata.columns:
                design_formula = args.contrast_factor
            elif 'treatment' in metadata.columns:
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
                    raise ValueError("No suitable grouping variable found in metadata.")
        
        print(f"Using design formula: {design_formula}\n")
        
        generate_qc_plots(count_matrix, metadata, output_dir)
        
        dds, stat_res, results_df, count_matrix_for_plots = run_pydeseq2_analysis(
            count_matrix, metadata, design_formula, output_dir,
            contrast_factor=args.contrast_factor,
            contrast_A=args.contrast_A,
            contrast_B=args.contrast_B
        )
        
        generate_ma_volcano_plots(results_df, output_dir)
        
        # Generate general DEG heatmap (all significant DEGs, not just CYP)
        generate_deg_heatmap(
            dds=dds,
            results_df=results_df,
            output_dir=output_dir,
            count_matrix_for_plots=count_matrix_for_plots,
            padj_cutoff=args.padj,
            lfc_cutoff=args.lfc,
            top_n=args.top_deg,
            scale_method=args.scale,
            contrast_A=args.contrast_A,
            contrast_B=args.contrast_B,
            contrast_factor=args.contrast_factor
        )
        
        if cyp_family_map is not None or cyp_genes is not None:
            generate_cyp_heatmap(
                dds=dds,
                results_df=results_df,
                output_dir=output_dir,
                count_matrix_for_plots=count_matrix_for_plots,
                cyp_family_map=cyp_family_map,
                cyp_genes=cyp_genes,
                padj_cutoff=args.padj,
                lfc_cutoff=args.lfc,
                root_up_only=args.root_up_only,
                sample_order=sample_order,
                scale_method=args.scale,
                row_cluster=args.row_cluster,
                col_cluster=args.col_cluster,
                contrast_A=args.contrast_A,
                contrast_B=args.contrast_B,
                contrast_factor=args.contrast_factor
            )
        else:
            print("\nNOTE: No CYP gene mapping provided. Skipping CYP heatmap generation.")
            print("      To generate Figure 6A-style heatmap, provide --cyp-family-map or --cyp-genes")
        
        generate_summary(args.count_matrix, args.metadata, results_df, dds, output_dir, args)
        
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
