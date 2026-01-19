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
    
    deg_df = deg_df[deg_df['padj'].notna() & (deg_df['padj'] < padj_cutoff)]
    n_padj = len(deg_df)
    print(f"  Genes with padj < {padj_cutoff}: {n_padj} (from {n_total})")
    
    if lfc_cutoff > 0:
        deg_df = deg_df[np.abs(deg_df['log2FoldChange']) > lfc_cutoff]
        n_lfc = len(deg_df)
        print(f"  Genes with |log2FC| > {lfc_cutoff}: {n_lfc}")
    
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
    
    # Step 9: Create heatmap plot
    print("\nStep 9: Creating heatmap plot...")
    
    n_genes = len(heatmap_data)
    n_samples = len(heatmap_data.columns)
    fig_width = max(10, 2 + n_samples * 0.6)
    fig_height = max(10, 2 + n_genes * 0.25)
    
    n_families = len(unique_families)
    if n_families <= 10:
        family_palette = sns.color_palette("tab10", n_families)
    elif n_families <= 20:
        family_palette = sns.color_palette("tab20", n_families)
    else:
        family_palette = sns.color_palette("husl", n_families)
    family_color_dict = dict(zip(unique_families, family_palette))
    
    condition_colors = {contrast_B: '#4393C3', contrast_A: '#D6604D'}
    
    row_colors = gene_families.map(family_color_dict)
    row_colors.name = 'CYP Family'
    
    col_colors = None
    if metadata_df is not None and contrast_factor in metadata_df.columns:
        col_conditions = metadata_df.loc[heatmap_data.columns, contrast_factor]
        col_colors = col_conditions.map(condition_colors)
        col_colors.name = contrast_factor.capitalize()
    
    try:
        g = sns.clustermap(
            heatmap_data,
            cmap='RdBu_r',
            center=0,
            vmin=-3 if scale_method == 'zscore' else None,
            vmax=3 if scale_method == 'zscore' else None,
            row_cluster=row_cluster,
            col_cluster=col_cluster,
            row_colors=row_colors,
            col_colors=col_colors,
            figsize=(fig_width, fig_height),
            dendrogram_ratio=(0.1, 0.05) if row_cluster else (0.0, 0.0),
            cbar_pos=(0.02, 0.8, 0.02, 0.15),
            cbar_kws={'label': scale_label, 'shrink': 0.5},
            yticklabels=True,
            xticklabels=True,
            linewidths=0.3,
            linecolor='white',
            method='ward' if row_cluster else 'single',
            metric='euclidean'
        )
        
        g.fig.suptitle(
            f'CYP Gene Expression Heatmap ({n_genes} genes)\n'
            f'{contrast_A.capitalize()} vs {contrast_B.capitalize()} (padj < {padj_cutoff}, |log2FC| > {lfc_cutoff})',
            fontsize=14, fontweight='bold', y=1.02
        )
        
        g.ax_heatmap.set_yticklabels(
            g.ax_heatmap.get_yticklabels(),
            rotation=0, fontsize=7, ha='right'
        )
        g.ax_heatmap.set_ylabel('Gene ID', fontsize=10)
        
        g.ax_heatmap.set_xticklabels(
            g.ax_heatmap.get_xticklabels(),
            rotation=45, fontsize=9, ha='right'
        )
        g.ax_heatmap.set_xlabel('Sample', fontsize=10)
        
        if hasattr(g, 'ax_row_dendrogram') and g.ax_row_dendrogram is not None:
            g.ax_row_dendrogram.set_visible(row_cluster)
        if hasattr(g, 'ax_col_dendrogram') and g.ax_col_dendrogram is not None:
            g.ax_col_dendrogram.set_visible(False)
        
        legend_elements = []
        
        for family in unique_families:
            n_in_family = (gene_families == family).sum()
            legend_elements.append(Patch(
                facecolor=family_color_dict[family], 
                edgecolor='black', 
                linewidth=0.5,
                label=f'{family} ({n_in_family})'
            ))
        
        legend_elements.append(Patch(facecolor='white', edgecolor='white', label=''))
        legend_elements.append(Patch(
            facecolor=condition_colors[contrast_B],
            edgecolor='black', 
            linewidth=0.5,
            label=f'{contrast_B.capitalize()}'
        ))
        legend_elements.append(Patch(
            facecolor=condition_colors[contrast_A],
            edgecolor='black',
            linewidth=0.5,
            label=f'{contrast_A.capitalize()}'
        ))
        
        g.fig.legend(
            handles=legend_elements,
            title='Legend',
            loc='center left', 
            bbox_to_anchor=(1.01, 0.5),
            fontsize=8,
            title_fontsize=10,
            frameon=True,
            edgecolor='black'
        )
        
        pdf_file = output_dir / "cyp_heatmap.pdf"
        plt.savefig(pdf_file, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"  Saved: {pdf_file}")
        
        png_file = output_dir / "cyp_heatmap.png"
        plt.savefig(png_file, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"  Saved: {png_file}")
        
        plt.close()
            
    except Exception as e:
        print(f"  ERROR creating clustermap: {e}")
        import traceback
        traceback.print_exc()
        
        print("  Attempting fallback heatmap...")
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        sns.heatmap(
            heatmap_data,
            cmap='RdBu_r',
            center=0,
            yticklabels=True,
            xticklabels=True,
            ax=ax,
            cbar_kws={'label': scale_label}
        )
        ax.set_title(f'CYP Gene Expression Heatmap ({n_genes} genes)')
        ax.set_ylabel('Gene ID')
        ax.set_xlabel('Sample')
        plt.tight_layout()
        
        pdf_file = output_dir / "cyp_heatmap.pdf"
        plt.savefig(pdf_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved fallback heatmap: {pdf_file}")
    
    print("\nCYP heatmap generation complete!")
    print(f"  Total CYP DEGs plotted: {n_genes}")
    print(f"  CYP families: {', '.join(unique_families)}")


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
                       help="Adjusted p-value cutoff (default: 0.05)")
    parser.add_argument("--lfc", type=float, default=2.0,
                       help="Absolute log2 fold change cutoff (default: 2.0, set 0 to disable)")
    parser.add_argument("--root-up-only", action="store_true", default=False,
                       help="Keep only genes upregulated in contrast-A (default: False)")
    parser.add_argument("--no-root-up-only", action="store_true", default=False,
                       help="Include both up and down-regulated genes")
    
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
