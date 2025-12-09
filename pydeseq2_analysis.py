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
import re
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
        # Use the stored count matrix - it's already a DataFrame with correct structure
        # count_matrix_for_plots is (samples × genes) as a pandas DataFrame
        if count_matrix_for_plots is None or not isinstance(count_matrix_for_plots, pd.DataFrame):
            raise ValueError("Count matrix for plots is not available or not a DataFrame")
        
        count_data = count_matrix_for_plots.copy()
        print(f"  Count matrix shape: {count_data.shape} (samples × genes)")
        
        # Normalize by size factors if available
        # count_data is (samples × genes) DataFrame
        if hasattr(dds, 'size_factors') and dds.size_factors is not None:
            # size_factors is a Series/array with one value per sample (row)
            # Divide each row (sample) by its size factor
            normalized_counts = count_data.div(dds.size_factors, axis=0)
            print("  Normalized counts using size factors")
        elif hasattr(dds, 'obs') and 'size_factors' in dds.obs.columns:
            # Size factors stored in obs
            normalized_counts = count_data.div(dds.obs['size_factors'], axis=0)
            print("  Normalized counts using size factors from obs")
        else:
            # If no size factors, use raw counts
            print("  Warning: No size factors found, using raw counts")
            normalized_counts = count_data
        
        # Calculate variance across samples (axis=0) to find variable genes
        # normalized_counts is (samples × genes), so var(axis=0) calculates variance per gene (column)
        print("  Calculating gene variance...")
        gene_vars = normalized_counts.var(axis=0)
        
        # Get top 50 most variable genes
        top_var_genes = gene_vars.nlargest(50).index.tolist()
        print(f"  Selected top {len(top_var_genes)} variable genes")
        
        # Select only the top variable genes (columns)
        # normalized_counts is (samples × genes), so we select gene columns
        heatmap_data = normalized_counts[top_var_genes]
        print(f"  Heatmap data shape before transpose: {heatmap_data.shape}")
        
        # Transpose for visualization: genes as rows, samples as columns
        # This makes it easier to see genes (rows) across samples (columns)
        heatmap_data = heatmap_data.T
        print(f"  Heatmap data shape after transpose: {heatmap_data.shape} (genes × samples)")
        
        # Center each gene by subtracting its mean across samples
        # After transpose: genes are rows, samples are columns
        # Subtract mean across columns (axis=1) for each gene (row)
        heatmap_data = heatmap_data.subtract(heatmap_data.mean(axis=1), axis=0)
        print("  Centered data (subtracted gene means)")
        
        # Create annotation for samples (columns)
        # In AnnData/PyDESeq2, metadata is stored in dds.obs, not dds.metadata
        print("  Preparing sample annotations...")
        metadata_df = None
        if hasattr(dds, 'obs'):
            metadata_df = dds.obs
            print(f"  Found metadata in dds.obs with shape: {metadata_df.shape}")
        elif hasattr(dds, 'metadata'):
            metadata_df = dds.metadata
            print(f"  Found metadata in dds.metadata with shape: {metadata_df.shape}")
        
        annotation_col = None
        if metadata_df is not None and len(metadata_df) > 0:
            # heatmap_data.columns are sample names (after transpose)
            # Ensure metadata index matches heatmap column names
            try:
                # Align metadata with heatmap columns (sample names)
                metadata_aligned = metadata_df.reindex(heatmap_data.columns)
                
                # Find annotation column
                if 'treatment' in metadata_aligned.columns:
                    annotation_col = metadata_aligned[['treatment']]
                    print("  Using 'treatment' for annotations")
                elif 'group' in metadata_aligned.columns:
                    annotation_col = metadata_aligned[['group']]
                    print("  Using 'group' for annotations")
                elif 'condition' in metadata_aligned.columns:
                    annotation_col = metadata_aligned[['condition']]
                    print("  Using 'condition' for annotations")
                else:
                    print("  No suitable annotation column found (treatment/group/condition)")
            except Exception as e:
                print(f"  Warning: Could not align metadata: {e}")
                annotation_col = None
        else:
            print("  No metadata available for annotations")
        
        # Create the heatmap
        print("  Creating heatmap plot...")
        
        # Intelligently identify gene families from actual gene names in the data
        print("  Analyzing gene names to identify families...")
        gene_families = {}
        family_patterns = {}
        
        # Get all gene names
        gene_names = [str(gene) for gene in heatmap_data.index]
        
        # Strategy 1: Extract common prefixes (before separators like _, -, |, .)
        # This works for patterns like: TRINITY_DN123, LOC108196229, AT1G01010, etc.
        prefix_patterns = {}
        for gene in gene_names:
            # Try different separators
            for sep in ['_', '-', '|', '.', ':']:
                if sep in gene:
                    parts = gene.split(sep)
                    if len(parts) > 1:
                        prefix = parts[0]
                        # Only use if prefix is meaningful (3+ chars, not just numbers)
                        if len(prefix) >= 3 and not prefix.isdigit():
                            if prefix not in prefix_patterns:
                                prefix_patterns[prefix] = []
                            prefix_patterns[prefix].append(gene)
                            break
        
        # Strategy 2: Extract patterns with numbers and letters (e.g., CYP71D, AT1G, LOC108)
        pattern_families = {}
        for gene in gene_names:
            # Pattern: letters followed by numbers (e.g., CYP71, LOC108, AT1G)
            match = re.search(r'^([A-Za-z]+\d+[A-Za-z]*)', gene)
            if match:
                pattern = match.group(1)
                if len(pattern) >= 3:  # Meaningful pattern
                    if pattern not in pattern_families:
                        pattern_families[pattern] = []
                    pattern_families[pattern].append(gene)
        
        # Strategy 3: Extract first N characters if they're consistent
        # This catches patterns where families share the same prefix length
        char_prefixes = {}
        for gene in gene_names:
            # Try different prefix lengths
            for prefix_len in [4, 5, 6, 7, 8]:
                if len(gene) >= prefix_len:
                    prefix = gene[:prefix_len]
                    # Check if this prefix appears in multiple genes
                    matching_genes = [g for g in gene_names if g.startswith(prefix)]
                    if len(matching_genes) >= 2:  # At least 2 genes share this prefix
                        if prefix not in char_prefixes:
                            char_prefixes[prefix] = set()
                        char_prefixes[prefix].update(matching_genes)
        
        # Combine strategies: prefer prefix patterns, then pattern families, then char prefixes
        all_families = {}
        
        # Use prefix patterns first (most specific)
        for prefix, genes in prefix_patterns.items():
            if len(genes) >= 2:  # At least 2 genes in family
                for gene in genes:
                    if gene not in all_families:
                        all_families[gene] = prefix
        
        # Fill in with pattern families
        for pattern, genes in pattern_families.items():
            if len(genes) >= 2:
                for gene in genes:
                    if gene not in all_families:
                        all_families[gene] = pattern
        
        # Fill in with char prefixes (less specific, use as fallback)
        for prefix, genes in char_prefixes.items():
            if len(genes) >= 2:
                for gene in genes:
                    if gene not in all_families:
                        all_families[gene] = prefix
        
        # Assign families
        for gene in gene_names:
            if gene in all_families:
                family = all_families[gene]
                gene_families[gene] = family
                if family not in family_patterns:
                    family_patterns[family] = []
                family_patterns[family].append(gene)
            else:
                gene_families[gene] = 'Other'
                if 'Other' not in family_patterns:
                    family_patterns['Other'] = []
                family_patterns['Other'].append(gene)
        
        # Filter: only keep families with at least 2 genes (or show top families)
        min_genes_per_family = 2
        valid_families = {fam: genes for fam, genes in family_patterns.items() 
                         if len(genes) >= min_genes_per_family or fam == 'Other'}
        
        # Create row colors for gene families
        unique_families = sorted([f for f in valid_families.keys() if f != 'Other']) + (['Other'] if 'Other' in valid_families else [])
        n_families = len(unique_families)
        
        if n_families > 1:
            # Use a color palette for gene families
            family_palette = sns.color_palette("tab20", n_families) if n_families <= 20 else sns.color_palette("husl", n_families)
            family_color_dict = dict(zip(unique_families, family_palette))
            row_colors_families = pd.Series([gene_families.get(gene, 'Other') for gene in heatmap_data.index], 
                                           index=heatmap_data.index).map(family_color_dict)
            
            # Report findings
            family_counts = {fam: len(valid_families.get(fam, [])) for fam in unique_families}
            print(f"  Found {n_families} gene families:")
            for fam in unique_families[:15]:  # Show top 15
                print(f"    - {fam}: {family_counts.get(fam, 0)} genes")
            if n_families > 15:
                print(f"    ... and {n_families - 15} more families")
        else:
            row_colors_families = None
            print("  Could not identify distinct gene families (all genes grouped as 'Other')")
        
        # Use clustermap for better visualization with dendrograms
        try:
            print("  Creating clustermap with dendrograms...")
            
            # Prepare column colors for samples (if available)
            col_colors_sample = None
            if annotation_col is not None:
                unique_vals = annotation_col.iloc[:, 0].unique()
                n_colors = len(unique_vals)
                palette = sns.color_palette("Set2", n_colors)
                color_dict = dict(zip(unique_vals, palette))
                col_colors_sample = annotation_col.iloc[:, 0].map(color_dict)
            
            # Combine row colors (families) and column colors (samples)
            row_colors_list = []
            if row_colors_families is not None:
                row_colors_list.append(row_colors_families)
            
            col_colors_list = []
            if col_colors_sample is not None:
                col_colors_list.append(col_colors_sample)
            
            # Create clustermap with both row and column colors
            # Match the style from the example: dendrograms on left and top, color bars
            g = sns.clustermap(
                heatmap_data,
                cmap='RdBu_r',
                center=0,
                annot=False,
                fmt='.2f',
                cbar_kws={'label': 'Centered Normalized Counts', 'shrink': 0.6, 'aspect': 20},
                yticklabels=True,  # Show gene names to see families
                xticklabels=True,  # Show sample names
                row_colors=row_colors_list if row_colors_list else None,
                col_colors=col_colors_list if col_colors_list else None,
                figsize=(18, 16),  # Even larger for better readability
                method='ward',
                metric='euclidean',
                dendrogram_ratio=(0.15, 0.1),  # Balanced space for dendrograms
                cbar_pos=(0.02, 0.7, 0.025, 0.25),  # Position colorbar on left, better sized
                row_cluster=True,  # Cluster genes
                col_cluster=True,  # Cluster samples
                linewidths=0.3,  # Very thin lines between cells for cleaner look
                linecolor='lightgray',  # Light gray lines
                rasterized=False  # Keep vector quality
            )
            
            # Set titles and labels (matching example style - cleaner and simpler)
            g.fig.suptitle('Top 50 Most Variable Genes (Clustered by Expression Similarity)', 
                          fontsize=14, fontweight='bold', y=0.98)
            
            # Clean up dendrograms - remove all labels and ticks for minimal look
            if hasattr(g, 'ax_row_dendrogram'):
                g.ax_row_dendrogram.axis('off')  # Completely remove axis
            if hasattr(g, 'ax_col_dendrogram'):
                g.ax_col_dendrogram.axis('off')  # Completely remove axis
            
            # Format gene labels (y-axis) - make them cleaner and match example style
            if hasattr(g, 'ax_heatmap'):
                # Get current labels and format them to be more readable
                ylabels = g.ax_heatmap.get_yticklabels()
                formatted_labels = []
                
                for label in ylabels:
                    gene_name = label.get_text()
                    # For LOC IDs, show full name but make it cleaner
                    # Example: LOC108196229 stays as LOC108196229 (full name is better for identification)
                    # Just ensure it's readable
                    if len(gene_name) > 20:
                        # Only truncate if extremely long
                        formatted = gene_name[:17] + '...'
                    else:
                        formatted = gene_name
                    formatted_labels.append(formatted)
                
                # Show all labels with appropriate size - match example style
                g.ax_heatmap.set_yticklabels(formatted_labels, rotation=0, fontsize=7, ha='right', va='center')
                
                # Remove tick marks for cleaner look
                g.ax_heatmap.tick_params(axis='y', which='major', length=0, pad=2)
                
                # Format sample labels (x-axis) - make them cleaner
                xlabels = g.ax_heatmap.get_xticklabels()
                g.ax_heatmap.set_xticklabels(xlabels, rotation=45, ha='right', fontsize=9, va='top')
                g.ax_heatmap.tick_params(axis='x', which='major', length=3, pad=5)
                
                # Remove axis labels (they're self-explanatory from context)
                g.ax_heatmap.set_xlabel('', fontsize=0)
                g.ax_heatmap.set_ylabel('', fontsize=0)
            
            # Add legend for gene families if available (like example with brackets)
            if row_colors_families is not None and len(unique_families) <= 20:
                from matplotlib.patches import Patch
                # Get family counts from the data we already calculated
                # Count genes per family in the heatmap data
                family_counts_dict = {}
                for gene in heatmap_data.index:
                    fam = gene_families.get(str(gene), 'Other')
                    family_counts_dict[fam] = family_counts_dict.get(fam, 0) + 1
                
                # Create legend with family names and colors
                legend_elements = []
                for family in unique_families:
                    if family != 'Other':
                        count = family_counts_dict.get(family, 0)
                        legend_elements.append(Patch(
                            facecolor=family_color_dict[family], 
                            edgecolor='black', 
                            linewidth=0.5,
                            label=f'{family} ({count} genes)'
                        ))
                
                # Add 'Other' if it exists
                if 'Other' in unique_families:
                    other_count = family_counts_dict.get('Other', 0)
                    legend_elements.append(Patch(
                        facecolor=family_color_dict.get('Other', 'gray'),
                        edgecolor='black', 
                        linewidth=0.5,
                        label=f'Other ({other_count} genes)'
                    ))
                
                # Position legend on the right side (like example) - make it cleaner
                g.fig.legend(handles=legend_elements, 
                           title='Gene Families',
                           loc='center left', 
                           bbox_to_anchor=(1.01, 0.5),
                           fontsize=8,
                           title_fontsize=10,
                           frameon=True,
                           fancybox=False,
                           shadow=False,
                           edgecolor='black',
                           framealpha=0.9)
            
            # Save the figure with high quality settings
            plt.savefig(output_dir / "heatmap_top_variable_genes.pdf", 
                       dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
            plt.close()
            print("  Saved high-quality heatmap PDF")
            
            if annotation_col is not None:
                print(f"  Sample annotations: {annotation_col.columns[0]}")
                print(f"  Annotation values: {unique_vals.tolist()}")
            print("  Heatmap created successfully with gene family clustering")
        except Exception as e:
            print(f"  Warning: clustermap failed ({e}), using regular heatmap...")
            # Fall back to regular heatmap
            fig, ax = plt.subplots(figsize=(12, 10))
            sns.heatmap(
                heatmap_data,
                cmap='RdBu_r',
                center=0,
                annot=False,
                fmt='.2f',
                cbar_kws={'label': 'Centered Normalized Counts'},
                yticklabels=False,
                xticklabels=True,
                ax=ax
            )
            ax.set_title('Top 50 Most Variable Genes', fontsize=14, fontweight='bold')
            ax.set_xlabel('Samples', fontsize=12)
            ax.set_ylabel('Genes', fontsize=12)
            plt.tight_layout()
            plt.savefig(output_dir / "heatmap_top_variable_genes.pdf", 
                       dpi=300, bbox_inches='tight')
            plt.close()
            print("  Heatmap created successfully (fallback mode)")
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

