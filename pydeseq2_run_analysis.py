#!/usr/bin/env python3
"""
PyDESeq2 Step 1: Run Statistical Analysis (NO FILTERING)
=========================================================

This script performs differential expression analysis using PyDESeq2
and saves ALL results without filtering.

Output: pydeseq2_results_UNFILTERED.tsv
  - Contains ALL genes with their statistics (including "bad" p-values)
  - Use this for BLAST and Geneious annotation
  - Filter later with pydeseq2_filter_results.py

Usage:
  python pydeseq2_run_analysis.py count_matrix.tsv metadata.tsv -o results \\
      --contrast-factor condition --contrast-A root --contrast-B leaf
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
    from pydeseq2.default_inference import DefaultInference
    import matplotlib.pyplot as plt
except ImportError as e:
    print(f"ERROR: Missing required package: {e}")
    print("Please install PyDESeq2 and dependencies:")
    print("  pip install pydeseq2 matplotlib")
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
    print("\nGenerating QC plots...")
    
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
    print(f"  Saved: {output_dir / 'qc_total_counts.pdf'}")


def generate_log2foldchange_histogram(results_df, output_dir):
    """Generate histogram of log2foldchange values."""
    print("Generating log2foldchange histogram...")
    plt.hist(results_df['log2FoldChange'], bins=100, edgecolor='black', alpha=0.5)
    plt.xlabel('log2foldchange')
    plt.ylabel('Number of Genes')
    plt.title('Distribution of log2foldchange values')
    plt.savefig(output_dir / "log2foldchange_histogram.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_dir / 'log2foldchange_histogram.pdf'}")

def generate_padj_histogram(results_df, output_dir):
    """Generate histogram of padj values."""
    print("Generating padj histogram...")
    plt.hist(results_df['padj'], bins=100, edgecolor='black', alpha=0.5)
    plt.xlabel('padj')
    plt.ylabel('Number of Genes')
    plt.title('Distribution of padj values')
    plt.savefig(output_dir / "padj_histogram.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_dir / 'padj_histogram.pdf'}")


def run_pydeseq2_analysis(count_matrix, metadata, design_formula, output_dir,
                          contrast_factor=None, contrast_A=None, contrast_B=None):
    """
    Run PyDESeq2 differential expression analysis.
    
    This function performs the 5 main PyDESeq2 steps:
    1. Size factor calculation (normalize for library size)
    2. Dispersion estimation (measure biological variability)
    3. Statistical testing (calculate p-values)
    4. Multiple testing correction (calculate padj)
    5. Log2 fold change calculation
    
    Returns ALL results without filtering.
    """
    print("\n" + "=" * 60)
    print("Running PyDESeq2 Differential Expression Analysis")
    print("=" * 60)
    print(f"Design formula: {design_formula}")
    
    # ========== STEP 0: Prepare data ==========
    print("\nStep 0: Preparing data...")
    count_matrix_t = count_matrix.T  # Transpose to (samples x genes)
    print(f"  Count matrix: {count_matrix_t.shape[0]} samples × {count_matrix_t.shape[1]} genes")
    
    metadata_indexed = metadata.copy()
    if metadata_indexed.index.name != 'sample' and 'sample' in metadata_indexed.columns:
        metadata_indexed = metadata_indexed.set_index('sample')
    metadata_indexed = metadata_indexed.loc[count_matrix_t.index]
    
    # Filter low count genes (< 10 total reads)
    min_counts = 10
    genes_to_keep = (count_matrix_t.sum(axis=0) >= min_counts)
    print(f"  Genes before filtering: {len(count_matrix_t.columns)}")
    print(f"  Genes with ≥ {min_counts} counts: {genes_to_keep.sum()}")
    
    count_matrix_filtered = count_matrix_t.loc[:, genes_to_keep]
    
    # Create DeseqDataSet
    print("\nCreating DeseqDataSet...")
    dds = DeseqDataSet(
        counts=count_matrix_filtered,
        metadata=metadata_indexed,
        design_factors=design_formula.split('+'),
        refit_cooks=True,
        n_cpus=1
    )
    
    # ========== STEP 1: Calculate Size Factors ==========
    print("\n" + "=" * 60)
    print("STEP 1: Calculating Size Factors (Library Size Normalization)")
    print("=" * 60)
    print("This adjusts for different sequencing depths between samples.")
    print("Example: If sample A has 2× more total reads than sample B,")
    print("         its size factor will be ~2.0")
    
    if hasattr(dds, 'deseq2'):
        dds.deseq2()
    elif hasattr(dds, 'fit'):
        dds.fit()
    else:
        # Step-by-step workflow
        dds.estimate_size_factors()
        print("\n  Size factors calculated:")
        if hasattr(dds, 'size_factors') and dds.size_factors is not None:
            for sample, sf in zip(dds.obs.index, dds.size_factors):
                print(f"    {sample}: {sf:.3f}")
        
        # ========== STEP 2: Estimate Dispersions ==========
        print("\n" + "=" * 60)
        print("STEP 2: Estimating Dispersions (Biological Variability)")
        print("=" * 60)
        print("Dispersion measures how much replicates vary from each other.")
        print("Lower dispersion = more consistent replicates = more reliable")
        
        dds.estimate_dispersions()
        
        # ========== STEP 3-5: Fit model and calculate statistics ==========
        print("\n" + "=" * 60)
        print("STEP 3-5: Fitting Model and Calculating Statistics")
        print("=" * 60)
        print("  - Statistical testing (p-values)")
        print("  - Multiple testing correction (adjusted p-values)")
        print("  - Log2 fold change calculation")
        
        dds.fit()
    
    # ========== Get Statistical Results ==========
    print("\nComputing statistical test results...")
    
    # Determine contrast
    # Contrast_B: baseline condition (leaf)
    # Contrast_A: measured condition (root)
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
        print(f"\nContrast: {factor}: {contrast_A} vs {contrast_B}")
        print(f"  (Positive log2FC = higher in {contrast_A})")
        print(f"  (Negative log2FC = higher in {contrast_B})")
    else:
        contrast = None
        print("\nUsing default contrast (auto-detected)")
    
    # Create DeseqStats object and get results
    stat_res = DeseqStats(dds, contrast=contrast, n_cpus=1)
    stat_res.summary()
    
    results_df = stat_res.results_df
    
    # ========== Save UNFILTERED Results ==========
    print("\n" + "=" * 60)
    print("Saving UNFILTERED Results (ALL genes)")
    print("=" * 60)
    
    results_file = output_dir / "pydeseq2_results_UNFILTERED.tsv"
    results_df.to_csv(results_file, sep='\t')
    print(f"\n✓ Saved: {results_file}")
    print(f"  Total genes: {len(results_df)}")

    # ========== Save Gene Candidate List (for BLAST) ==========
    # One gene ID per line — genes with valid padj only (no padj/log2FC filter yet)
    valid_genes = results_df[results_df['padj'].notna()].index
    gene_list_file = output_dir / "all_gene_ids.txt"
    pd.Series(valid_genes).to_csv(gene_list_file, index=False, header=False)
    print(f"\n✓ Saved gene candidate list: {gene_list_file}")
    print(f"  Genes with valid padj: {len(valid_genes)}")
    
    # Print summary statistics
    print("\n" + "-" * 60)
    print("RESULTS SUMMARY (before filtering)")
    print("-" * 60)
    
    print(f"\nTotal genes analyzed: {len(results_df)}")
    
    valid_results = results_df[results_df['padj'].notna()]
    print(f"Genes with valid padj: {len(valid_results)}")
    
    sig_005 = valid_results[valid_results['padj'] < 0.05]
    print(f"\nGenes with padj < 0.05: {len(sig_005)}")
    if len(sig_005) > 0:
        up = sig_005[sig_005['log2FoldChange'] > 0]
        down = sig_005[sig_005['log2FoldChange'] < 0]
        print(f"  Upregulated in {contrast_A if contrast_A else 'A'}: {len(up)}")
        print(f"  Upregulated in {contrast_B if contrast_B else 'B'}: {len(down)}")
    
    sig_001 = valid_results[valid_results['padj'] < 0.01]
    print(f"\nGenes with padj < 0.01: {len(sig_001)}")
    
    high_lfc = valid_results[np.abs(valid_results['log2FoldChange']) > 2]
    print(f"Genes with |log2FC| > 2: {len(high_lfc)}")
    
    # baseMean statistics
    print(f"\nExpression level (baseMean) distribution:")
    print(f"  Min: {results_df['baseMean'].min():.1f}")
    print(f"  Q1:  {results_df['baseMean'].quantile(0.25):.1f}")
    print(f"  Median: {results_df['baseMean'].median():.1f}")
    print(f"  Q3:  {results_df['baseMean'].quantile(0.75):.1f}")
    print(f"  Max: {results_df['baseMean'].max():.1f}")
    
    print("\n" + "=" * 60)
    print("NEXT STEPS:")
    print("=" * 60)
    print("1. Use this file for BLAST annotation:")
    print(f"   Input: {results_file}")
    print("")
    print("2. Filter results by statistical criteria:")
    print(f"   python pydeseq2_filter_results.py {results_file} \\")
    print("       --padj 0.05 --lfc 2.0 -o filtered_results.tsv")
    print("")
    print("3. Generate plots:")
    print("   python pydeseq2_generate_plots.py filtered_results.tsv ...")
    
    return dds, stat_res, results_df


def main():
    parser = argparse.ArgumentParser(
        description="PyDESeq2 Step 1: Run statistical analysis (no filtering)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python pydeseq2_run_analysis.py count_matrix.tsv metadata.tsv \\
      -o results \\
      --contrast-factor condition \\
      --contrast-A root \\
      --contrast-B leaf

Output:
  results/pydeseq2_results_UNFILTERED.tsv - ALL genes with statistics
  results/qc_total_counts.pdf - Quality control plot
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
    parser.add_argument("--contrast-A", default="R",
                       help="Numerator condition - positive log2FC (default: R = root)")
    parser.add_argument("--contrast-B", default="L",
                       help="Denominator condition (default: L = leaf)")
    
    args = parser.parse_args()
    
    # Validate inputs
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
    print("PyDESeq2 Statistical Analysis (Step 1: No Filtering)")
    print("=" * 60)
    print(f"Count matrix: {args.count_matrix}")
    print(f"Metadata: {args.metadata}")
    print(f"Output directory: {output_dir}")
    print()
    
    try:
        # Read data
        count_matrix, metadata = read_data(args.count_matrix, args.metadata)
        
        # Determine design formula
        if args.design:
            design_formula = args.design
        else:
            if args.contrast_factor in metadata.columns:
                design_formula = args.contrast_factor
            elif 'condition' in metadata.columns:
                design_formula = 'condition'
            elif 'treatment' in metadata.columns:
                design_formula = 'treatment'
            else:
                categorical_cols = metadata.select_dtypes(include=['object', 'category']).columns
                if len(categorical_cols) > 0:
                    design_formula = categorical_cols[0]
                    print(f"Auto-detected design factor: {design_formula}")
                else:
                    raise ValueError("No suitable grouping variable found in metadata.")
        
        print(f"Using design formula: {design_formula}\n")
        
        # Generate QC plots
        generate_qc_plots(count_matrix, metadata, output_dir)
        
        # Run PyDESeq2 analysis
        dds, stat_res, results_df = run_pydeseq2_analysis(
            count_matrix, metadata, design_formula, output_dir,
            contrast_factor=args.contrast_factor,
            contrast_A=args.contrast_A,
            contrast_B=args.contrast_B
        )

        # Generate histograms to visualize padj and log2FC distributions
        generate_log2foldchange_histogram(results_df, output_dir)
        generate_padj_histogram(results_df, output_dir)
        
        print("\n" + "=" * 60)
        print("✓ Analysis complete!")
        print("=" * 60)
        
        return 0
        
    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
