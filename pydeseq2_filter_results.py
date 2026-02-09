#!/usr/bin/env python3
"""
PyDESeq2 Step 2: Filter Results by Statistical Criteria
========================================================

This script filters PyDESeq2 results based on:
  - Adjusted p-value (padj) cutoff
  - Log2 fold change (log2FC) cutoff
  - Optional: direction filter (root-up only, leaf-up only, or both)

Input: pydeseq2_results_UNFILTERED.tsv (from Step 1)
Output: pydeseq2_results_FILTERED.tsv (genes passing filters)

Usage:
  # Filter for significant DEGs with large effect size
  python pydeseq2_filter_results.py pydeseq2_results_UNFILTERED.tsv \\
      --padj 0.05 --lfc 2.0 -o pydeseq2_results_FILTERED.tsv

  # Root-upregulated genes only
  python pydeseq2_filter_results.py pydeseq2_results_UNFILTERED.tsv \\
      --padj 0.05 --lfc 2.0 --root-up-only -o root_up_genes.tsv

  # No log2FC filter (any fold change)
  python pydeseq2_filter_results.py pydeseq2_results_UNFILTERED.tsv \\
      --padj 0.05 --lfc 0 -o all_significant_genes.tsv
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
from pathlib import Path


def filter_results(results_file, padj_cutoff, lfc_cutoff, root_up_only, leaf_up_only):
    """
    Filter PyDESeq2 results by statistical criteria.
    
    Args:
        results_file: Path to unfiltered results TSV
        padj_cutoff: Adjusted p-value cutoff (e.g., 0.05)
        lfc_cutoff: Absolute log2 fold change cutoff (e.g., 2.0)
        root_up_only: Keep only root-upregulated genes (log2FC > 0)
        leaf_up_only: Keep only leaf-upregulated genes (log2FC < 0)
    
    Returns:
        Filtered DataFrame
    """
    print("=" * 60)
    print("Filtering PyDESeq2 Results")
    print("=" * 60)
    
    # Read results
    print(f"\nReading results: {results_file}")
    results_df = pd.read_csv(results_file, sep='\t', index_col=0)
    print(f"Total genes in input: {len(results_df)}")
    
    # ========== STEP 1: Remove genes with missing padj (QC filter) ==========
    print("\n" + "-" * 60)
    print("STEP 1: Quality Control Filter (Remove NA padj values)")
    print("-" * 60)
    print("Genes with NA padj values are unreliable (usually very low counts)")
    
    n_before = len(results_df)
    results_df = results_df[results_df['padj'].notna()]
    n_after = len(results_df)
    n_removed = n_before - n_after
    
    print(f"\nGenes before: {n_before}")
    print(f"Genes removed (NA padj): {n_removed}")
    print(f"Genes remaining: {n_after}")
    
    if n_after == 0:
        print("\nERROR: No genes remaining after QC filter!")
        print("Check your input file for valid padj values.")
        sys.exit(1)
    
    # ========== STEP 2: Filter by adjusted p-value (statistical significance) ==========
    print("\n" + "-" * 60)
    print("STEP 2: Statistical Significance Filter (padj cutoff)")
    print("-" * 60)
    
    if padj_cutoff < 1.0:
        print(f"Keeping genes with padj < {padj_cutoff}")
        print("\nWhat this means:")
        print(f"  - padj < 0.05: < 5% false discovery rate (standard)")
        print(f"  - padj < 0.01: < 1% false discovery rate (stringent)")
        print(f"  - padj < 0.10: < 10% false discovery rate (lenient)")
        
        n_before = len(results_df)
        results_df = results_df[results_df['padj'] < padj_cutoff]
        n_after = len(results_df)
        n_removed = n_before - n_after
        
        print(f"\nGenes before: {n_before}")
        print(f"Genes passing padj filter: {n_after}")
        print(f"Genes removed: {n_removed}")
    else:
        print(f"padj filter DISABLED (cutoff = {padj_cutoff} ≥ 1.0)")
        print("Keeping all genes regardless of statistical significance")
    
    if len(results_df) == 0:
        print(f"\nWARNING: No genes pass padj < {padj_cutoff} filter!")
        print("Consider:")
        print("  1. Using a more lenient padj cutoff (e.g., 0.1)")
        print("  2. Checking if your samples have real biological differences")
        print("  3. Verifying your sample groupings in metadata")
        return results_df
    
    # ========== STEP 3: Filter by log2 fold change (biological significance) ==========
    print("\n" + "-" * 60)
    print("STEP 3: Biological Significance Filter (log2FC cutoff)")
    print("-" * 60)
    
    if lfc_cutoff > 0:
        print(f"Keeping genes with |log2FoldChange| > {lfc_cutoff}")
        print("\nWhat this means:")
        print(f"  - log2FC > {lfc_cutoff}:  ≥ {2**lfc_cutoff:.1f}× higher in root")
        print(f"  - log2FC < -{lfc_cutoff}: ≥ {2**lfc_cutoff:.1f}× higher in leaf")
        
        n_before = len(results_df)
        results_df = results_df[np.abs(results_df['log2FoldChange']) > lfc_cutoff]
        n_after = len(results_df)
        n_removed = n_before - n_after
        
        print(f"\nGenes before: {n_before}")
        print(f"Genes passing |log2FC| filter: {n_after}")
        print(f"Genes removed: {n_removed}")
    else:
        print(f"log2FC filter DISABLED (cutoff = {lfc_cutoff})")
        print("Keeping genes with any fold change (even small changes)")
    
    if len(results_df) == 0:
        print(f"\nWARNING: No genes pass |log2FC| > {lfc_cutoff} filter!")
        print("Consider using a lower log2FC cutoff (e.g., 1.0 for 2× change)")
        return results_df
    
    # ========== STEP 4: Filter by direction (optional) ==========
    print("\n" + "-" * 60)
    print("STEP 4: Direction Filter (Up/Down regulation)")
    print("-" * 60)
    
    if root_up_only:
        print("Keeping only ROOT-UPREGULATED genes (log2FC > 0)")
        n_before = len(results_df)
        results_df = results_df[results_df['log2FoldChange'] > 0]
        n_after = len(results_df)
        print(f"\nGenes before: {n_before}")
        print(f"Root-upregulated genes: {n_after}")
        print(f"Leaf-upregulated genes (removed): {n_before - n_after}")
    elif leaf_up_only:
        print("Keeping only LEAF-UPREGULATED genes (log2FC < 0)")
        n_before = len(results_df)
        results_df = results_df[results_df['log2FoldChange'] < 0]
        n_after = len(results_df)
        print(f"\nGenes before: {n_before}")
        print(f"Leaf-upregulated genes: {n_after}")
        print(f"Root-upregulated genes (removed): {n_before - n_after}")
    else:
        print("No direction filter applied")
        print("Keeping both root-upregulated and leaf-upregulated genes")
        n_up = sum(results_df['log2FoldChange'] > 0)
        n_down = sum(results_df['log2FoldChange'] < 0)
        print(f"\n  Root-upregulated: {n_up}")
        print(f"  Leaf-upregulated: {n_down}")
    
    # ========== SUMMARY ==========
    print("\n" + "=" * 60)
    print("FINAL FILTERED RESULTS SUMMARY")
    print("=" * 60)
    
    print(f"\nTotal genes passing all filters: {len(results_df)}")
    
    if len(results_df) > 0:
        # Direction breakdown
        n_up = sum(results_df['log2FoldChange'] > 0)
        n_down = sum(results_df['log2FoldChange'] < 0)
        print(f"\nDirection:")
        print(f"  Root-upregulated (log2FC > 0): {n_up}")
        print(f"  Leaf-upregulated (log2FC < 0): {n_down}")
        
        # Log2FC distribution
        print(f"\nLog2 Fold Change distribution:")
        print(f"  Min:    {results_df['log2FoldChange'].min():.2f}")
        print(f"  Q1:     {results_df['log2FoldChange'].quantile(0.25):.2f}")
        print(f"  Median: {results_df['log2FoldChange'].median():.2f}")
        print(f"  Q3:     {results_df['log2FoldChange'].quantile(0.75):.2f}")
        print(f"  Max:    {results_df['log2FoldChange'].max():.2f}")
        
        # padj distribution
        print(f"\nAdjusted p-value (padj) distribution:")
        print(f"  Min:    {results_df['padj'].min():.2e}")
        print(f"  Median: {results_df['padj'].median():.2e}")
        print(f"  Max:    {results_df['padj'].max():.2e}")
        
        # Expression level (baseMean)
        print(f"\nExpression level (baseMean) distribution:")
        print(f"  Min:    {results_df['baseMean'].min():.1f}")
        print(f"  Median: {results_df['baseMean'].median():.1f}")
        print(f"  Max:    {results_df['baseMean'].max():.1f}")
        
        # Top 10 genes
        print(f"\nTop 10 genes by |log2FoldChange|:")
        top_genes = results_df.copy()
        top_genes['abs_log2FC'] = np.abs(top_genes['log2FoldChange'])
        top_genes = top_genes.sort_values('abs_log2FC', ascending=False).head(10)
        
        print("\n{:<20} {:>12} {:>12} {:>12}".format("Gene", "log2FC", "padj", "baseMean"))
        print("-" * 60)
        for gene, row in top_genes.iterrows():
            print("{:<20} {:>12.2f} {:>12.2e} {:>12.1f}".format(
                gene, row['log2FoldChange'], row['padj'], row['baseMean']))
    
    return results_df


def main():
    parser = argparse.ArgumentParser(
        description="PyDESeq2 Step 2: Filter results by statistical criteria",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Standard filtering (padj < 0.05, |log2FC| > 2)
  python pydeseq2_filter_results.py pydeseq2_results_UNFILTERED.tsv \\
      --padj 0.05 --lfc 2.0 -o filtered_results.tsv

  # More stringent (padj < 0.01, |log2FC| > 3)
  python pydeseq2_filter_results.py pydeseq2_results_UNFILTERED.tsv \\
      --padj 0.01 --lfc 3.0 -o highly_significant.tsv

  # Root-upregulated genes only
  python pydeseq2_filter_results.py pydeseq2_results_UNFILTERED.tsv \\
      --padj 0.05 --lfc 2.0 --root-up-only -o root_up_genes.tsv

  # Any fold change (just statistical significance)
  python pydeseq2_filter_results.py pydeseq2_results_UNFILTERED.tsv \\
      --padj 0.05 --lfc 0 -o all_significant.tsv

Understanding the cutoffs:
  padj cutoffs:
    0.05 = 5% false discovery rate (standard)
    0.01 = 1% false discovery rate (stringent)
    0.10 = 10% false discovery rate (lenient)
  
  log2FC cutoffs:
    1.0 = 2× fold change (2^1 = 2)
    2.0 = 4× fold change (2^2 = 4)
    3.0 = 8× fold change (2^3 = 8)
        """
    )
    
    parser.add_argument("results_file", 
                       help="Path to unfiltered PyDESeq2 results TSV")
    
    parser.add_argument("-o", "--output", required=True,
                       help="Output file path for filtered results")
    
    parser.add_argument("--padj", type=float, default=0.05,
                       help="Adjusted p-value cutoff (default: 0.05, set 1.0 to disable)")
    
    parser.add_argument("--lfc", type=float, default=2.0,
                       help="Absolute log2 fold change cutoff (default: 2.0, set 0 to disable)")
    
    parser.add_argument("--root-up-only", action="store_true",
                       help="Keep only root-upregulated genes (log2FC > 0)")
    
    parser.add_argument("--leaf-up-only", action="store_true",
                       help="Keep only leaf-upregulated genes (log2FC < 0)")
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.results_file):
        print(f"ERROR: Results file not found: {args.results_file}")
        return 1
    
    if args.root_up_only and args.leaf_up_only:
        print("ERROR: Cannot use both --root-up-only and --leaf-up-only")
        return 1
    
    # Run filtering
    try:
        filtered_df = filter_results(
            args.results_file,
            args.padj,
            args.lfc,
            args.root_up_only,
            args.leaf_up_only
        )
        
        # Save filtered results
        if len(filtered_df) > 0:
            output_path = Path(args.output)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            
            filtered_df.to_csv(output_path, sep='\t')
            print("\n" + "=" * 60)
            print("✓ Filtered results saved!")
            print("=" * 60)
            print(f"\nOutput: {output_path}")
            print(f"Genes: {len(filtered_df)}")
            
            print("\nNext steps:")
            print("  1. Review filtered genes")
            print("  2. Generate plots:")
            print(f"     python pydeseq2_generate_plots.py {output_path} ...")
        else:
            print("\n" + "=" * 60)
            print("⚠ WARNING: No genes passed filters!")
            print("=" * 60)
            print("\nNo output file created.")
            print("Consider using more lenient cutoffs:")
            print("  --padj 0.1 (instead of 0.05)")
            print("  --lfc 1.0 (instead of 2.0)")
        
        return 0
        
    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
