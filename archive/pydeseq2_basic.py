#!/usr/bin/env python3
"""
=============================================================
  PyDESeq2 Basic Analysis  
=============================================================

WHAT THIS SCRIPT 
  This script performs a differential expression analysis using PyDESeq2.
  It takes two input files:
  1. gene_count_matrix.tsv  — a table of how many times each gene was
                          detected in each sample
  2. sample_metadata.tsv      — a table that says which sample belongs
                          to which condition (root or leaf)
  
WHAT IT PRODUCES:
  results_unfiltered.tsv  — one row per gene with these columns:
      baseMean        = average expression across all samples
      log2FoldChange  = how much higher/lower in A vs B
      pvalue          = raw statistical significance
      padj            = adjusted p-value (corrected for many tests)

THE KEY CONCEPT  --  "contrast":
  Two conditions to compare are specified:
      A  = the condition you are INTERESTED IN   (e.g. root)
      B  = the BASELINE you compare against      (e.g. leaf)

  The math is:  log2FoldChange = log2( A / B )

      Positive number  -->  gene is HIGHER in A (root)
      Negative number  -->  gene is HIGHER in B (leaf)
      Zero             -->  no difference

HOW TO RUN:
  python pydeseq2_basic.py count_matrix.tsv metadata.tsv \\
      --contrast-A R --contrast-B L
"""

# ---------------------------------------------------------------
#  Load the tools (libraries) we need
# ---------------------------------------------------------------
import sys                          # lets us exit cleanly on error
import argparse                     # reads command-line arguments
import pandas as pd                 # works with tables (DataFrames)
from pathlib import Path            # handles file/folder paths

from pydeseq2.dds import DeseqDataSet   # builds the statistical model
from pydeseq2.ds  import DeseqStats     # runs the statistical test
from pydeseq2.default_inference import DefaultInference


# ---------------------------------------------------------------
# STEP 1 :  Read the input files
# ---------------------------------------------------------------
def read_data(count_file, metadata_file):
    """
    Reads two TSV files and makes sure the sample names match.

    Count matrix:  rows = genes, columns = samples
    (reference: build_count_matrix.py )

    PyDESeq2 wants the OPPOSITE:   rows = samples, columns = genes
    So we transpose it here before returning.


    metadata_file must have a 'sample' column and a 'condition' column.
    """

    

    # Read the count matrix — comes in as genes (rows) × samples (columns)
    counts = pd.read_csv(count_file, sep='\t', index_col=0)
    meta   = pd.read_csv(metadata_file, sep='\t')

   
    
    # Find sample names that appear in BOTH files
    samples_in_counts   = set(counts.columns)
    samples_in_metadata = set(meta['sample'])
    common_samples      = list(samples_in_counts & samples_in_metadata)

    if len(common_samples) == 0:
        print("ERROR: No sample names match between your two files.")
        print(f"  Count matrix columns : {list(counts.columns)}")
        print(f"  Metadata 'sample' col: {list(meta['sample'])}")
        sys.exit(1)

    # Keep only the matching samples, in the same order
    counts = counts[common_samples]
    meta   = (meta[meta['sample'].isin(common_samples)]
                  .set_index('sample')
                  .loc[common_samples])

    # Transpose: genes×samples  -->  samples×genes (what PyDESeq2 expects)
    counts = counts.T

    print(counts.shape)                          # should be (samples, genes)
    print(counts.index.equals(meta.index))       # should be True
  


    return counts, meta

# ---------------------------------------------------------------
# STEP 2 :  Run the differential expression analysis
# ---------------------------------------------------------------
def run_analysis(counts, meta, design, contrast_factor, contrast_A, contrast_B):
    """
    Asks PyDESeq2:  "For each gene, is expression different
    between condition A and condition B?"

    Returns a table with one row per gene and columns:
        baseMean, log2FoldChange, pvalue, padj
    """

    # Remove genes if they have fewer than 10 total reads across all samples
    # (counts is already samples×genes from read_data, so axis=0 sums across samples)
    keep = counts.sum(axis=0) >= 10              # axis = 0 sums across samples
    counts = counts.loc[:, keep] # keep only the genes with at least 10 reads
    print(f"  Genes with enough data (>= 10 reads): {keep.sum()}")

    # Creating the DESeqDataSet 
    inference = DefaultInference(n_cpus=4) # default inference method
    dds = DeseqDataSet(
        counts=counts_filtered, # filtered counts matrix
        metadata=meta, # metadata table
        design_factors=design, # design formula
        refit_cooks=True, # refit cooks parameter
        inference=inference, # inference method
    )

    print("  Running DESeq2 (this may take a minute)...")
    # fit size factors and dispersions and GLM coefficients. 
    dds.deseq2() 

    # --- 2e.  Set up the comparison (contrast) ---
    #   We are comparing A vs B :
    #       A = the condition we are interested in   (numerator)
    #       B = the baseline / reference             (denominator)
    contrast = (contrast_factor, contrast_A, contrast_B)

    print(f"\n  Comparison: {contrast_A}  vs  {contrast_B}")
    print(f"    Positive log2FC  -->  gene is HIGHER in {contrast_A}")
    print(f"    Negative log2FC  -->  gene is HIGHER in {contrast_B} (baseline)")

    # --- 2f.  Calculate p-values for every gene ---
    stat_res = DeseqStats(dds, contrast=contrast, n_cpus=1)
    stat_res.summary()

    return stat_res.results_df


# ---------------------------------------------------------------
# STEP 3 :  Save the results and print a summary
# ---------------------------------------------------------------
def save_and_summarise(results, contrast_A, contrast_B, output_dir):
    """Write results to a file and print a human-readable summary."""

    outfile = output_dir / "results_unfiltered.tsv"
    results.to_csv(outfile, sep='\t')
    print(f"\n  Saved all {len(results)} genes to: {outfile}")

    # --- Which genes are statistically significant? ---
    # "padj < 0.05" is the standard cutoff:
    #   it means there is less than a 5% chance the result is a false positive
    significant = results[results['padj'] < 0.05]

    print(f"\n  Significant genes (padj < 0.05): {len(significant)} out of {len(results)}")

    if len(significant) > 0:
        higher_in_A = (significant['log2FoldChange'] > 0).sum()
        higher_in_B = (significant['log2FoldChange'] < 0).sum()
        print(f"    Higher in {contrast_A} : {higher_in_A} genes")
        print(f"    Higher in {contrast_B} : {higher_in_B} genes")
    else:
        print("    (No genes passed the significance threshold)")


# ---------------------------------------------------------------
# MAIN :  Tie everything together
# ---------------------------------------------------------------
def main():
    # --- Parse command-line arguments ---
    parser = argparse.ArgumentParser(
        description="Beginner-friendly PyDESeq2 differential expression analysis"
    )

    parser.add_argument("count_matrix",
        help="Path to your gene count matrix (.tsv)")

    parser.add_argument("metadata",
        help="Path to your sample metadata (.tsv, must have 'sample' column)")

    parser.add_argument("-o", "--output",
        default="pydeseq2_basic_results",
        help="Folder to save results in (default: pydeseq2_basic_results)")

    parser.add_argument("--design",
        default="condition",
        help="Which metadata column defines your groups (default: condition)")

    parser.add_argument("--contrast-factor",
        default="condition",
        help="Same as --design for simple experiments (default: condition)")

    parser.add_argument("--contrast-A",
        default="R",
        help="Condition you are interested in / numerator (default: R = root)")

    parser.add_argument("--contrast-B",
        default="L",
        help="Baseline condition / denominator (default: L = leaf)")

    args = parser.parse_args()

    # --- Create output folder ---
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- Run the pipeline ---
    print("=" * 50)
    print("  PyDESeq2 Differential Expression Analysis")
    print("=" * 50)

    print("\n[1/3]  Reading input files...")
    counts, meta = read_data(args.count_matrix, args.metadata)

    print("\n[2/3]  Running statistical analysis...")
    results = run_analysis(
        counts, meta,
        design=args.design,
        contrast_factor=args.contrast_factor,
        contrast_A=args.contrast_A,
        contrast_B=args.contrast_B,
    )

    print("\n[3/3]  Saving results...")
    save_and_summarise(results, args.contrast_A, args.contrast_B, output_dir)

    print("\n" + "=" * 50)
    print("  Done!")
    print("=" * 50)


if __name__ == "__main__":
    sys.exit(main())
