#!/usr/bin/env python3
"""
Combine BLAST results with PyDESeq2 statistics

This script merges:
  - Expression statistics from PyDESeq2 (log2FC, padj, baseMean)
  - Functional annotations from BLAST (best hit per gene)

Usage:
  python scripts/combine_blast_deseq.py

Or with custom paths:
  python scripts/combine_blast_deseq.py \
      --deseq 06_analysis/pydeseq2_step1_unfiltered/pydeseq2_results_UNFILTERED.tsv \
      --blast 06_analysis/blast_input/blast_results.txt \
      --output 06_analysis/combined_results.tsv
"""

import pandas as pd
import numpy as np
import argparse
import sys
from pathlib import Path


def parse_blast_stitle(stitle):
    """
    Extract useful information from BLAST stitle field.
    
    Examples:
      "cytochrome P450 71A1 [Daucus carota]" 
        -> gene_name="CYP71A1", description="cytochrome P450 71A1"
      
      "PREDICTED: cytochrome P450 family 72 protein"
        -> gene_name="CYP72", description="cytochrome P450 family 72 protein"
    """
    if pd.isna(stitle):
        return None, None, None
    
    # Extract species (between square brackets)
    species = None
    if '[' in stitle and ']' in stitle:
        start = stitle.rfind('[')
        end = stitle.rfind(']')
        species = stitle[start+1:end]
        description = stitle[:start].strip()
    else:
        description = stitle.strip()
    
    # Extract gene name (if it looks like CYP or P450)
    gene_name = None
    parts = description.upper().split()
    for i, part in enumerate(parts):
        if 'CYP' in part or 'P450' in part:
            # Try to get the specific identifier
            if i + 1 < len(parts):
                next_part = parts[i + 1]
                # CYP71A1, CYP71, etc.
                if next_part[0].isdigit():
                    gene_name = f"CYP{next_part}"
                    break
            if 'CYP' in part:
                gene_name = part
                break
    
    return gene_name, description, species


def combine_results(deseq_file, blast_file, output_file):
    """
    Combine PyDESeq2 and BLAST results.
    """
    
    print("=" * 70)
    print("COMBINING BLAST + PyDESeq2 RESULTS")
    print("=" * 70)
    print()
    
    # ========== Step 1: Load PyDESeq2 results ==========
    print("Step 1: Loading PyDESeq2 results...")
    print(f"  File: {deseq_file}")
    
    deseq = pd.read_csv(deseq_file, sep='\t', index_col=0)
    print(f"  ✓ Loaded {len(deseq)} genes")
    
    # Basic statistics
    valid = deseq[deseq['padj'].notna()]
    print(f"    - Valid padj: {len(valid)}")
    print(f"    - NA padj: {len(deseq) - len(valid)}")
    print()
    
    # ========== Step 2: Load BLAST results ==========
    print("Step 2: Loading BLAST results...")
    print(f"  File: {blast_file}")
    
    if not Path(blast_file).exists():
        print(f"  ERROR: BLAST file not found: {blast_file}")
        print()
        print("Have you run BLAST yet?")
        print("  bash scripts/extract_sequences_for_blast.sh")
        print("  sbatch scripts/run_cyp_blast.sbatch")
        return 1
    
    blast_cols = ['qseqid', 'sseqid', 'pident', 'length', 'evalue', 'bitscore', 'stitle']
    
    try:
        blast = pd.read_csv(blast_file, sep='\t', names=blast_cols, comment='#')
    except pd.errors.EmptyDataError:
        print(f"  ERROR: BLAST file is empty")
        print("  Check your BLAST run - did it complete successfully?")
        return 1
    
    print(f"  ✓ Loaded {len(blast)} BLAST hits")
    print(f"    - Unique query genes: {blast['qseqid'].nunique()}")
    print()
    
    # ========== Step 3: Keep best hit per gene ==========
    print("Step 3: Selecting best BLAST hit per gene (lowest e-value)...")
    
    blast_best = blast.sort_values('evalue').groupby('qseqid').first().reset_index()
    print(f"  ✓ Best hits selected: {len(blast_best)} genes")
    print()
    
    # ========== Step 4: Parse BLAST annotations ==========
    print("Step 4: Parsing BLAST annotations...")
    
    blast_best[['gene_name', 'blast_description', 'blast_species']] = \
        blast_best['stitle'].apply(lambda x: pd.Series(parse_blast_stitle(x)))
    
    cyp_genes = blast_best[blast_best['stitle'].str.contains('cytochrome|P450|CYP', 
                                                              case=False, na=False)]
    print(f"  ✓ Identified {len(cyp_genes)} CYP genes from BLAST")
    print()
    
    # ========== Step 5: Merge with PyDESeq2 ==========
    print("Step 5: Merging BLAST annotations with expression data...")
    
    # Rename qseqid to gene_id for clarity
    blast_best = blast_best.rename(columns={'qseqid': 'gene_id'})
    
    # Merge (keep all PyDESeq2 genes, add BLAST where available)
    merged = deseq.merge(blast_best, left_index=True, right_on='gene_id', how='left')
    
    # Clean up columns
    merged = merged[[
        'gene_id', 'baseMean', 'log2FoldChange', 'pvalue', 'padj',
        'gene_name', 'blast_description', 'blast_species',
        'sseqid', 'pident', 'evalue', 'bitscore'
    ]]
    
    # Reorder by padj (most significant first)
    merged = merged.sort_values('padj', na_position='last')
    
    print(f"  ✓ Merged dataset: {len(merged)} genes")
    print()
    
    # ========== Step 6: Save results ==========
    print("Step 6: Saving combined results...")
    
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    merged.to_csv(output_file, sep='\t', index=False)
    print(f"  ✓ Saved: {output_file}")
    print()
    
    # ========== Step 7: Generate summary ==========
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print()
    
    total = len(merged)
    with_blast = merged['sseqid'].notna().sum()
    without_blast = total - with_blast
    
    print(f"Total genes: {total}")
    print(f"  With BLAST annotation:    {with_blast:5d} ({100*with_blast/total:5.1f}%)")
    print(f"  Without BLAST annotation: {without_blast:5d} ({100*without_blast/total:5.1f}%)")
    print()
    
    # CYP genes
    cyp_in_merged = merged[merged['gene_name'].str.contains('CYP', case=False, na=False)]
    print(f"CYP genes identified: {len(cyp_in_merged)}")
    
    if len(cyp_in_merged) > 0:
        cyp_sig = cyp_in_merged[cyp_in_merged['padj'].notna() & (cyp_in_merged['padj'] < 0.05)]
        print(f"  Significantly DE (padj < 0.05): {len(cyp_sig)}")
        
        if len(cyp_sig) > 0:
            cyp_up = cyp_sig[cyp_sig['log2FoldChange'] > 0]
            cyp_down = cyp_sig[cyp_sig['log2FoldChange'] < 0]
            print(f"    - Upregulated:   {len(cyp_up)}")
            print(f"    - Downregulated: {len(cyp_down)}")
    print()
    
    # Top annotations
    print("Top 10 most common BLAST descriptions:")
    top_desc = merged['blast_description'].value_counts().head(10)
    for desc, count in top_desc.items():
        print(f"  {count:4d}  {desc}")
    print()
    
    # Save CYP genes separately
    if len(cyp_in_merged) > 0:
        cyp_file = output_path.parent / f"{output_path.stem}_CYP_only.tsv"
        cyp_in_merged.to_csv(cyp_file, sep='\t', index=False)
        print(f"✓ CYP genes saved separately: {cyp_file}")
        print()
    
    # Save genes without BLAST hits
    no_blast = merged[merged['sseqid'].isna()]
    if len(no_blast) > 0:
        no_blast_file = output_path.parent / f"{output_path.stem}_no_blast_hits.tsv"
        no_blast[['gene_id', 'baseMean', 'log2FoldChange', 'padj']].to_csv(
            no_blast_file, sep='\t', index=False
        )
        print(f"✓ Genes without BLAST hits saved: {no_blast_file}")
        print(f"  (These might be novel genes or need different database)")
        print()
    
    print("=" * 70)
    print("NEXT STEPS")
    print("=" * 70)
    print()
    print(f"1. Review combined results:")
    print(f"   less {output_file}")
    print()
    print(f"2. Import into Geneious:")
    print(f"   - Import {output_file}")
    print(f"   - Import sequences from 06_analysis/blast_input/all_genes.faa")
    print()
    print(f"3. Filter for specific gene families:")
    print(f"   grep -i 'CYP71' {output_file}")
    print()
    print(f"4. Continue to PyDESeq2 Step 2 (filtering with annotations):")
    print(f"   sbatch scripts/run_pydeseq2_step2_filter.sbatch")
    print()
    
    return 0


def main():
    parser = argparse.ArgumentParser(
        description="Combine BLAST and PyDESeq2 results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python scripts/combine_blast_deseq.py
  
  python scripts/combine_blast_deseq.py \\
      --deseq 06_analysis/pydeseq2_step1_unfiltered/pydeseq2_results_UNFILTERED.tsv \\
      --blast 06_analysis/blast_input/blast_results.txt \\
      --output 06_analysis/combined_results.tsv
        """
    )
    
    parser.add_argument(
        "--deseq",
        default="06_analysis/pydeseq2_step1_unfiltered/pydeseq2_results_UNFILTERED.tsv",
        help="PyDESeq2 results file (default: %(default)s)"
    )
    
    parser.add_argument(
        "--blast",
        default="06_analysis/blast_input/blast_results.txt",
        help="BLAST results file (default: %(default)s)"
    )
    
    parser.add_argument(
        "-o", "--output",
        default="06_analysis/combined_results.tsv",
        help="Output file (default: %(default)s)"
    )
    
    args = parser.parse_args()
    
    # Check PyDESeq2 file exists
    if not Path(args.deseq).exists():
        print(f"ERROR: PyDESeq2 file not found: {args.deseq}")
        print()
        print("Run PyDESeq2 Step 1 first:")
        print("  sbatch scripts/run_pydeseq2_step1_analysis.sbatch")
        return 1
    
    return combine_results(args.deseq, args.blast, args.output)


if __name__ == "__main__":
    sys.exit(main())
