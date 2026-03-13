#!/usr/bin/env python3
"""
Filter combined DESeq2 + BLAST results using multiple criteria

This script filters the output from combine_blast_deseq.py using:
  - DESeq2 statistics (padj, log2FoldChange, baseMean)
  - BLAST quality metrics (pident, qcovhsp, evalue)
  - Protein family filters (CYP, kinase, etc.)
  - BLAST hit presence/absence

Usage:
  # Standard filtering (recommended for most analyses)
  python scripts/filter_combined_results.py \
      --input 06_analysis/DC_ALL_annotated.tsv \
      --output 06_analysis/DC_filtered_standard.tsv

  # Strict filtering (high-confidence only)
  python scripts/filter_combined_results.py \
      --input 06_analysis/DC_ALL_annotated.tsv \
      --output 06_analysis/DC_filtered_strict.tsv \
      --mode strict

  # Custom filtering
  python scripts/filter_combined_results.py \
      --input 06_analysis/DC_ALL_annotated.tsv \
      --output 06_analysis/DC_filtered_custom.tsv \
      --padj 0.01 --lfc 3.0 --basemean 50 \
      --pident 70 --qcovhsp 60 --evalue 1e-10 \
      --require-blast

  # Root-upregulated genes only
  python scripts/filter_combined_results.py \
      --input 06_analysis/DC_ALL_annotated.tsv \
      --output 06_analysis/DC_root_up.tsv \
      --direction root-up
"""

import pandas as pd
import numpy as np
import argparse
import sys
from pathlib import Path


# ============================================================================
# FILTER PRESETS
# ============================================================================

FILTER_PRESETS = {
    'strict': {
        'description': 'High-confidence hits only (publication-ready)',
        'padj': 0.01,
        'lfc': 2.5,
        'basemean': 50,
        'pident': 70,
        'qcovhsp': 60,
        'evalue': 1e-10,
        'require_blast': False
    },
    'standard': {
        'description': 'Balanced filtering (recommended)',
        'padj': 0.05,
        'lfc': 2.0,
        'basemean': 10,
        'pident': 50,
        'qcovhsp': 40,
        'evalue': 1e-5,
        'require_blast': False
    },
    'lenient': {
        'description': 'Permissive filtering (exploratory analysis)',
        'padj': 0.1,
        'lfc': 1.0,
        'basemean': 5,
        'pident': 40,
        'qcovhsp': 30,
        'evalue': 1e-3,
        'require_blast': False
    },
    'blast_only': {
        'description': 'Only genes with BLAST hits (any DE significance)',
        'padj': 1.0,      # Accept all
        'lfc': 0.0,       # Accept all
        'basemean': 0,    # Accept all
        'pident': 50,
        'qcovhsp': 40,
        'evalue': 1e-5,
        'require_blast': True
    }
}


def load_combined_results(input_file):
    """Load combined DESeq2 + BLAST results."""
    print("=" * 70)
    print("LOADING DATA")
    print("=" * 70)
    print()
    print(f"Reading: {input_file}")
    
    if not Path(input_file).exists():
        print(f"ERROR: Input file not found: {input_file}")
        sys.exit(1)
    
    # Read file, skipping comment lines
    df = pd.read_csv(input_file, sep='\t', comment='#')
    
    print(f"  ✓ Loaded {len(df)} genes")
    print()
    
    return df


def apply_deseq_filters(df, padj_cutoff, lfc_cutoff, basemean_cutoff, direction=None):
    """Apply DESeq2 statistical filters."""
    print("Applying DESeq2 filters...")
    print(f"  padj < {padj_cutoff}")
    print(f"  |log2FC| > {lfc_cutoff}")
    print(f"  baseMean > {basemean_cutoff}")
    if direction:
        print(f"  Direction: {direction}")
    
    mask = (
        (df['padj'] <= padj_cutoff) &
        (df['log2FoldChange'].abs() >= lfc_cutoff) &
        (df['baseMean'] >= basemean_cutoff)
    )
    
    # Direction filter
    if direction == 'root-up':
        mask = mask & (df['log2FoldChange'] > 0)
        print("    (keeping only positive log2FC - higher in root)")
    elif direction == 'leaf-up':
        mask = mask & (df['log2FoldChange'] < 0)
        print("    (keeping only negative log2FC - higher in leaf)")
    
    filtered = df[mask].copy()
    print(f"  → {len(filtered)} genes pass DESeq2 filters")
    print()
    
    return filtered


def apply_blast_filters(df, pident_cutoff, qcovhsp_cutoff, evalue_cutoff, require_blast=False):
    """Apply BLAST quality filters."""
    print("Applying BLAST filters...")
    print(f"  pident > {pident_cutoff}%")
    print(f"  qcovhsp > {qcovhsp_cutoff}%")
    print(f"  evalue < {evalue_cutoff}")
    print(f"  Require BLAST hit: {require_blast}")
    
    if require_blast:
        # Only keep genes WITH BLAST hits
        has_blast = df['sseqid'].notna()
        df_with_blast = df[has_blast].copy()
        print(f"  → {len(df_with_blast)} genes have BLAST hits")
        
        # Apply quality filters
        mask = (
            (df_with_blast['pident'] >= pident_cutoff) &
            (df_with_blast['qcovhsp'] >= qcovhsp_cutoff) &
            (df_with_blast['evalue'] <= evalue_cutoff)
        )
        filtered = df_with_blast[mask].copy()
        print(f"  → {len(filtered)} genes pass BLAST quality filters")
    else:
        # Keep genes WITHOUT BLAST hits + genes WITH good BLAST hits
        no_blast = df['sseqid'].isna()
        has_good_blast = (
            df['sseqid'].notna() &
            (df['pident'] >= pident_cutoff) &
            (df['qcovhsp'] >= qcovhsp_cutoff) &
            (df['evalue'] <= evalue_cutoff)
        )
        filtered = df[no_blast | has_good_blast].copy()
        
        n_no_blast = no_blast.sum()
        n_good_blast = has_good_blast.sum()
        print(f"  → {n_no_blast} genes without BLAST hits (kept)")
        print(f"  → {n_good_blast} genes with good BLAST hits (kept)")
        print(f"  → {len(filtered)} total genes")
    
    print()
    return filtered


def save_results(df, output_file, filters_applied):
    """Save filtered results with metadata."""
    print("Saving results...")
    
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Write header with filter info
    from datetime import datetime
    with open(output_file, 'w') as f:
        f.write("# ========================================================================\n")
        f.write("# Filtered Combined DESeq2 + BLAST Results\n")
        f.write("# ========================================================================\n")
        f.write(f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("#\n")
        f.write("# FILTERS APPLIED:\n")
        for key, value in filters_applied.items():
            f.write(f"#   {key}: {value}\n")
        f.write("#\n")
        f.write(f"# Total genes passing filters: {len(df)}\n")
        f.write("# ========================================================================\n")
        f.write("#\n")
    
    # Append data
    df.to_csv(output_file, sep='\t', index=False, mode='a')
    print(f"  ✓ Saved: {output_file}")
    print(f"    {len(df)} genes")
    print()
    
    # Also save genes without BLAST hits separately
    no_blast = df[df['sseqid'].isna()]
    if len(no_blast) > 0:
        no_blast_file = output_path.parent / f"{output_path.stem}_no_blast_hits.tsv"
        no_blast.to_csv(no_blast_file, sep='\t', index=False)
        print(f"  ✓ Saved genes without BLAST hits: {no_blast_file}")
        print(f"    {len(no_blast)} genes")
        print()


def print_summary(df_original, df_filtered):
    """Print summary statistics."""
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print()
    
    print(f"Original genes:  {len(df_original)}")
    print(f"Filtered genes:  {len(df_filtered)}")
    print(f"Retained:        {100*len(df_filtered)/len(df_original):.1f}%")
    print()
    
    # BLAST annotation rate
    with_blast = df_filtered['sseqid'].notna().sum()
    without_blast = len(df_filtered) - with_blast
    print(f"With BLAST hit:     {with_blast} ({100*with_blast/len(df_filtered):.1f}%)")
    print(f"Without BLAST hit:  {without_blast} ({100*without_blast/len(df_filtered):.1f}%)")
    print()
    
    # Direction
    up = (df_filtered['log2FoldChange'] > 0).sum()
    down = (df_filtered['log2FoldChange'] < 0).sum()
    print(f"Upregulated (root):    {up}")
    print(f"Downregulated (leaf):  {down}")
    print()
    
    # Top protein families
    if with_blast > 0:
        print("Top 10 most common protein descriptions:")
        top_desc = df_filtered['blast_description'].value_counts().head(10)
        for desc, count in top_desc.items():
            if pd.notna(desc):
                print(f"  {count:3d}  {desc[:60]}")
        print()
    
    # Species distribution
    if with_blast > 0:
        print("Top 5 species:")
        top_species = df_filtered['blast_species'].value_counts().head(5)
        for species, count in top_species.items():
            if pd.notna(species):
                print(f"  {count:3d}  {species}")
        print()


def main():
    parser = argparse.ArgumentParser(
        description='Filter combined DESeq2 + BLAST results',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Input/output
    parser.add_argument('--input', '-i', required=True,
                        help='Input file (output from combine_blast_deseq.py)')
    parser.add_argument('--output', '-o', required=True,
                        help='Output file for filtered results')
    
    # Filter mode
    parser.add_argument('--mode', choices=['strict', 'standard', 'lenient', 'blast_only', 'custom'],
                        default='standard',
                        help='Filter preset (default: standard). Use "custom" for manual thresholds.')
    
    # DESeq2 filters
    parser.add_argument('--padj', type=float, default=None,
                        help='Adjusted p-value cutoff (default: depends on mode)')
    parser.add_argument('--lfc', type=float, default=None,
                        help='Log2 fold change cutoff (absolute value, default: depends on mode)')
    parser.add_argument('--basemean', type=float, default=None,
                        help='Minimum baseMean (default: depends on mode)')
    parser.add_argument('--direction', choices=['root-up', 'leaf-up'],
                        help='Keep only one direction (root-up = positive log2FC, leaf-up = negative)')
    
    # BLAST filters
    parser.add_argument('--pident', type=float, default=None,
                        help='Minimum percent identity (default: depends on mode)')
    parser.add_argument('--qcovhsp', type=float, default=None,
                        help='Minimum query coverage (default: depends on mode)')
    parser.add_argument('--evalue', type=float, default=None,
                        help='Maximum e-value (default: depends on mode)')
    parser.add_argument('--require-blast', action='store_true',
                        help='Only keep genes WITH BLAST hits (exclude genes without annotations)')
    
    args = parser.parse_args()
    
    # Load preset if not custom
    if args.mode != 'custom':
        preset = FILTER_PRESETS[args.mode]
        print("=" * 70)
        print(f"FILTER MODE: {args.mode.upper()}")
        print("=" * 70)
        print(f"{preset['description']}")
        print()
        
        # Use preset values if not overridden
        padj_cutoff = args.padj if args.padj is not None else preset['padj']
        lfc_cutoff = args.lfc if args.lfc is not None else preset['lfc']
        basemean_cutoff = args.basemean if args.basemean is not None else preset['basemean']
        pident_cutoff = args.pident if args.pident is not None else preset['pident']
        qcovhsp_cutoff = args.qcovhsp if args.qcovhsp is not None else preset['qcovhsp']
        evalue_cutoff = args.evalue if args.evalue is not None else preset['evalue']
        require_blast = args.require_blast or preset['require_blast']
    else:
        # Custom mode - require all parameters
        if any(x is None for x in [args.padj, args.lfc, args.basemean, args.pident, args.qcovhsp, args.evalue]):
            print("ERROR: Custom mode requires all filter parameters:")
            print("  --padj, --lfc, --basemean, --pident, --qcovhsp, --evalue")
            sys.exit(1)
        
        padj_cutoff = args.padj
        lfc_cutoff = args.lfc
        basemean_cutoff = args.basemean
        pident_cutoff = args.pident
        qcovhsp_cutoff = args.qcovhsp
        evalue_cutoff = args.evalue
        require_blast = args.require_blast
    
    # Load data
    df = load_combined_results(args.input)
    df_original = df.copy()
    
    # Apply filters
    df = apply_deseq_filters(df, padj_cutoff, lfc_cutoff, basemean_cutoff, args.direction)
    df = apply_blast_filters(df, pident_cutoff, qcovhsp_cutoff, evalue_cutoff, require_blast)
    
    # Save results
    filters_applied = {
        'mode': args.mode,
        'padj': padj_cutoff,
        'log2FC': lfc_cutoff,
        'baseMean': basemean_cutoff,
        'direction': args.direction or 'both',
        'pident': pident_cutoff,
        'qcovhsp': qcovhsp_cutoff,
        'evalue': evalue_cutoff,
        'require_blast': require_blast
    }
    
    save_results(df, args.output, filters_applied)
    
    # Print summary
    print_summary(df_original, df)
    
    print("=" * 70)
    print("FILTERING COMPLETE!")
    print("=" * 70)
    print()
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
