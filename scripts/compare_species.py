#!/usr/bin/env python3
"""
Compare differential expression between two species

Merges annotated results from two species into a single table, categorizes
each gene by its expression pattern across species, and generates summary
statistics.

Both species must have been aligned to the same reference genome so that
gene IDs (LOC numbers) are directly comparable.

Usage:
  python scripts/compare_species.py \
      --sp1 06_analysis/combined_DC/DC_swissprot_discovery_annotated.tsv \
      --sp2 06_analysis/combined_DG/DG_swissprot_discovery_annotated.tsv \
      --sp1-name DC --sp2-name DG \
      --output 06_analysis/DC_vs_DG_comparison.tsv

  # With custom significance thresholds:
  python scripts/compare_species.py \
      --sp1 DC_annotated.tsv --sp2 DG_annotated.tsv \
      --sp1-name DC --sp2-name DG \
      --output comparison.tsv \
      --padj 0.01 --lfc 2.0
"""

import pandas as pd
import numpy as np
import argparse
import sys
from pathlib import Path
from datetime import datetime


def load_annotated(filepath, species_name):
    """Load annotated results (output of combine_blast_deseq.py)."""
    print(f"Loading {species_name}: {filepath}")

    if not Path(filepath).exists():
        print(f"ERROR: File not found: {filepath}")
        sys.exit(1)

    df = pd.read_csv(filepath, sep='\t', comment='#')
    print(f"  Loaded {len(df)} genes")

    required = ['gene_id', 'log2FoldChange', 'padj']
    missing = [c for c in required if c not in df.columns]
    if missing:
        print(f"ERROR: Missing columns in {species_name}: {missing}")
        print(f"  Available: {list(df.columns)}")
        sys.exit(1)

    return df


def categorize_gene(row, padj_col_1, padj_col_2, lfc_col_1, lfc_col_2,
                    padj_cutoff, lfc_cutoff, sp1_name, sp2_name):
    """Assign a category to each gene based on DE status in both species."""
    sig_1 = (pd.notna(row[padj_col_1])
             and row[padj_col_1] <= padj_cutoff
             and abs(row[lfc_col_1]) >= lfc_cutoff)
    sig_2 = (pd.notna(row[padj_col_2])
             and row[padj_col_2] <= padj_cutoff
             and abs(row[lfc_col_2]) >= lfc_cutoff)

    if sig_1 and sig_2:
        same_dir = (row[lfc_col_1] > 0) == (row[lfc_col_2] > 0)
        if same_dir:
            return "shared_same_direction"
        return "shared_opposite_direction"

    if sig_1 and not sig_2:
        return f"{sp1_name}_only"

    if not sig_1 and sig_2:
        return f"{sp2_name}_only"

    return "not_significant"


def compare_species(sp1_file, sp2_file, sp1_name, sp2_name, output_file,
                    padj_cutoff, lfc_cutoff):
    """Merge and compare two species' annotated results."""

    print("=" * 70)
    print(f"CROSS-SPECIES COMPARISON: {sp1_name} vs {sp2_name}")
    print("=" * 70)
    print()

    # Load data
    df1 = load_annotated(sp1_file, sp1_name)
    df2 = load_annotated(sp2_file, sp2_name)
    print()

    # Rename columns with species prefix before merging
    expr_cols = ['baseMean', 'log2FoldChange', 'pvalue', 'padj']
    blast_cols = ['gene_name', 'blast_description', 'blast_species', 'sseqid',
                  'pident', 'qcovhsp', 'evalue', 'bitscore']

    rename_1 = {}
    rename_2 = {}
    for col in expr_cols:
        if col in df1.columns:
            rename_1[col] = f"{sp1_name}_{col}"
        if col in df2.columns:
            rename_2[col] = f"{sp2_name}_{col}"
    for col in blast_cols:
        if col in df1.columns:
            rename_1[col] = f"{sp1_name}_{col}"
        if col in df2.columns:
            rename_2[col] = f"{sp2_name}_{col}"

    df1 = df1.rename(columns=rename_1)
    df2 = df2.rename(columns=rename_2)

    # Keep only renamed + gene_id columns (drop alignment detail columns)
    keep_1 = ['gene_id'] + list(rename_1.values())
    keep_2 = ['gene_id'] + list(rename_2.values())
    df1 = df1[[c for c in keep_1 if c in df1.columns]]
    df2 = df2[[c for c in keep_2 if c in df2.columns]]

    # Merge on gene_id (outer join to capture all genes from both species)
    print("Merging species by gene_id (outer join)...")
    merged = df1.merge(df2, on='gene_id', how='outer')
    print(f"  Total genes in merged table: {len(merged)}")

    genes_both = merged[f"{sp1_name}_padj"].notna() & merged[f"{sp2_name}_padj"].notna()
    genes_sp1_only = merged[f"{sp1_name}_padj"].notna() & merged[f"{sp2_name}_padj"].isna()
    genes_sp2_only = merged[f"{sp1_name}_padj"].isna() & merged[f"{sp2_name}_padj"].notna()
    print(f"  Genes in both species:  {genes_both.sum()}")
    print(f"  Genes in {sp1_name} only:      {genes_sp1_only.sum()}")
    print(f"  Genes in {sp2_name} only:      {genes_sp2_only.sum()}")
    print()

    # Categorize each gene
    print(f"Categorizing genes (padj < {padj_cutoff}, |log2FC| > {lfc_cutoff})...")

    padj_1 = f"{sp1_name}_padj"
    padj_2 = f"{sp2_name}_padj"
    lfc_1 = f"{sp1_name}_log2FoldChange"
    lfc_2 = f"{sp2_name}_log2FoldChange"

    merged['category'] = merged.apply(
        categorize_gene,
        axis=1,
        padj_col_1=padj_1, padj_col_2=padj_2,
        lfc_col_1=lfc_1, lfc_col_2=lfc_2,
        padj_cutoff=padj_cutoff, lfc_cutoff=lfc_cutoff,
        sp1_name=sp1_name, sp2_name=sp2_name
    )

    # Add direction columns for readability
    def get_direction(lfc, padj, cutoff_padj, cutoff_lfc):
        if pd.isna(padj) or pd.isna(lfc):
            return "no_data"
        if padj > cutoff_padj or abs(lfc) < cutoff_lfc:
            return "ns"
        return "root_up" if lfc > 0 else "leaf_up"

    merged[f"{sp1_name}_direction"] = merged.apply(
        lambda r: get_direction(r[lfc_1], r[padj_1], padj_cutoff, lfc_cutoff), axis=1
    )
    merged[f"{sp2_name}_direction"] = merged.apply(
        lambda r: get_direction(r[lfc_2], r[padj_2], padj_cutoff, lfc_cutoff), axis=1
    )

    # Use whichever species has the BLAST annotation (prefer sp1)
    sp1_desc = f"{sp1_name}_blast_description"
    sp2_desc = f"{sp2_name}_blast_description"
    if sp1_desc in merged.columns and sp2_desc in merged.columns:
        merged['blast_description'] = merged[sp1_desc].fillna(merged[sp2_desc])
    elif sp1_desc in merged.columns:
        merged['blast_description'] = merged[sp1_desc]
    elif sp2_desc in merged.columns:
        merged['blast_description'] = merged[sp2_desc]

    sp1_gn = f"{sp1_name}_gene_name"
    sp2_gn = f"{sp2_name}_gene_name"
    if sp1_gn in merged.columns and sp2_gn in merged.columns:
        merged['gene_name'] = merged[sp1_gn].fillna(merged[sp2_gn])
    elif sp1_gn in merged.columns:
        merged['gene_name'] = merged[sp1_gn]
    elif sp2_gn in merged.columns:
        merged['gene_name'] = merged[sp2_gn]

    # Organize columns
    out_cols = ['gene_id', 'category', 'blast_description', 'gene_name']
    for sp in [sp1_name, sp2_name]:
        out_cols.extend([
            f"{sp}_log2FoldChange", f"{sp}_padj", f"{sp}_baseMean", f"{sp}_direction"
        ])
    # Add BLAST detail columns
    for sp in [sp1_name, sp2_name]:
        for bc in ['blast_species', 'sseqid', 'pident', 'evalue']:
            col = f"{sp}_{bc}"
            if col in merged.columns:
                out_cols.append(col)

    out_cols = [c for c in out_cols if c in merged.columns]
    merged = merged[out_cols]

    # Sort: significant shared genes first, then species-specific, then ns
    cat_order = {
        'shared_same_direction': 0,
        'shared_opposite_direction': 1,
        f'{sp1_name}_only': 2,
        f'{sp2_name}_only': 3,
        'not_significant': 4
    }
    merged['_sort'] = merged['category'].map(cat_order).fillna(5)
    merged = merged.sort_values(
        ['_sort', f'{sp1_name}_padj', f'{sp2_name}_padj'],
        na_position='last'
    ).drop(columns='_sort')

    # Save
    print("Saving results...")
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, 'w') as f:
        f.write("# ========================================================================\n")
        f.write(f"# Cross-Species Comparison: {sp1_name} vs {sp2_name}\n")
        f.write("# ========================================================================\n")
        f.write(f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"# {sp1_name} input: {sp1_file}\n")
        f.write(f"# {sp2_name} input: {sp2_file}\n")
        f.write(f"# Significance: padj < {padj_cutoff}, |log2FC| > {lfc_cutoff}\n")
        f.write("#\n")
        f.write("# CATEGORIES:\n")
        f.write("#   shared_same_direction     - DE in both species, same direction\n")
        f.write("#   shared_opposite_direction  - DE in both species, opposite direction\n")
        f.write(f"#   {sp1_name}_only                   - DE only in {sp1_name}\n")
        f.write(f"#   {sp2_name}_only                   - DE only in {sp2_name}\n")
        f.write("#   not_significant            - Not DE in either species\n")
        f.write("#\n")
        f.write("# DIRECTIONS:\n")
        f.write("#   root_up  = positive log2FC (higher in root)\n")
        f.write("#   leaf_up  = negative log2FC (higher in leaf)\n")
        f.write("#   ns       = not significant\n")
        f.write("#   no_data  = gene not present in that species' analysis\n")
        f.write("# ========================================================================\n")

    merged.to_csv(output_file, sep='\t', index=False, mode='a')
    print(f"  Saved: {output_file}")
    print()

    # Save subsets
    for cat in ['shared_same_direction', 'shared_opposite_direction',
                f'{sp1_name}_only', f'{sp2_name}_only']:
        subset = merged[merged['category'] == cat]
        if len(subset) > 0:
            subset_file = output_path.parent / f"{output_path.stem}_{cat}.tsv"
            subset.to_csv(subset_file, sep='\t', index=False)
            print(f"  Saved: {subset_file} ({len(subset)} genes)")

    # Summary
    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print()

    counts = merged['category'].value_counts()
    total = len(merged)

    for cat in ['shared_same_direction', 'shared_opposite_direction',
                f'{sp1_name}_only', f'{sp2_name}_only', 'not_significant']:
        n = counts.get(cat, 0)
        pct = 100 * n / total if total > 0 else 0
        print(f"  {cat:<30s}  {n:6d}  ({pct:5.1f}%)")

    print()
    print(f"  Total genes:                  {total:6d}")
    print()

    # Shared genes detail
    shared = merged[merged['category'].str.startswith('shared')]
    if len(shared) > 0:
        print("-" * 70)
        print(f"SHARED DE GENES ({len(shared)} genes)")
        print("-" * 70)
        print()

        # Direction breakdown
        same_dir = merged[merged['category'] == 'shared_same_direction']
        if len(same_dir) > 0:
            both_root = same_dir[same_dir[f'{sp1_name}_direction'] == 'root_up']
            both_leaf = same_dir[same_dir[f'{sp1_name}_direction'] == 'leaf_up']
            print(f"  Same direction:")
            print(f"    Both root-up:   {len(both_root)}")
            print(f"    Both leaf-up:   {len(both_leaf)}")

        opp_dir = merged[merged['category'] == 'shared_opposite_direction']
        if len(opp_dir) > 0:
            print(f"  Opposite direction: {len(opp_dir)}")
        print()

        # Top shared genes
        print("  Top 15 shared DE genes (by combined significance):")
        print(f"  {'gene_id':<20s} {sp1_name+'_lfc':>10s} {sp2_name+'_lfc':>10s} {'description'}")
        print(f"  {'-'*20} {'-'*10} {'-'*10} {'-'*40}")

        display = shared.head(15)
        for _, row in display.iterrows():
            gid = str(row['gene_id'])[:20]
            lfc1 = f"{row[lfc_1]:.2f}" if pd.notna(row[lfc_1]) else "NA"
            lfc2 = f"{row[lfc_2]:.2f}" if pd.notna(row[lfc_2]) else "NA"
            desc = str(row.get('blast_description', ''))[:40] if pd.notna(row.get('blast_description', np.nan)) else "unknown"
            print(f"  {gid:<20s} {lfc1:>10s} {lfc2:>10s} {desc}")
        print()

    # CYP genes
    if 'gene_name' in merged.columns:
        cyp = merged[merged['gene_name'].str.contains('CYP', case=False, na=False)]
        if len(cyp) > 0:
            print("-" * 70)
            print(f"CYP GENES ({len(cyp)} found)")
            print("-" * 70)
            cyp_cats = cyp['category'].value_counts()
            for cat, n in cyp_cats.items():
                print(f"  {cat}: {n}")
            print()

    print("=" * 70)
    print("COMPARISON COMPLETE!")
    print("=" * 70)
    print()

    return 0


def main():
    parser = argparse.ArgumentParser(
        description='Compare differential expression between two species',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('--sp1', required=True,
                        help='Species 1 annotated results (from combine_blast_deseq.py)')
    parser.add_argument('--sp2', required=True,
                        help='Species 2 annotated results (from combine_blast_deseq.py)')
    parser.add_argument('--sp1-name', required=True,
                        help='Species 1 short name (e.g., DC)')
    parser.add_argument('--sp2-name', required=True,
                        help='Species 2 short name (e.g., DG)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output file for comparison table')
    parser.add_argument('--padj', type=float, default=0.05,
                        help='Adjusted p-value cutoff (default: 0.05)')
    parser.add_argument('--lfc', type=float, default=2.0,
                        help='Log2 fold change cutoff (default: 2.0)')

    args = parser.parse_args()

    return compare_species(
        args.sp1, args.sp2,
        args.sp1_name, args.sp2_name,
        args.output,
        args.padj, args.lfc
    )


if __name__ == '__main__':
    sys.exit(main())
