#!/usr/bin/env python3
"""
Extract and verify CYP + OMT gene families from TWO species combined

Takes BLAST and HMMER results from both species and produces a single merged
table showing each gene's expression and evidence across both species.
Both species must share the same reference genome (same gene IDs).

HMMER is optional for either species -- the script works with whatever
evidence is available and scores confidence accordingly.

Usage:
  # BLAST only:
  python scripts/extract_gene_families_combined.py \
      --sp1-blast 06_analysis/combined_DC/DC_swissprot_discovery_annotated.tsv \
      --sp2-blast 06_analysis/combined_DG/DG_swissprot_discovery_annotated.tsv \
      --sp1-name DC --sp2-name DG \
      --output 06_analysis/gene_families_DC_DG

  # BLAST + HMMER:
  python scripts/extract_gene_families_combined.py \
      --sp1-blast 06_analysis/combined_DC/DC_swissprot_discovery_annotated.tsv \
      --sp2-blast 06_analysis/combined_DG/DG_swissprot_discovery_annotated.tsv \
      --sp1-hmmer 06_analysis/hmmer_DC/all_genes_protein_pfam_domains.txt \
      --sp2-hmmer 06_analysis/hmmer_DG/all_genes_protein_pfam_domains.txt \
      --sp1-name DC --sp2-name DG \
      --output 06_analysis/gene_families_DC_DG
"""

import pandas as pd
import numpy as np
import argparse
import re
import sys
from pathlib import Path
from datetime import datetime

# Import the single-species search functions
sys.path.insert(0, str(Path(__file__).parent))
from extract_gene_families import (
    GENE_FAMILIES, search_blast, search_hmmer
)


def load_annotated(filepath, species_name):
    """Load annotated BLAST+DESeq2 results."""
    print(f"  Loading {species_name}: {filepath}")
    if not Path(filepath).exists():
        print(f"  ERROR: File not found: {filepath}")
        sys.exit(1)
    df = pd.read_csv(filepath, sep='\t', comment='#')
    print(f"    {len(df)} genes loaded")
    return df


def extract_family_combined(sp1_blast, sp2_blast,
                             sp1_hmmer, sp2_hmmer,
                             sp1_name, sp2_name,
                             family_name, family_def,
                             padj_cutoff=0.05, lfc_cutoff=1.0):
    """Extract a gene family from both species and merge into one table."""
    print(f"\n{'='*70}")
    print(f"Extracting: {family_name} ({family_def['full_name']}) from {sp1_name} + {sp2_name}")
    print(f"{'='*70}")

    # Search each source for each species
    print(f"\n  --- {sp1_name} ---")
    sp1_blast_ids = search_blast(sp1_blast, family_name, family_def)
    sp1_hmmer_ids, sp1_hmmer_domains = search_hmmer(sp1_hmmer, family_name, family_def)

    print(f"\n  --- {sp2_name} ---")
    sp2_blast_ids = search_blast(sp2_blast, family_name, family_def)
    sp2_hmmer_ids, sp2_hmmer_domains = search_hmmer(sp2_hmmer, family_name, family_def)

    # Union of all gene IDs across both species
    all_ids = (sp1_blast_ids | sp1_hmmer_ids |
               sp2_blast_ids | sp2_hmmer_ids)
    print(f"\n  Union across both species: {len(all_ids)} unique genes")

    if not all_ids:
        return pd.DataFrame()

    # Build lookup tables for expression data
    sp1_lookup = sp1_blast.set_index('gene_id') if 'gene_id' in sp1_blast.columns else sp1_blast
    sp2_lookup = sp2_blast.set_index('gene_id') if 'gene_id' in sp2_blast.columns else sp2_blast

    rows = []
    for gid in sorted(all_ids):
        row = {'gene_id': gid, 'gene_family': family_name}

        # Species 1 expression stats
        if gid in sp1_lookup.index:
            gene_data = sp1_lookup.loc[gid]
            if isinstance(gene_data, pd.DataFrame):
                gene_data = gene_data.iloc[0]
            for col in ['baseMean', 'log2FoldChange', 'pvalue', 'padj']:
                if col in gene_data.index:
                    row[f'{sp1_name}_{col}'] = gene_data[col]
            for col in ['blast_description', 'blast_species', 'pident', 'evalue']:
                if col in gene_data.index:
                    row[f'{sp1_name}_{col}'] = gene_data[col]

        # Species 2 expression stats
        if gid in sp2_lookup.index:
            gene_data = sp2_lookup.loc[gid]
            if isinstance(gene_data, pd.DataFrame):
                gene_data = gene_data.iloc[0]
            for col in ['baseMean', 'log2FoldChange', 'pvalue', 'padj']:
                if col in gene_data.index:
                    row[f'{sp2_name}_{col}'] = gene_data[col]
            for col in ['blast_description', 'blast_species', 'pident', 'evalue']:
                if col in gene_data.index:
                    row[f'{sp2_name}_{col}'] = gene_data[col]

        # Consolidated BLAST description (prefer sp1)
        desc1 = row.get(f'{sp1_name}_blast_description', np.nan)
        desc2 = row.get(f'{sp2_name}_blast_description', np.nan)
        row['blast_description'] = desc1 if pd.notna(desc1) else desc2

        # Evidence tracking per species
        sp1_sources = []
        if gid in sp1_blast_ids:
            sp1_sources.append('BLAST')
        if gid in sp1_hmmer_ids:
            sp1_sources.append('HMMER')

        sp2_sources = []
        if gid in sp2_blast_ids:
            sp2_sources.append('BLAST')
        if gid in sp2_hmmer_ids:
            sp2_sources.append('HMMER')

        row[f'{sp1_name}_evidence'] = '+'.join(sp1_sources) if sp1_sources else 'none'
        row[f'{sp2_name}_evidence'] = '+'.join(sp2_sources) if sp2_sources else 'none'

        all_sources = set(sp1_sources + sp2_sources)
        row['evidence_sources'] = '+'.join(sorted(all_sources)) if all_sources else 'none'
        row['evidence_count'] = len(all_sources)
        row['confidence'] = 'high' if len(all_sources) >= 2 else 'low'

        # Domain details (merged across species)
        hmmer_parts = []
        if gid in sp1_hmmer_domains:
            hmmer_parts.append(sp1_hmmer_domains[gid])
        if gid in sp2_hmmer_domains:
            hmmer_parts.append(sp2_hmmer_domains[gid])
        row['hmmer_domains'] = '; '.join(set('; '.join(hmmer_parts).split('; '))) if hmmer_parts else ''

        # Direction per species
        for sp_name in [sp1_name, sp2_name]:
            lfc = row.get(f'{sp_name}_log2FoldChange', np.nan)
            padj = row.get(f'{sp_name}_padj', np.nan)
            if pd.notna(padj) and pd.notna(lfc) and padj <= padj_cutoff and abs(lfc) >= lfc_cutoff:
                row[f'{sp_name}_direction'] = 'root_up' if lfc > 0 else 'leaf_up'
            elif pd.notna(padj):
                row[f'{sp_name}_direction'] = 'ns'
            else:
                row[f'{sp_name}_direction'] = 'no_data'

        # Category (shared vs species-specific)
        sp1_sig = row.get(f'{sp1_name}_direction', 'no_data') in ('root_up', 'leaf_up')
        sp2_sig = row.get(f'{sp2_name}_direction', 'no_data') in ('root_up', 'leaf_up')

        if sp1_sig and sp2_sig:
            same_dir = (row[f'{sp1_name}_direction'] == row[f'{sp2_name}_direction'])
            row['category'] = 'shared_same_direction' if same_dir else 'shared_opposite_direction'
        elif sp1_sig:
            row['category'] = f'{sp1_name}_only'
        elif sp2_sig:
            row['category'] = f'{sp2_name}_only'
        else:
            row['category'] = 'not_significant'

        rows.append(row)

    result = pd.DataFrame(rows)

    # Sort: shared first, then species-specific, then ns
    cat_order = {
        'shared_same_direction': 0,
        'shared_opposite_direction': 1,
        f'{sp1_name}_only': 2,
        f'{sp2_name}_only': 3,
        'not_significant': 4
    }
    result['_sort'] = result['category'].map(cat_order).fillna(5)
    result = result.sort_values(
        ['_sort', f'{sp1_name}_padj', f'{sp2_name}_padj'],
        na_position='last'
    ).drop(columns='_sort')

    # Summary
    print(f"\n  Results: {len(result)} {family_name} genes")
    for cat, count in result['category'].value_counts().items():
        print(f"    {cat}: {count}")
    for conf in ['high', 'medium', 'low']:
        n = (result['confidence'] == conf).sum()
        if n > 0:
            print(f"    {conf} confidence: {n}")

    return result


def main():
    parser = argparse.ArgumentParser(
        description='Extract CYP/OMT gene families from two species combined',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('--sp1-blast', required=True,
                        help='Species 1 combined annotated file')
    parser.add_argument('--sp2-blast', required=True,
                        help='Species 2 combined annotated file')
    parser.add_argument('--sp1-hmmer', default=None,
                        help='Species 1 HMMER domtblout file (optional)')
    parser.add_argument('--sp2-hmmer', default=None,
                        help='Species 2 HMMER domtblout file (optional)')
    parser.add_argument('--sp1-name', required=True, help='Species 1 code (e.g., DC)')
    parser.add_argument('--sp2-name', required=True, help='Species 2 code (e.g., DG)')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    parser.add_argument('--families', nargs='+', default=['CYP', 'OMT'],
                        choices=list(GENE_FAMILIES.keys()),
                        help='Which families to extract (default: CYP OMT)')
    parser.add_argument('--padj', type=float, default=0.05,
                        help='padj cutoff for DE significance (default: 0.05)')
    parser.add_argument('--lfc', type=float, default=1.0,
                        help='|log2FC| cutoff for DE significance (default: 1.0)')

    args = parser.parse_args()

    sp1 = args.sp1_name
    sp2 = args.sp2_name

    print("=" * 70)
    print(f"COMBINED GENE FAMILY EXTRACTION: {sp1} + {sp2}")
    print("=" * 70)
    print()

    # Load BLAST files
    sp1_blast = load_annotated(args.sp1_blast, sp1)
    sp2_blast = load_annotated(args.sp2_blast, sp2)

    # Check optional inputs
    def check_optional(filepath, label):
        if filepath and Path(filepath).exists():
            print(f"  {label}: found")
            return filepath
        if filepath:
            print(f"  {label}: not found (skipping)")
        else:
            print(f"  {label}: not provided")
        return None

    sp1_hmmer = check_optional(args.sp1_hmmer, f"{sp1} HMMER")
    sp2_hmmer = check_optional(args.sp2_hmmer, f"{sp2} HMMER")

    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Extract each family
    all_results = []
    for family_name in args.families:
        family_def = GENE_FAMILIES[family_name]
        result = extract_family_combined(
            sp1_blast, sp2_blast,
            sp1_hmmer, sp2_hmmer,
            sp1, sp2,
            family_name, family_def,
            padj_cutoff=args.padj, lfc_cutoff=args.lfc
        )

        if not result.empty:
            family_file = out_dir / f"{sp1}_{sp2}_{family_name}_verified.tsv"
            result.to_csv(family_file, sep='\t', index=False)
            print(f"\n  Saved: {family_file} ({len(result)} genes)")
            all_results.append(result)

    # Save combined
    if all_results:
        combined = pd.concat(all_results, ignore_index=True)
        combined_file = out_dir / f"{sp1}_{sp2}_CYP_OMT_combined.tsv"

        with open(combined_file, 'w') as f:
            f.write("# ========================================================================\n")
            f.write(f"# Combined Gene Family Extraction: {sp1} + {sp2}\n")
            f.write("# ========================================================================\n")
            f.write(f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"# {sp1} BLAST: {args.sp1_blast}\n")
            f.write(f"# {sp2} BLAST: {args.sp2_blast}\n")
            f.write(f"# {sp1} HMMER: {sp1_hmmer or 'not provided'}\n")
            f.write(f"# {sp2} HMMER: {sp2_hmmer or 'not provided'}\n")
            f.write(f"# DE thresholds: padj < {args.padj}, |log2FC| > {args.lfc}\n")
            f.write("#\n")
            f.write("# CATEGORIES:\n")
            f.write("#   shared_same_direction     - DE in both species, same direction\n")
            f.write("#   shared_opposite_direction  - DE in both species, opposite direction\n")
            f.write(f"#   {sp1}_only                - DE only in {sp1}\n")
            f.write(f"#   {sp2}_only                - DE only in {sp2}\n")
            f.write("#   not_significant            - Not DE in either species\n")
            f.write("#\n")
            f.write("# CONFIDENCE: high (BLAST+HMMER agree), low (1 source only)\n")
            f.write("# NOTE: blast_species shows the closest characterized homolog,\n")
            f.write("#   not contamination. Arabidopsis hits = your gene matches\n")
            f.write("#   a known Arabidopsis protein in the database.\n")
            f.write("# ========================================================================\n")

        combined.to_csv(combined_file, sep='\t', index=False, mode='a')
        print(f"\n  Saved combined: {combined_file} ({len(combined)} genes)")

        # Save subsets
        for cat in ['shared_same_direction', 'shared_opposite_direction',
                     f'{sp1}_only', f'{sp2}_only']:
            subset = combined[combined['category'] == cat]
            if len(subset) > 0:
                subset_file = out_dir / f"{sp1}_{sp2}_CYP_OMT_{cat}.tsv"
                subset.to_csv(subset_file, sep='\t', index=False)
                print(f"  Saved: {subset_file} ({len(subset)} genes)")

    # Final summary
    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print()

    if not all_results:
        print("  No gene family members found.")
        return 0

    combined = pd.concat(all_results, ignore_index=True)

    for fam in args.families:
        fam_df = combined[combined['gene_family'] == fam]
        if fam_df.empty:
            continue

        print(f"  {fam} ({GENE_FAMILIES[fam]['full_name']}): {len(fam_df)} genes")
        for cat, count in fam_df['category'].value_counts().items():
            print(f"    {cat}: {count}")
        print()

    # Summary file
    summary_file = out_dir / f"{sp1}_{sp2}_gene_family_summary.txt"
    with open(summary_file, 'w') as f:
        f.write(f"Combined Gene Family Summary: {sp1} + {sp2}\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        for fam in args.families:
            fam_df = combined[combined['gene_family'] == fam]
            f.write(f"{fam} ({GENE_FAMILIES[fam]['full_name']}): {len(fam_df)} genes\n")
            for cat, count in fam_df['category'].value_counts().items():
                f.write(f"  {cat}: {count}\n")
            for conf in ['high', 'medium', 'low']:
                n = (fam_df['confidence'] == conf).sum()
                f.write(f"  {conf} confidence: {n}\n")
            f.write("\n")

    print(f"  Summary: {summary_file}")
    print()
    print("=" * 70)
    print("EXTRACTION COMPLETE!")
    print("=" * 70)

    return 0


if __name__ == '__main__':
    sys.exit(main())
