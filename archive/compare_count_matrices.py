#!/usr/bin/env python3
"""
Compare three count matrices: DC-only, DG-only, and DC_DG combined.

Checks whether the split matrices are consistent with the combined one by
comparing sample lists, gene lists, and actual count values.

Usage:
  python scripts/compare_count_matrices.py \
      --dc  03_count_tables/00_1_DC/gene_count_matrix.tsv \
      --dg  03_count_tables/00_2_DG/gene_count_matrix.tsv \
      --combined 03_count_tables/DC_DG_Count_table/gene_count_matrix.tsv

  # Or with metadata comparison:
  python scripts/compare_count_matrices.py \
      --dc  03_count_tables/00_1_DC/gene_count_matrix.tsv \
      --dg  03_count_tables/00_2_DG/gene_count_matrix.tsv \
      --combined 03_count_tables/DC_DG_Count_table/gene_count_matrix.tsv \
      --dc-meta  03_count_tables/00_1_DC/sample_metadata.tsv \
      --dg-meta  03_count_tables/00_2_DG/sample_metadata.tsv \
      --combined-meta 03_count_tables/DC_DG_Count_table/sample.tsv
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
import numpy as np


def load_matrix(filepath, label):
    """Load a count matrix and print basic info."""
    path = Path(filepath)
    if not path.exists():
        print(f"  ERROR: {label} file not found: {filepath}")
        return None

    df = pd.read_csv(filepath, sep='\t', index_col=0)
    print(f"  {label}: {df.shape[0]} genes x {df.shape[1]} samples")
    print(f"    Samples: {sorted(df.columns.tolist())}")
    print(f"    Total counts: {df.sum().sum():,.0f}")
    print()
    return df


def compare_samples(dc, dg, combined):
    """Compare sample lists across the three matrices."""
    print("=" * 70)
    print("1. SAMPLE COMPARISON")
    print("=" * 70)
    print()

    dc_samples = set(dc.columns)
    dg_samples = set(dg.columns)
    combined_samples = set(combined.columns)
    merged_split = dc_samples | dg_samples

    print(f"  DC samples ({len(dc_samples)}):       {sorted(dc_samples)}")
    print(f"  DG samples ({len(dg_samples)}):       {sorted(dg_samples)}")
    print(f"  Combined samples ({len(combined_samples)}): {sorted(combined_samples)}")
    print()

    overlap = dc_samples & dg_samples
    if overlap:
        print(f"  WARNING: DC and DG share {len(overlap)} samples: {sorted(overlap)}")
    else:
        print(f"  OK: DC and DG have no overlapping samples")

    in_combined_not_splits = combined_samples - merged_split
    in_splits_not_combined = merged_split - combined_samples

    if in_combined_not_splits:
        print(f"  MISMATCH: Samples in combined but NOT in DC+DG: {sorted(in_combined_not_splits)}")
    if in_splits_not_combined:
        print(f"  MISMATCH: Samples in DC+DG but NOT in combined: {sorted(in_splits_not_combined)}")

    if not in_combined_not_splits and not in_splits_not_combined and not overlap:
        print(f"  PERFECT: DC ({len(dc_samples)}) + DG ({len(dg_samples)}) "
              f"= Combined ({len(combined_samples)}) samples")

    print()
    return dc_samples, dg_samples, combined_samples


def compare_genes(dc, dg, combined):
    """Compare gene lists across the three matrices."""
    print("=" * 70)
    print("2. GENE COMPARISON")
    print("=" * 70)
    print()

    dc_genes = set(dc.index)
    dg_genes = set(dg.index)
    combined_genes = set(combined.index)

    print(f"  DC genes:       {len(dc_genes)}")
    print(f"  DG genes:       {len(dg_genes)}")
    print(f"  Combined genes: {len(combined_genes)}")
    print()

    shared_dc_dg = dc_genes & dg_genes
    dc_only_genes = dc_genes - dg_genes
    dg_only_genes = dg_genes - dc_genes
    print(f"  Genes in both DC and DG:  {len(shared_dc_dg)}")
    print(f"  Genes in DC only:         {len(dc_only_genes)}")
    print(f"  Genes in DG only:         {len(dg_only_genes)}")
    print()

    all_split_genes = dc_genes | dg_genes
    in_combined_not_splits = combined_genes - all_split_genes
    in_splits_not_combined = all_split_genes - combined_genes

    if in_combined_not_splits:
        print(f"  MISMATCH: {len(in_combined_not_splits)} genes in combined but NOT in DC or DG")
        if len(in_combined_not_splits) <= 20:
            print(f"    {sorted(in_combined_not_splits)}")
        else:
            print(f"    (first 20): {sorted(in_combined_not_splits)[:20]}")
    if in_splits_not_combined:
        print(f"  MISMATCH: {len(in_splits_not_combined)} genes in DC/DG but NOT in combined")
        if len(in_splits_not_combined) <= 20:
            print(f"    {sorted(in_splits_not_combined)}")
        else:
            print(f"    (first 20): {sorted(in_splits_not_combined)[:20]}")

    if not in_combined_not_splits and not in_splits_not_combined:
        print(f"  PERFECT: Gene lists match between splits and combined")
    elif not in_combined_not_splits and in_splits_not_combined:
        print(f"\n  NOTE: The split matrices have extra genes not in combined.")
        print(f"  This can happen if build_count_matrix.py was run on different input files.")
    elif in_combined_not_splits and not in_splits_not_combined:
        print(f"\n  NOTE: The combined matrix has extra genes not in either split.")
        print(f"  This can happen if the combined was built from a different source.")

    print()
    return dc_genes, dg_genes, combined_genes


def compare_counts(dc, dg, combined):
    """Compare actual count values for shared samples/genes."""
    print("=" * 70)
    print("3. COUNT VALUE COMPARISON")
    print("=" * 70)
    print()

    issues = []

    for label, split_df in [("DC", dc), ("DG", dg)]:
        shared_samples = sorted(set(split_df.columns) & set(combined.columns))
        shared_genes = sorted(set(split_df.index) & set(combined.index))

        if not shared_samples or not shared_genes:
            print(f"  {label}: No overlapping samples/genes to compare")
            continue

        split_sub = split_df.loc[shared_genes, shared_samples]
        combined_sub = combined.loc[shared_genes, shared_samples]

        diff = split_sub - combined_sub
        mismatches = (diff != 0).sum().sum()
        total_cells = diff.shape[0] * diff.shape[1]

        print(f"  {label} vs Combined:")
        print(f"    Compared {len(shared_genes)} genes x {len(shared_samples)} samples "
              f"= {total_cells:,} values")

        if mismatches == 0:
            print(f"    PERFECT: All count values match exactly")
        else:
            pct = 100 * mismatches / total_cells
            print(f"    MISMATCH: {mismatches:,} values differ ({pct:.2f}%)")

            mismatched_genes = (diff != 0).any(axis=1)
            mismatched_samples = (diff != 0).any(axis=0)
            print(f"    Affected genes:   {mismatched_genes.sum()}")
            print(f"    Affected samples: {mismatched_samples.sum()}")
            print(f"      Samples: {sorted(mismatched_samples[mismatched_samples].index.tolist())}")

            print(f"\n    First 10 mismatches:")
            print(f"    {'gene_id':<20s} {'sample':<12s} {'split':>10s} {'combined':>10s} {'diff':>10s}")
            print(f"    {'-'*20} {'-'*12} {'-'*10} {'-'*10} {'-'*10}")
            shown = 0
            for gene in diff.index:
                for sample in diff.columns:
                    if diff.loc[gene, sample] != 0 and shown < 10:
                        sv = split_sub.loc[gene, sample]
                        cv = combined_sub.loc[gene, sample]
                        d = diff.loc[gene, sample]
                        print(f"    {str(gene):<20s} {sample:<12s} {sv:>10} {cv:>10} {d:>10}")
                        shown += 1

            issues.append((label, mismatches, total_cells))

        print()

    return issues


def compare_per_sample_totals(dc, dg, combined):
    """Compare total counts per sample across matrices."""
    print("=" * 70)
    print("4. PER-SAMPLE TOTAL COUNTS")
    print("=" * 70)
    print()

    print(f"  {'Sample':<12s} {'Split':>14s} {'Combined':>14s} {'Diff':>12s} {'Match':>8s}")
    print(f"  {'-'*12} {'-'*14} {'-'*14} {'-'*12} {'-'*8}")

    all_ok = True
    for label, split_df in [("DC", dc), ("DG", dg)]:
        shared_samples = sorted(set(split_df.columns) & set(combined.columns))
        for sample in shared_samples:
            # Use shared genes to compare apples-to-apples
            shared_genes = sorted(set(split_df.index) & set(combined.index))
            s_total = split_df.loc[shared_genes, sample].sum()
            c_total = combined.loc[shared_genes, sample].sum()
            diff = s_total - c_total
            match = "OK" if diff == 0 else "DIFF"
            if diff != 0:
                all_ok = False
            print(f"  {sample:<12s} {s_total:>14,} {c_total:>14,} {diff:>12,} {match:>8s}")

    if all_ok:
        print(f"\n  PERFECT: All per-sample totals match")
    print()


def compare_zero_genes(dc, dg, combined):
    """Compare how zero-count genes are handled."""
    print("=" * 70)
    print("5. ZERO-COUNT GENE ANALYSIS")
    print("=" * 70)
    print()

    for label, split_df in [("DC", dc), ("DG", dg)]:
        shared_genes = sorted(set(split_df.index) & set(combined.index))
        shared_samples = sorted(set(split_df.columns) & set(combined.columns))

        split_zeros = (split_df.loc[shared_genes, shared_samples] == 0).sum().sum()
        combined_zeros = (combined.loc[shared_genes, shared_samples] == 0).sum().sum()
        total = len(shared_genes) * len(shared_samples)

        print(f"  {label}:")
        print(f"    Zero values in split:    {split_zeros:>10,} / {total:,} "
              f"({100*split_zeros/total:.1f}%)")
        print(f"    Zero values in combined: {combined_zeros:>10,} / {total:,} "
              f"({100*combined_zeros/total:.1f}%)")

        genes_only_split = set(split_df.index) - set(combined.index)
        genes_only_combined = set(combined.index) - set(split_df.index)
        if genes_only_split:
            all_zero = sum(1 for g in genes_only_split
                          if split_df.loc[g].sum() == 0)
            print(f"    Genes in {label} but not combined: {len(genes_only_split)} "
                  f"({all_zero} are all-zero)")
        if genes_only_combined:
            sp_samples = sorted(set(split_df.columns) & set(combined.columns))
            if sp_samples:
                all_zero = sum(1 for g in genes_only_combined
                              if combined.loc[g, sp_samples].sum() == 0)
                print(f"    Genes in combined but not {label}: {len(genes_only_combined)} "
                      f"({all_zero} are zero for {label} samples)")
        print()


def compare_metadata(dc_meta_file, dg_meta_file, combined_meta_file):
    """Compare metadata files if provided."""
    print("=" * 70)
    print("6. METADATA COMPARISON")
    print("=" * 70)
    print()

    metas = {}
    for label, filepath in [("DC", dc_meta_file), ("DG", dg_meta_file),
                             ("Combined", combined_meta_file)]:
        if filepath and Path(filepath).exists():
            df = pd.read_csv(filepath, sep='\t')
            metas[label] = df
            print(f"  {label} metadata: {len(df)} samples, columns: {list(df.columns)}")
        elif filepath:
            print(f"  {label} metadata: FILE NOT FOUND ({filepath})")

    if not metas:
        print("  No metadata files provided or found.")
        print()
        return

    print()

    if "Combined" in metas:
        combined_meta = metas["Combined"]
        for label in ["DC", "DG"]:
            if label not in metas:
                continue
            split_meta = metas[label]

            shared_cols = set(split_meta.columns) & set(combined_meta.columns)
            print(f"  {label} vs Combined metadata:")
            print(f"    Shared columns: {sorted(shared_cols)}")

            if 'sample' in shared_cols:
                split_samples = set(split_meta['sample'])
                combined_samples = set(combined_meta['sample'])
                overlap = split_samples & combined_samples
                print(f"    Matching samples: {len(overlap)} / {len(split_samples)}")

                if overlap and len(shared_cols) > 1:
                    for col in sorted(shared_cols - {'sample'}):
                        split_vals = split_meta.set_index('sample').loc[
                            sorted(overlap), col
                        ]
                        combined_vals = combined_meta.set_index('sample').loc[
                            sorted(overlap), col
                        ]
                        mismatches = (split_vals.astype(str) != combined_vals.astype(str)).sum()
                        if mismatches > 0:
                            print(f"    Column '{col}': {mismatches} mismatches")
                        else:
                            print(f"    Column '{col}': all match")
            print()


def main():
    parser = argparse.ArgumentParser(
        description="Compare DC, DG, and DC_DG combined count matrices",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument('--dc', required=True,
                        help='DC-only count matrix (gene_count_matrix.tsv)')
    parser.add_argument('--dg', required=True,
                        help='DG-only count matrix (gene_count_matrix.tsv)')
    parser.add_argument('--combined', required=True,
                        help='Combined DC_DG count matrix (gene_count_matrix.tsv)')
    parser.add_argument('--dc-meta', default=None,
                        help='DC metadata (sample_metadata.tsv)')
    parser.add_argument('--dg-meta', default=None,
                        help='DG metadata (sample_metadata.tsv)')
    parser.add_argument('--combined-meta', default=None,
                        help='Combined metadata (sample.tsv)')

    args = parser.parse_args()

    print()
    print("#" * 70)
    print("#  COUNT MATRIX COMPARISON: DC vs DG vs DC_DG Combined")
    print("#" * 70)
    print()

    print("Loading matrices...")
    dc = load_matrix(args.dc, "DC")
    dg = load_matrix(args.dg, "DG")
    combined = load_matrix(args.combined, "Combined")

    if dc is None or dg is None or combined is None:
        print("ERROR: Could not load all three matrices. Exiting.")
        return 1

    compare_samples(dc, dg, combined)
    compare_genes(dc, dg, combined)
    compare_counts(dc, dg, combined)
    compare_per_sample_totals(dc, dg, combined)
    compare_zero_genes(dc, dg, combined)

    if args.dc_meta or args.dg_meta or args.combined_meta:
        compare_metadata(args.dc_meta, args.dg_meta, args.combined_meta)

    print("=" * 70)
    print("COMPARISON COMPLETE")
    print("=" * 70)
    print()
    print("If everything says PERFECT, your split matrices are consistent")
    print("with the combined one and you can use any of them confidently.")
    print()
    print("If there are MISMATCH warnings, check:")
    print("  - Were all three built from the same STAR count files?")
    print("  - Were the input directories set up correctly?")
    print("  - Did you re-run STAR alignment between builds?")
    print()

    return 0


if __name__ == '__main__':
    sys.exit(main())
