#!/usr/bin/env python3
"""
Merge count matrices and metadata from multiple species directories.

PURPOSE:
  Combines gene_count_matrix.tsv and sample_metadata.tsv from two (or more)
  species directories into a single set of files. This lets you run PyDESeq2
  on the combined dataset (e.g., DC + DG together for Root vs Leaf).

  Both species must have been aligned to the SAME genome so gene IDs match.
  Genes present in one matrix but not the other get filled with 0 counts.

USAGE:
  python merge_count_matrices.py \
      --dirs 03_count_tables/00_1_DC 03_count_tables/00_2_DG \
      --output-dir 03_count_tables/00_1_2_DC_DG

INPUT FILES (per directory):
  gene_count_matrix.tsv   — genes (rows) × samples (columns), tab-separated
  sample_metadata.tsv     — sample name + condition columns, tab-separated

OUTPUT FILES:
  gene_count_matrix.tsv   — merged count matrix (all samples, union of genes)
  sample_metadata.tsv     — merged metadata (all samples)

EXAMPLE INPUT (DC, 12 samples):
  Geneid     DC1L1  DC1L2  DC1R1  DC1R2 ...
  LOC108192  150    200    300    400   ...

EXAMPLE INPUT (DG, 6 samples):
  Geneid     DG1L1  DG1L2  DG1R1  DG1R2 ...
  LOC108192  120    180    280    350   ...

EXAMPLE OUTPUT (combined, 18 samples):
  Geneid     DC1L1  DC1L2  DC1R1  DC1R2 ... DG1L1  DG1L2  DG1R1 ...
  LOC108192  150    200    300    400   ... 120    180    280   ...
"""

import argparse
import sys
from pathlib import Path

try:
    import pandas as pd
except ImportError:
    print("ERROR: pandas is required. Install with: pip install pandas")
    sys.exit(1)


def merge_counts(dirs, output_dir):
    """
    Merge gene_count_matrix.tsv files from multiple directories.

    Uses an OUTER join on gene IDs — genes missing from one dataset
    get 0 counts (they weren't expressed / detected in that dataset).
    """
    print("=" * 60)
    print("Merging count matrices")
    print("=" * 60)

    all_counts = []
    all_metadata = []

    for d in dirs:
        d = Path(d)
        counts_file = d / "gene_count_matrix.tsv"
        meta_file = d / "sample_metadata.tsv"

        # ── Load count matrix ──
        if not counts_file.exists():
            print(f"  ERROR: {counts_file} not found — skipping {d.name}")
            continue

        counts = pd.read_csv(counts_file, sep='\t', index_col=0)

        # Drop featureCounts annotation columns if present
        annotation_cols = {'Chr', 'Start', 'End', 'Strand', 'Length'}
        drop_cols = [c for c in counts.columns if c in annotation_cols]
        if drop_cols:
            counts = counts.drop(columns=drop_cols)

        counts = counts.apply(pd.to_numeric, errors='coerce').fillna(0).astype(int)
        print(f"  {d.name}: {counts.shape[0]} genes × {counts.shape[1]} samples")
        print(f"    Samples: {list(counts.columns)}")
        all_counts.append(counts)

        # ── Load metadata ──
        if not meta_file.exists():
            print(f"  WARNING: {meta_file} not found — building from column names")
            meta = pd.DataFrame({'sample': counts.columns})
            meta['condition'] = meta['sample'].apply(
                lambda s: 'R' if 'R' in s.upper() else 'L'
            )
            meta = meta.set_index('sample')
        else:
            meta = pd.read_csv(meta_file, sep='\t')
            if 'sample' in meta.columns:
                meta = meta.set_index('sample')
        all_metadata.append(meta)

    if len(all_counts) < 2:
        print("ERROR: Need at least 2 directories with valid count matrices.")
        sys.exit(1)

    # ── Merge counts (outer join — keep ALL genes from both) ──
    merged = all_counts[0]
    for other in all_counts[1:]:
        # Check for overlapping sample names
        overlap_samples = set(merged.columns) & set(other.columns)
        if overlap_samples:
            print(f"  WARNING: Overlapping sample names: {overlap_samples}")
            print(f"  These samples appear in multiple directories.")

        merged = merged.join(other, how='outer')

    # Fill NaN with 0 (genes not in one dataset get zero counts)
    merged = merged.fillna(0).astype(int)

    # ── Merge metadata ──
    merged_meta = pd.concat(all_metadata)

    # Only keep samples that exist in the count matrix
    common = sorted(set(merged.columns) & set(merged_meta.index))
    merged = merged[common]
    merged_meta = merged_meta.loc[common]

    # ── Report ──
    print()
    print(f"  Combined matrix: {merged.shape[0]} genes × {merged.shape[1]} samples")

    gene_sets = [set(c.index) for c in all_counts]
    shared = gene_sets[0]
    all_genes = gene_sets[0]
    for gs in gene_sets[1:]:
        shared = shared & gs
        all_genes = all_genes | gs
    print(f"  Genes shared by all: {len(shared)}")
    print(f"  Genes in union:      {len(all_genes)}")
    only_counts = [len(gs - shared) for gs in gene_sets]
    for i, d in enumerate(dirs):
        if only_counts[i] > 0:
            print(f"  Genes unique to {Path(d).name}: {only_counts[i]} (filled with 0)")

    print(f"  Samples: {list(merged.columns)}")
    print()
    print(f"  Conditions in metadata:")
    if 'condition' in merged_meta.columns:
        for cond, count in merged_meta['condition'].value_counts().items():
            print(f"    {cond}: {count} samples")
    print()

    # ── Save ──
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    counts_out = output_dir / "gene_count_matrix.tsv"
    meta_out = output_dir / "sample_metadata.tsv"

    merged.to_csv(counts_out, sep='\t')
    merged_meta.to_csv(meta_out, sep='\t')

    print(f"  Saved: {counts_out}")
    print(f"  Saved: {meta_out}")
    print("=" * 60)

    return counts_out, meta_out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Merge count matrices from multiple species directories",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--dirs", nargs='+', required=True,
                        help="Paths to count table directories (each must have "
                             "gene_count_matrix.tsv)")
    parser.add_argument("--output-dir", required=True,
                        help="Directory to write merged files to")

    args = parser.parse_args()
    merge_counts(args.dirs, args.output_dir)
