#!/usr/bin/env python3
"""
Build gene count matrix from STAR ReadsPerGene.out.tab files
and prepare for statistical analysis with DESeq2/edgeR.

WHAT THIS SCRIPT DOES (big picture):
  STAR (the aligner) produces one count file per sample. Each file lists
  every gene and how many reads mapped to it. This script takes all those
  individual files and combines them into ONE table (a "count matrix")
  where rows = genes, columns = samples, and each cell = raw read count.

  That matrix is exactly what DESeq2 needs as input.

INPUT:  A folder full of files like DC1L1_ReadsPerGene.out.tab
OUTPUT: gene_count_matrix.tsv  (the matrix DESeq2 will read)
        sample_metadata.tsv    (which sample belongs to which condition)
        count_summary.txt      (quick stats about the data)
"""

import os
import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import re


# =============================================================================
# STEP 1: HELPER FUNCTIONS — extract info from filenames
# =============================================================================

def parse_sample_name(filename):
    """
    Strip the path and suffix to get just the sample name.
    Example: "/path/to/DC1L1_ReadsPerGene.out.tab" → "DC1L1"
    """
    basename = os.path.basename(filename)
    sample_name = basename.replace('_ReadsPerGene.out.tab', '')
    return sample_name


def extract_sample_info(sample_name):
    """
    Break a sample name like "DC1L1" into its biological meaning using regex.

    Naming convention:  DC  1  L  1
                        ^^  ^  ^  ^
                        |   |  |  └─ replicate number (1, 2, 3...)
                        |   |  └──── condition: L = leaf, R = root
                        |   └─────── group number (subspecies/line)
                        └─────────── species code (DC = D. carota, DG = D. glaber)

    This metadata ends up in sample_metadata.tsv so DESeq2 knows which
    samples to compare against each other.
    """
    # Regex captures each piece: ([A-Z]+) = species, (\d*) = group number,
    # ([LR]) = tissue condition, (\d+) = replicate
    pattern = r'([A-Z]+)(\d*)([LR])(\d+)'
    match = re.match(pattern, sample_name)

    if match:
        group = match.group(1)       # e.g. "DC"
        group_num = match.group(2) if match.group(2) else ""  # e.g. "1"
        condition = match.group(3)   # "L" or "R"
        replicate = match.group(4)   # e.g. "1"

        return {
            'sample': sample_name,             # "DC1L1" — full name
            'group': group,                    # "DC" — species code
            'group_number': group_num,         # "1" — subspecies/line
            'condition': condition,            # "L" — leaf or "R" — root
            'replicate': replicate,            # "1" — biological replicate
            'treatment': f"{group}{group_num}",          # "DC1" — species+line
            'full_condition': f"{group}{group_num}_{condition}"  # "DC1_L"
        }
    else:
        # If the name doesn't match the expected pattern, store what we can
        return {
            'sample': sample_name,
            'group': sample_name[:2] if len(sample_name) >= 2 else sample_name,
            'group_number': "",
            'condition': "",
            'replicate': "",
            'treatment': sample_name,
            'full_condition': sample_name
        }


# =============================================================================
# STEP 2: READ ONE STAR COUNT FILE
# =============================================================================

def read_star_counts(filepath):
    """
    Read a single STAR ReadsPerGene.out.tab file.

    STAR's output looks like this (tab-separated):
        N_unmapped       1234    1234    1234
        N_multimapping   5678    5678    5678
        N_noFeature      9012    9012    9012
        N_ambiguous       345     345     345
        gene_00001        100      50      50    ← actual gene counts start here
        gene_00002        200     150      50
        ...

    First 4 rows = summary stats (unmapped, multimapping, etc.) — we skip these.
    After that, each row is one gene with 4 columns:
        col 1: gene ID
        col 2: unstranded count  ← we use this one
        col 3: sense-strand count
        col 4: antisense-strand count

    We use the UNSTRANDED count (col 2) because our library prep is unstranded.
    """
    try:
        # Read the tab-separated file; there's no header row so header=None
        df = pd.read_csv(filepath, sep='\t', header=None,
                        names=['gene_id', 'unstranded', 'first_strand', 'second_strand'])

        # Skip the first 4 rows (N_unmapped, N_multimapping, N_noFeature, N_ambiguous)
        df = df.iloc[4:].reset_index(drop=True)

        # Make sure counts are integers (some might parse as strings)
        df['unstranded'] = pd.to_numeric(df['unstranded'], errors='coerce').fillna(0).astype(int)

        # Return only gene_id and the unstranded count — that's all we need
        return df[['gene_id', 'unstranded']]
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return None


# =============================================================================
# STEP 3: COMBINE ALL SAMPLES INTO ONE MATRIX
# =============================================================================

def build_count_matrix(count_dir, output_dir="count_matrices"):
    """
    This is the main function. It:
      1. Finds all STAR count files in the input directory
      2. Reads each one (one per sample)
      3. Merges them into a single table: rows = genes, columns = samples
      4. Saves the matrix + metadata + summary stats

    The result is the "count matrix" that DESeq2 expects:
        gene_id     DC1L1   DC1L2   DC1R1   DC1R2   ...
        gene_001      150     130     200     180
        gene_002       42      38      10      12
        ...
    """
    count_dir = Path(count_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    # --- Find all STAR count files in the directory ---
    count_files = list(count_dir.glob("*ReadsPerGene.out.tab"))

    if not count_files:
        raise ValueError(f"No ReadsPerGene.out.tab files found in {count_dir}")

    print(f"Found {len(count_files)} count files")

    # --- Read each file and extract sample metadata ---
    # count_data will be a dict like {"DC1L1": DataFrame, "DC1L2": DataFrame, ...}
    count_data = {}
    sample_info = []

    for filepath in count_files:
        sample_name = parse_sample_name(str(filepath))
        print(f"Processing {sample_name}...")

        counts = read_star_counts(filepath)
        if counts is not None:
            count_data[sample_name] = counts
            sample_info.append(extract_sample_info(sample_name))

    if not count_data:
        raise ValueError("No valid count files found")

    # --- Collect every gene ID across all samples ---
    # Using a set ensures no duplicates; all samples *should* have the same genes,
    # but the set handles any edge cases where a gene appears in one file but not another
    print("Building count matrix...")

    all_genes = set()
    for counts in count_data.values():
        all_genes.update(counts['gene_id'].tolist())

    all_genes = sorted(list(all_genes))
    print(f"Found {len(all_genes)} genes")

    # --- Build the matrix: start with an empty DataFrame with genes as rows ---
    count_matrix = pd.DataFrame(index=all_genes)

    # Add each sample as a new column
    for sample_name, counts in count_data.items():
        counts_indexed = counts.set_index('gene_id')
        count_matrix[sample_name] = counts_indexed['unstranded']

    # If a gene was missing in a sample's file, it gets NaN — fill those with 0
    count_matrix = count_matrix.fillna(0).astype(int)

    # --- Build the metadata table (tells DESeq2 which sample = which condition) ---
    metadata = pd.DataFrame(sample_info)

    # --- Save everything to files ---
    print("Saving outputs...")

    # The count matrix — this is the main input for DESeq2
    count_matrix_file = output_dir / "gene_count_matrix.tsv"
    count_matrix.to_csv(count_matrix_file, sep='\t')
    print(f"Count matrix saved: {count_matrix_file}")

    # The metadata — DESeq2 needs this to know leaf vs root, species, etc.
    metadata_file = output_dir / "sample_metadata.tsv"
    metadata.to_csv(metadata_file, sep='\t', index=False)
    print(f"Sample metadata saved: {metadata_file}")

    # --- Write a human-readable summary so you can sanity-check the data ---
    summary_file = output_dir / "count_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("Gene Count Matrix Summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total samples: {len(count_data)}\n")
        f.write(f"Total genes: {len(all_genes)}\n")
        f.write(f"Matrix dimensions: {count_matrix.shape}\n\n")

        f.write("Sample information:\n")
        f.write("-" * 30 + "\n")
        for info in sample_info:
            f.write(f"{info['sample']}: {info['treatment']} - {info['condition']} (rep {info['replicate']})\n")

        # Per-sample totals — big differences here could indicate a problem
        f.write(f"\nCount statistics:\n")
        f.write("-" * 30 + "\n")
        f.write(f"Total reads: {count_matrix.sum().sum():,}\n")
        f.write(f"Mean reads per sample: {count_matrix.sum().mean():,.0f}\n")
        f.write(f"Median reads per sample: {count_matrix.sum().median():,.0f}\n")
        f.write(f"Min reads per sample: {count_matrix.sum().min():,}\n")
        f.write(f"Max reads per sample: {count_matrix.sum().max():,}\n")

        # Gene-level stats — tells you how many genes have almost no expression,
        # which helps decide filtering thresholds for DESeq2
        f.write(f"\nGene statistics:\n")
        f.write("-" * 30 + "\n")
        f.write(f"Genes with zero counts: {(count_matrix.sum(axis=1) == 0).sum()}\n")
        f.write(f"Genes with < 10 total counts: {(count_matrix.sum(axis=1) < 10).sum()}\n")
        f.write(f"Genes with < 100 total counts: {(count_matrix.sum(axis=1) < 100).sum()}\n")

    print(f"Summary saved: {summary_file}")

    return count_matrix, metadata


# =============================================================================
# STEP 4: COMMAND-LINE ENTRY POINT
# =============================================================================

def main():
    """
    Usage from the terminal:
      python build_count_matrix.py /path/to/count_files -o /path/to/output

    Or via the sbatch wrapper:
      sbatch build_count_matrix.sbatch DC
    """
    parser = argparse.ArgumentParser(description="Build gene count matrix from STAR output")
    parser.add_argument("count_dir", help="Directory containing ReadsPerGene.out.tab files")
    parser.add_argument("-o", "--output", default="count_matrices",
                       help="Output directory (default: count_matrices)")

    args = parser.parse_args()

    try:
        count_matrix, metadata = build_count_matrix(args.count_dir, args.output)
        print("\n" + "="*60)
        print("SUCCESS: Gene count matrix created!")
        print("="*60)
        print(f"Matrix shape: {count_matrix.shape}")
        print(f"Samples: {list(count_matrix.columns)}")
        print(f"\nNext steps:")
        print(f"1. Review the count matrix: {args.output}/gene_count_matrix.tsv")
        print(f"2. Review sample metadata: {args.output}/sample_metadata.tsv")
        print(f"3. Run DESeq2/edgeR analysis using the provided R script")

    except Exception as e:
        print(f"ERROR: {e}")
        return 1

    return 0


if __name__ == "__main__":
    exit(main())
