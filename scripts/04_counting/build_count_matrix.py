#!/usr/bin/env python3
"""
Build gene count matrix from STAR or featureCounts output.

Supports two counting paths:
  PATH A — STAR GeneCounts (your original pipeline)
    STAR --quantMode GeneCounts → one ReadsPerGene.out.tab per sample → merged here
  PATH B — featureCounts (previous student's pipeline)
    samtools sort → featureCounts → single featurecounts.txt → parsed here

Both paths produce the same outputs:
  gene_count_matrix.tsv  — genes as rows, samples as columns (DESeq2 input)
  sample_metadata.tsv    — biological info per sample (DESeq2 design)
  count_summary.txt      — quick sanity-check stats
"""

import os
import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import re


# =============================================================================
# HELPER: EXTRACT SAMPLE NAME FROM FILENAMES
# =============================================================================

def parse_sample_name(filename, file_type='star'):
    """
    Strip path and suffix to get the sample name.

    Examples:
      star:          DC1L1_ReadsPerGene.out.tab          → DC1L1
      featurecounts: /path/DC1L1_Aligned.sortedBySamtools.bam → DC1L1
    """
    basename = os.path.basename(filename)

    if file_type == 'star':
        return basename.replace('_ReadsPerGene.out.tab', '')

    if file_type == 'featurecounts':
        for suffix in ['_Aligned.sortedByCoord.out.bam',
                       '_Aligned.sortedBySamtools.bam',
                       '_Aligned.out.bam',
                       '.sorted.bam',
                       '.bam']:
            if basename.endswith(suffix):
                return basename[:-len(suffix)]
        return basename

    return basename


def extract_sample_info(sample_name):
    """
    Break a sample name like 'DC1L1' into biological metadata.

    Naming convention:  DC  1  L  1
                        ^^  ^  ^  ^
                        |   |  |  └─ replicate number
                        |   |  └──── condition: L = leaf, R = root
                        |   └─────── group number (subspecies/line)
                        └─────────── species code (DC, DG, MF …)
    """
    pattern = r'([A-Z]+)(\d*)([LR])(\d+)'
    match = re.match(pattern, sample_name)

    if match:
        group     = match.group(1)
        group_num = match.group(2) if match.group(2) else ""
        condition = match.group(3)
        replicate = match.group(4)
        variety   = f"{group}{group_num}"
        return {
            'sample': sample_name,
            'group': group,
            'group_number': group_num,
            'condition': condition,
            'replicate': replicate,
            'variety': variety,
            'treatment': variety,
            'full_condition': f"{variety}_{condition}",
        }

    return {
        'sample': sample_name,
        'group': sample_name[:2] if len(sample_name) >= 2 else sample_name,
        'group_number': "",
        'condition': "",
        'replicate': "",
        'variety': sample_name,
        'treatment': sample_name,
        'full_condition': sample_name,
    }


def parse_condition_map(condition_map_str):
    """
    Turn a condition map string into a dictionary

    What it does:
      Takes a string like "_L_=L _R_=R" (from config.sh) and turns it into:
        {"_L_": "L", "_R_": "R"}

      This means: "if a sample name contains '_L_', its condition is L"

    How it works:
      1. Split the string by spaces → ["_L_=L", "_R_=R"]
      2. Split each piece by '='   → ("_L_", "L") and ("_R_", "R")
      3. Store in a dictionary     → {"_L_": "L", "_R_": "R"}
    """
    # This dictionary will hold our substring-to-condition mapping
    mapping = {}

    # Split the string by spaces to get each "substring=condition" pair
    # Example: "_L_=L _R_=R" → ["_L_=L", "_R_=R"]
    pairs = condition_map_str.strip().split()

    for pair in pairs:
        # Each pair should have an '=' sign. If not, skip it and warn.
        if '=' not in pair:
            print(f"  WARNING: Skipping '{pair}' — expected format: substring=condition")
            print(f"  Example: _L_=L  or  Ctrl=C")
            continue

        # Split on the '=' to get the two parts
        # Example: "_L_=L" → substring="_L_", condition="L"
        substring, condition = pair.split('=', 1)

        # Add to our dictionary
        mapping[substring] = condition

    return mapping


def extract_sample_info_from_map(sample_names, condition_map_str):
    """
    Figure out which condition each sample belongs to, using substring matching.

    WHY this exists:
      The default extract_sample_info() uses a regex that only works for names
      like "DC1L1". For non-standard names like "T_L_R1", it fails and produces
      empty metadata. This function uses a simpler approach: check if the sample
      name CONTAINS a known substring (like "_L_") to assign the condition.

    Example:
      sample_names = ["T_L_R1", "T_L_R2", "T_R_R1", "T_R_R2"]
      condition_map_str = "_L_=L _R_=R"

      Result:
        T_L_R1 → contains "_L_" → condition = L, replicate = 1
        T_L_R2 → contains "_L_" → condition = L, replicate = 2
        T_R_R1 → contains "_R_" → condition = R, replicate = 1
        T_R_R2 → contains "_R_" → condition = R, replicate = 2
    """

    # ── Step 1: Parse the condition map string into a dictionary ──
    # "_L_=L _R_=R" becomes {"_L_": "L", "_R_": "R"}
    mapping = parse_condition_map(condition_map_str)

    # ── Step 2: For each sample, find which condition it belongs to ──
    # We check if the sample name contains any of the substrings.
    # Store results in a dictionary: {"T_L_R1": "L", "T_R_R1": "R", ...}
    sample_to_condition = {}

    for name in sample_names:
        # Start with no match
        found_condition = None

        # Check each substring from the mapping
        for substring, condition in mapping.items():
            # Does this sample name contain this substring?
            if substring in name:
                found_condition = condition

        # If no substring matched, something is wrong
        if found_condition is None:
            print(f"  ERROR: Sample '{name}' does not match any condition!")
            print(f"  The condition map is: {condition_map_str}")
            print(f"  Expected the sample name to contain one of: "
                  f"{list(mapping.keys())}")
            found_condition = "UNKNOWN"

        # Save the condition for this sample
        sample_to_condition[name] = found_condition

    # ── Step 3: Assign replicate numbers (1, 2, 3...) within each condition ──
    # Example: T_L_R1 is the 1st Leaf sample, T_L_R2 is the 2nd Leaf sample
    replicate_counter = {}  # keeps track of how many we've seen per condition
    results = []            # the final list of metadata dictionaries

    for name in sample_names:
        condition = sample_to_condition[name]

        # Count up: first sample in this condition = 1, second = 2, etc.
        if condition not in replicate_counter:
            replicate_counter[condition] = 0
        replicate_counter[condition] = replicate_counter[condition] + 1
        replicate_number = replicate_counter[condition]

        # Extract the group name (the part before the condition substring)
        # For "T_L_R1" with substring "_L_": everything before "_L_" is "T"
        group = name.split('_')[0]

        # Build the metadata dictionary (same format as extract_sample_info)
        sample_info = {
            'sample': name,               # full sample name: T_L_R1
            'group': group,                # group/species: T
            'group_number': "",            # not used for these names
            'condition': condition,         # L or R (from the condition map)
            'replicate': str(replicate_number),  # 1, 2, 3...
            'variety': group,                               # T (same as group for non-standard names)
            'treatment': group + "_" + condition,           # T_L or T_R
            'full_condition': group + "_" + condition,      # T_L or T_R
        }
        results.append(sample_info)

    # ── Step 4: Print a summary so the user can verify it looks right ──
    print(f"  Condition map applied: {condition_map_str}")
    for condition in sorted(set(sample_to_condition.values())):
        # Collect all sample names that have this condition
        names_in_condition = []
        for name, cond in sample_to_condition.items():
            if cond == condition:
                names_in_condition.append(name)
        print(f"    {condition}: {names_in_condition}")

    return results


# =============================================================================
# PATH A — READ STAR ReadsPerGene.out.tab (one file per sample)
# =============================================================================

def read_star_counts(filepath):
    """
    Read a single STAR ReadsPerGene.out.tab file.

    Layout (tab-separated, no header):
      Row 1-4: summary stats (N_unmapped, N_multimapping, …) — skipped
      Row 5+:  gene_id   unstranded   sense   antisense

    Returns a DataFrame with columns [gene_id, unstranded].
    """
    try:
        df = pd.read_csv(
            filepath, sep='\t', header=None,
            names=['gene_id', 'unstranded', 'first_strand', 'second_strand'],
        )
        df = df.iloc[4:].reset_index(drop=True)
        df['unstranded'] = pd.to_numeric(df['unstranded'], errors='coerce').fillna(0).astype(int)
        return df[['gene_id', 'unstranded']]
    except Exception as e:
        print(f"  ERROR reading {filepath}: {e}")
        return None


# =============================================================================
# PATH B — READ featureCounts OUTPUT (single combined file, all samples)
# =============================================================================

def read_featurecounts(filepath):
    """
    Parse a featureCounts output file into a gene × sample DataFrame.

    featureCounts format:
      Line 1:  # comment (command used to generate the file)
      Line 2:  Geneid  Chr  Start  End  Strand  Length  <bam1>  <bam2>  …
      Line 3+: gene rows with integer counts per BAM

    The BAM column headers are full paths; we strip them to sample names.
    Returns a DataFrame indexed by Geneid with one column per sample.
    """
    try:
        df = pd.read_csv(filepath, sep='\t', comment='#')

        if 'Geneid' not in df.columns:
            raise ValueError(f"Missing 'Geneid' column — is this a featureCounts file?")

        annotation_cols = {'Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length'}
        bam_cols = [c for c in df.columns if c not in annotation_cols]

        if not bam_cols:
            raise ValueError("No sample (BAM) columns found")

        rename_map = {col: parse_sample_name(col, 'featurecounts') for col in bam_cols}
        counts = df.set_index('Geneid')[bam_cols].rename(columns=rename_map)
        counts = counts.apply(pd.to_numeric, errors='coerce').fillna(0).astype(int)

        print(f"  Parsed featureCounts: {counts.shape[0]} genes × {counts.shape[1]} samples")
        print(f"  Samples: {list(counts.columns)}")
        return counts

    except Exception as e:
        print(f"  ERROR reading {filepath}: {e}")
        return None


# =============================================================================
# MAIN: BUILD THE COUNT MATRIX
# =============================================================================

def build_count_matrix(count_dir, output_dir="count_matrices", count_type="auto",
                       condition_map=None):
    """
    Build gene count matrix from STAR or featureCounts output.

    Args:
        count_dir:     Directory containing count files
        output_dir:    Where to write gene_count_matrix.tsv + metadata
        count_type:    'star', 'featurecounts', or 'auto'
        condition_map: Optional string like "_L_=L _R_=R" for non-standard
                       sample names. When provided, uses substring matching
                       instead of regex to assign conditions.
    """
    count_dir  = Path(count_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    count_matrix_file = output_dir / "gene_count_matrix.tsv"
    metadata_file     = output_dir / "sample_metadata.tsv"

    if count_matrix_file.exists() and metadata_file.exists():
        print(f"\nOutput files already exist in {output_dir}/")
        print(f"  - {count_matrix_file.name}")
        print(f"  - {metadata_file.name}")
        print("Skipping. Delete them or use a different -o directory to regenerate.")
        return (
            pd.read_csv(count_matrix_file, sep='\t', index_col=0),
            pd.read_csv(metadata_file, sep='\t'),
        )

    # ── Auto-detect count type ────────────────────────────────────────────
    if count_type == "auto":
        star_hits = list(count_dir.glob("*ReadsPerGene.out.tab"))
        fc_hits   = [f for f in count_dir.glob("*featurecounts*")
                     if not f.name.endswith('.summary')]
        fc_hits   = sorted(set(fc_hits))

        if star_hits and fc_hits:
            print("Found both STAR and featureCounts files — using featureCounts.")
            count_type = "featurecounts"
        elif fc_hits:
            count_type = "featurecounts"
        elif star_hits:
            count_type = "star"
        else:
            raise ValueError(
                f"No count files found in {count_dir}\n"
                "Expected *ReadsPerGene.out.tab (STAR) or *featurecounts*."
            )

    print(f"Count source: {count_type.upper()}")

    # ── PATH B: featureCounts (single combined file) ──────────────────────
    if count_type == "featurecounts":
        fc_files = [f for f in count_dir.glob("*featurecounts*")
                    if not f.name.endswith('.summary')]
        fc_files = sorted(set(fc_files))

        if not fc_files:
            raise ValueError(f"No featureCounts files found in {count_dir}")
        if len(fc_files) > 1:
            print(f"  Warning: {len(fc_files)} featureCounts files found, using: {fc_files[0].name}")

        count_matrix = read_featurecounts(fc_files[0])
        if count_matrix is None or count_matrix.empty:
            raise ValueError("Failed to parse featureCounts file")

        if condition_map:
            metadata = pd.DataFrame(
                extract_sample_info_from_map(list(count_matrix.columns), condition_map)
            )
        else:
            metadata = pd.DataFrame([extract_sample_info(s) for s in count_matrix.columns])

    # ── PATH A: STAR GeneCounts (one file per sample) ─────────────────────
    elif count_type == "star":
        count_files = list(count_dir.glob("*ReadsPerGene.out.tab"))
        if not count_files:
            raise ValueError(f"No ReadsPerGene.out.tab files in {count_dir}")

        print(f"Found {len(count_files)} STAR count file(s)")

        count_data = {}
        sample_info = []

        for filepath in sorted(count_files):
            sample_name = parse_sample_name(str(filepath), 'star')
            print(f"  Reading {sample_name}...")
            counts = read_star_counts(filepath)
            if counts is not None:
                count_data[sample_name] = counts
                sample_info.append(extract_sample_info(sample_name))

        if not count_data:
            raise ValueError("No valid STAR count files found")

        all_genes = sorted({g for df in count_data.values() for g in df['gene_id']})
        print(f"  {len(all_genes)} genes across {len(count_data)} samples")

        count_matrix = pd.DataFrame(index=all_genes)
        for name, df in count_data.items():
            count_matrix[name] = df.set_index('gene_id')['unstranded']
        count_matrix = count_matrix.fillna(0).astype(int)

        if condition_map:
            metadata = pd.DataFrame(
                extract_sample_info_from_map(list(count_matrix.columns), condition_map)
            )
        else:
            metadata = pd.DataFrame(sample_info)

    else:
        raise ValueError(f"Unknown count_type: {count_type}")

    # ── Save outputs ──────────────────────────────────────────────────────
    print("\nSaving outputs...")

    count_matrix.to_csv(count_matrix_file, sep='\t')
    print(f"  {count_matrix_file}")

    metadata.to_csv(metadata_file, sep='\t', index=False)
    print(f"  {metadata_file}")

    summary_file = output_dir / "count_summary.txt"
    info_rows = metadata.to_dict('records')

    with open(summary_file, 'w') as f:
        f.write("Gene Count Matrix Summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Count source:    {count_type.upper()}\n")
        f.write(f"Total samples:   {count_matrix.shape[1]}\n")
        f.write(f"Total genes:     {count_matrix.shape[0]}\n")
        f.write(f"Matrix shape:    {count_matrix.shape}\n\n")

        f.write("Sample information:\n")
        f.write("-" * 40 + "\n")
        for info in info_rows:
            f.write(f"  {info['sample']}: {info['treatment']} — "
                    f"{info['condition']} (rep {info['replicate']})\n")

        col_sums = count_matrix.sum()
        row_sums = count_matrix.sum(axis=1)
        f.write(f"\nPer-sample count totals:\n")
        f.write("-" * 40 + "\n")
        f.write(f"  Total reads:    {col_sums.sum():,}\n")
        f.write(f"  Mean / sample:  {col_sums.mean():,.0f}\n")
        f.write(f"  Median / sample:{col_sums.median():,.0f}\n")
        f.write(f"  Min:            {col_sums.min():,}\n")
        f.write(f"  Max:            {col_sums.max():,}\n")

        f.write(f"\nGene-level filtering stats:\n")
        f.write("-" * 40 + "\n")
        f.write(f"  Zero counts:     {(row_sums == 0).sum()}\n")
        f.write(f"  < 10 total:      {(row_sums < 10).sum()}\n")
        f.write(f"  < 100 total:     {(row_sums < 100).sum()}\n")

    print(f"  {summary_file}")

    return count_matrix, metadata


# =============================================================================
# CLI ENTRY POINT
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Build gene count matrix from STAR or featureCounts output",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Pipeline paths:

  PATH A — STAR GeneCounts (your original pipeline):
    1. STAR --quantMode GeneCounts  → *_ReadsPerGene.out.tab
    2. python build_count_matrix.py counts/ --type star

  PATH B — featureCounts (previous student's pipeline):
    1. samtools sort                → *_Aligned.sortedBySamtools.bam
    2. featureCounts                → featurecounts.txt
    3. python build_count_matrix.py counts/ --type featurecounts

Examples:
  python build_count_matrix.py 03_count_tables/00_1_DC/ -o 03_count_tables/00_1_DC/ --type star
  python build_count_matrix.py 03_count_tables/00_1_DC/ -o 03_count_tables/00_1_DC/ --type featurecounts
  python build_count_matrix.py 03_count_tables/00_1_DC/   # auto-detect

  # For non-standard sample names (e.g., T_L_R1 instead of DC1L1):
  python build_count_matrix.py counts/ -o counts/ --condition-map "_L_=L _R_=R"
        """,
    )
    parser.add_argument("count_dir", help="Directory containing count files")
    parser.add_argument("-o", "--output", default="count_matrices",
                        help="Output directory (default: count_matrices)")
    parser.add_argument("--type",
                        choices=["star", "featurecounts", "auto"],
                        default="auto",
                        help="Count source: 'star' (ReadsPerGene.out.tab), "
                             "'featurecounts' (featureCounts output), "
                             "or 'auto' (default: auto-detect)")
    parser.add_argument("--condition-map",
                        default=None,
                        help="How to extract conditions from sample names. "
                             "Format: 'substring1=COND1 substring2=COND2'. "
                             "Example: '_L_=L _R_=R' means names containing "
                             "'_L_' are Leaf, '_R_' are Root. "
                             "If not provided, uses regex auto-detection "
                             "(works for DC1L1-style names).")

    args = parser.parse_args()

    try:
        count_matrix, metadata = build_count_matrix(
            args.count_dir, args.output, args.type,
            condition_map=args.condition_map,
        )
        print("\n" + "=" * 60)
        print("SUCCESS — Gene count matrix created!")
        print("=" * 60)
        print(f"Shape:   {count_matrix.shape}")
        print(f"Samples: {list(count_matrix.columns)}")
        print(f"\nNext step — run PyDESeq2:")
        print(f"  python pydeseq2_run_analysis.py "
              f"{args.output}/gene_count_matrix.tsv "
              f"{args.output}/sample_metadata.tsv "
              f"-o pydeseq2_results")

    except Exception as e:
        print(f"\nERROR: {e}")
        return 1

    return 0


if __name__ == "__main__":
    exit(main())
