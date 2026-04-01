#!/usr/bin/env python3
"""
==========================================================================
WHAT THIS SCRIPT DOES:
  Compares YOUR RNA-seq pipeline results against:
    - An official P450 gene list
    - A previous student's results from the same experiment

  It runs 6 checks to see if your results make sense and agree with
  the previous student's findings.

INPUTS (passed as command-line arguments):
  --results          YOUR DE results file (TSV)
  --gene-list        The official P450 gene list (P450_list_RefSeq.txt)
  --prev-de          Previous student's DE gene list with log2fold
  --prev-expression  Previous student's expression data
  --prev-log2fold    Previous student's log2fold data
  --prev-upregulated Previous student's upregulated genes
  --prev-downregulated Previous student's downregulated genes
  -o / --output      Where to save the comparison table

OUTPUTS:
  1. Printed report in the terminal (all 6 checks)
  2. A TSV comparison table (side-by-side data)
  3. A text summary file

BIOLOGY CONTEXT:
  Both you and the previous student compared Root vs Leaf (R vs L).
    - Positive log2FoldChange = gene is MORE active in Root
    - Negative log2FoldChange = gene is MORE active in Leaf
==========================================================================
"""

# -----------------------------------------------------------------------
# STEP 0: Import the tools (libraries) we need
# -----------------------------------------------------------------------

# argparse: lets us read command-line arguments like --results, --gene-list
import argparse

# sys: lets us exit the script with a success/failure code
import sys

# Path: makes it easy to work with file paths
from pathlib import Path

# numpy (np): math library — we use it to check for numeric columns
import numpy as np

# pandas (pd): the main library for working with tables of data
import pandas as pd


# =======================================================================
# HELPER FUNCTIONS: These read different file formats
# =======================================================================
# Why do we need so many?  Because the input files come in different
# formats (.txt, .tsv, .csv) and have different column names.
# These functions figure out how to read each one.
# =======================================================================


def load_gene_ids_from_list_file(filepath):
    """
    Read gene IDs from a gene list file.

    This handles three file types:
      - .txt files: one gene ID per line (like LOC12345)
      - .tsv files: tab-separated with columns
      - .csv files: comma-separated with columns

    Returns a SET of gene IDs (a set is a collection with no duplicates).
    """
    # Convert the filepath to a string so we can check the extension
    path_string = str(filepath)

    # --- CASE 1: Plain text file (.txt) ---
    if path_string.endswith(".txt"):
        gene_ids = set()

        # Open the file and read it line by line
        file_handle = open(filepath)
        for line in file_handle:
            # strip() removes whitespace/newlines from the ends
            cleaned_line = line.strip()

            # Skip empty lines
            if cleaned_line == "":
                continue

            # Skip comment lines (lines starting with #)
            if cleaned_line.startswith("#"):
                continue

            # This is a gene ID — add it to our set
            gene_ids.add(cleaned_line)

        file_handle.close()
        return gene_ids

    # --- CASE 2: TSV or CSV file ---
    # Figure out the separator character
    if path_string.endswith(".tsv") or path_string.endswith(".tab"):
        separator = "\t"   # Tab character for TSV files
    else:
        separator = ","    # Comma for CSV files

    # Read the file into a pandas table (DataFrame)
    data_table = pd.read_csv(filepath, sep=separator)

    # Look for a column named "gene_id"
    if "gene_id" in data_table.columns:
        # Convert the gene_id column to strings and return as a set
        gene_id_column = data_table["gene_id"].astype(str)
        return set(gene_id_column)

    # If no gene_id column, check if the first column has gene IDs
    first_column_name = data_table.columns[0]
    number_of_columns = len(data_table.columns)

    # If the first column starts with "LOC" or there's only 1 column, use it
    if first_column_name.startswith("LOC") or number_of_columns == 1:
        first_column = data_table.iloc[:, 0].astype(str)
        return set(first_column)

    # If we get here, we can't figure out which column has gene IDs
    raise ValueError(
        "Cannot find gene_id column in " + str(filepath)
        + ". Columns found: " + str(list(data_table.columns))
    )


def load_gene_ids_from_results_file(filepath):
    """
    Read gene IDs from YOUR results TSV file.

    Your results file might have gene IDs as:
      - A column named "gene_id"
      - The first column (unnamed)
      - The row index

    Returns a SET of gene IDs.
    """
    # Read the TSV file, ignoring lines that start with #
    data_table = pd.read_csv(filepath, sep="\t", comment="#")

    # Check if there's a column called "gene_id"
    if "gene_id" in data_table.columns:
        return set(data_table["gene_id"].astype(str))

    # Check the first column — it might be unnamed or have LOC-style IDs
    first_column_name = data_table.columns[0]
    first_column_values = data_table.iloc[:, 0].astype(str)

    # Check if any values in the first column start with "LOC"
    has_loc_ids = False
    for value in first_column_values:
        if value.startswith("LOC"):
            has_loc_ids = True
            break

    if first_column_name == "Unnamed: 0" or has_loc_ids:
        return set(first_column_values)

    # Last resort: try reading with the first column as the index
    data_table_v2 = pd.read_csv(filepath, sep="\t", index_col=0, comment="#")
    index_values = data_table_v2.index.astype(str)

    has_loc_in_index = False
    for value in index_values:
        if value.startswith("LOC"):
            has_loc_in_index = True
            break

    if has_loc_in_index:
        return set(index_values)

    raise ValueError(
        "Cannot find gene IDs in " + str(filepath)
        + ". Columns: " + str(list(data_table.columns))
    )


def load_results_as_table(filepath):
    """
    Read YOUR results file as a full pandas table (DataFrame).
    The gene_id becomes the row label (index).

    We need the full table (not just IDs) so we can access columns
    like log2FoldChange, baseMean, padj, etc.
    """
    data_table = pd.read_csv(filepath, sep="\t", comment="#")

    # If there's no gene_id column, figure out which column has the IDs
    if "gene_id" not in data_table.columns:
        first_column_name = data_table.columns[0]
        first_column_values = data_table.iloc[:, 0].astype(str)

        has_loc_ids = False
        for value in first_column_values:
            if value.startswith("LOC"):
                has_loc_ids = True
                break

        if first_column_name == "Unnamed: 0" or has_loc_ids:
            # Rename the first column to "gene_id"
            data_table = data_table.rename(columns={first_column_name: "gene_id"})
        else:
            # Try reading with first column as index
            data_table_v2 = pd.read_csv(filepath, sep="\t", index_col=0, comment="#")
            data_table_v2.index.name = "gene_id"
            data_table = data_table_v2.reset_index()

    # Make sure gene_id values are strings (not numbers)
    data_table["gene_id"] = data_table["gene_id"].astype(str)

    # Set gene_id as the row labels (index)
    data_table = data_table.set_index("gene_id")

    return data_table


def load_previous_student_data(filepath):
    """
    Read the previous student's data file.
    These files can be in many formats, so we have to be flexible.

    Returns a DataFrame with gene_id as the index.
    """
    path_string = str(filepath)

    # Peek at the first line to understand the file format
    file_handle = open(filepath)
    first_line = file_handle.readline()
    file_handle.close()

    # Does this file use tabs to separate columns?
    has_tabs = "\t" in first_line

    # Does the first line look like a header (column names)?
    # If it starts with "LOC", it's probably data, not a header
    first_line_cleaned = first_line.strip()
    looks_like_header = True
    if first_line_cleaned.startswith("LOC"):
        looks_like_header = False

    # ---- Format 1: Tab-separated file with a header row ----
    if has_tabs and looks_like_header:
        data_table = pd.read_csv(filepath, sep="\t", index_col=0)

        if "gene_id" in data_table.columns:
            # Check if the index already has LOC-style IDs
            index_has_loc = False
            for value in data_table.index.astype(str):
                if value.startswith("LOC"):
                    index_has_loc = True
                    break

            if index_has_loc:
                data_table.index = data_table.index.astype(str)
                data_table.index.name = "gene_id"
                # Remove the extra gene_id column since it's already the index
                if "gene_id" in data_table.columns:
                    data_table = data_table.drop(columns=["gene_id"])
            else:
                data_table["gene_id"] = data_table["gene_id"].astype(str)
                data_table = data_table.set_index("gene_id")
        else:
            data_table.index = data_table.index.astype(str)
            data_table.index.name = "gene_id"

        return data_table

    # ---- Format 2: Plain text file (one gene ID per line, maybe with extra values) ----
    if path_string.endswith(".txt") and not has_tabs:
        rows = []

        file_handle = open(filepath)
        for line in file_handle:
            cleaned_line = line.strip()

            # Skip empty lines and comments
            if cleaned_line == "":
                continue
            if cleaned_line.startswith("#"):
                continue

            # Split the line by spaces
            parts = cleaned_line.split()

            # The first part is always the gene ID
            row = {"gene_id": parts[0]}

            # If there are additional values after the gene ID, save them too
            for column_number in range(1, len(parts)):
                value = parts[column_number]
                # Try to convert to a number
                try:
                    row["col_" + str(column_number)] = float(value)
                except ValueError:
                    # If it's not a number, keep it as text
                    row["col_" + str(column_number)] = value

            rows.append(row)
        file_handle.close()

        # If the file was empty, return an empty table
        if len(rows) == 0:
            empty_table = pd.DataFrame(columns=["gene_id"])
            empty_table = empty_table.set_index("gene_id")
            return empty_table

        data_table = pd.DataFrame(rows)
        data_table["gene_id"] = data_table["gene_id"].astype(str)
        data_table = data_table.set_index("gene_id")
        return data_table

    # ---- Format 3: CSV or TSV as a fallback ----
    if path_string.endswith(".tsv") or path_string.endswith(".tab"):
        separator = "\t"
    else:
        separator = ","

    data_table = pd.read_csv(filepath, sep=separator)

    # Try to find which column has the gene IDs
    if "gene_id" not in data_table.columns:
        first_col_name = data_table.columns[0]

        if first_col_name.startswith("LOC"):
            data_table = data_table.rename(columns={first_col_name: "gene_id"})
        elif len(data_table.columns) == 1:
            data_table = data_table.rename(columns={first_col_name: "gene_id"})
        else:
            # Check if first column values look like gene IDs
            first_col_values = data_table.iloc[:, 0].astype(str)
            has_loc = False
            for value in first_col_values:
                if value.startswith("LOC"):
                    has_loc = True
                    break
            if has_loc:
                data_table = data_table.rename(columns={first_col_name: "gene_id"})

    data_table["gene_id"] = data_table["gene_id"].astype(str)
    data_table = data_table.set_index("gene_id")
    return data_table


def find_column_name(data_table, possible_names):
    """
    Search for a column in a table by trying multiple possible names.

    Why? Different files might call the same thing different names:
      - "log2FoldChange" vs "log2fc" vs "LFC"
      - "baseMean" vs "mean_counts" vs "expression"

    This function tries each name (case-insensitive) and returns the
    first match it finds. Returns None if no match.
    """
    # Build a dictionary mapping lowercase column names -> actual column names
    lowercase_to_actual = {}
    for column_name in data_table.columns:
        lowercase_to_actual[column_name.lower()] = column_name

    # Try each candidate name
    for candidate in possible_names:
        if candidate.lower() in lowercase_to_actual:
            return lowercase_to_actual[candidate.lower()]

    # No match found
    return None


# =======================================================================
# CHECK FUNCTIONS: Each one runs a specific verification
# =======================================================================


def check1_gene_list(results_path, list_path):
    """
    CHECK 1: Gene List Membership
    Question: Are all genes in my results also in the official P450 gene list?

    WHY: If a gene shows up in your results but is NOT in the official
    P450 list, something went wrong (maybe a wrong gene snuck in).
    """
    # Load gene IDs from both files
    results_gene_ids = load_gene_ids_from_results_file(results_path)
    list_gene_ids = load_gene_ids_from_list_file(list_path)

    # Find genes that appear in BOTH sets (the overlap)
    in_both = results_gene_ids & list_gene_ids

    # Find genes that are ONLY in results (not in the official list — bad!)
    results_only = results_gene_ids - list_gene_ids

    # Find genes that are ONLY in the list (not in results — filtered out)
    list_only = list_gene_ids - results_gene_ids

    # Print the report
    print("=" * 60)
    print("  CHECK 1: Gene List Membership")
    print("=" * 60)
    print("  Results file: " + str(results_path))
    print("    Genes: " + str(len(results_gene_ids)))
    print("  Gene list:    " + str(list_path))
    print("    Genes: " + str(len(list_gene_ids)))
    print()
    print("  In BOTH (verified):       " + str(len(in_both)))
    print("  In results ONLY:          " + str(len(results_only)))
    print("  In gene list ONLY:        " + str(len(list_only)))
    print()

    # If there are genes in results that shouldn't be there, warn us
    if len(results_only) > 0:
        print("  *** UNEXPECTED: genes in results but NOT in gene list ***")
        for gene in sorted(results_only):
            print("    - " + gene)
        print()

    # Show genes that were in the list but didn't make it into results
    if len(list_only) > 0:
        print("  Genes in gene list but NOT in results (" + str(len(list_only)) + "):")
        print("  (filtered out by DE cutoffs, or not in count matrix)")
        sorted_missing = sorted(list_only)

        if len(sorted_missing) <= 30:
            for gene in sorted_missing:
                print("    - " + gene)
        else:
            # Only show first 15 if there are too many
            for gene in sorted_missing[:15]:
                print("    - " + gene)
            remaining = len(sorted_missing) - 15
            print("    ... and " + str(remaining) + " more")
        print()

    # Final verdict
    if len(results_only) == 0:
        # Calculate what percentage of the gene list made it through
        if len(list_gene_ids) > 0:
            percentage = len(in_both) / len(list_gene_ids) * 100
        else:
            percentage = 0
        print("  PASS: All " + str(len(results_gene_ids)) + " genes in results are in the gene list")
        print("        (" + str(round(percentage, 1)) + "% of the gene list made it through filtering)")
    else:
        print("  FAIL: " + str(len(results_only)) + " gene(s) in results are NOT in the gene list")
        print("        Something unexpected happened — investigate above genes")
    print()

    # Return True if passed (no unexpected genes), False if failed
    did_pass = (len(results_only) == 0)
    return did_pass


def check2_prev_de_overlap(results_path, prev_de_path):
    """
    CHECK 2: Previous Student DE Gene Overlap
    Question: How many of the previous student's DE genes appear in my results too?

    WHY: If both of you found the same genes as differentially expressed,
    that's a good sign — it means the results are reproducible.
    """
    print("=" * 60)
    print("  CHECK 2: Previous Student DE Gene Overlap")
    print("  File: " + str(prev_de_path))
    print("=" * 60)

    # Load both datasets
    previous_data = load_previous_student_data(prev_de_path)
    results_data = load_results_as_table(results_path)

    # Get the gene IDs from each
    previous_gene_ids = set(previous_data.index)
    results_gene_ids = set(results_data.index)

    # Find overlaps
    shared_genes = sorted(previous_gene_ids & results_gene_ids)
    previous_only = sorted(previous_gene_ids - results_gene_ids)
    yours_only = sorted(results_gene_ids - previous_gene_ids)

    print("  Previous student DE genes: " + str(len(previous_gene_ids)))
    print("  Your DE genes:             " + str(len(results_gene_ids)))
    print("  In common:                 " + str(len(shared_genes)))
    print()

    if len(shared_genes) > 0:
        percentage = len(shared_genes) / len(previous_gene_ids) * 100
        print("  " + str(round(percentage)) + "% of previous student's DE genes reproduced in your results")
        print()

    if len(previous_only) > 0:
        print("  In previous student ONLY (" + str(len(previous_only)) + "):")
        print("  (may not pass your DE cutoffs, or not in count matrix)")
        # Show up to 20 genes
        show_count = min(20, len(previous_only))
        for i in range(show_count):
            print("    - " + previous_only[i])
        if len(previous_only) > 20:
            print("    ... and " + str(len(previous_only) - 20) + " more")
        print()

    if len(yours_only) > 0:
        print("  In YOUR results ONLY (" + str(len(yours_only)) + "):")
        print("  (new DE genes not in previous student's list)")
        show_count = min(20, len(yours_only))
        for i in range(show_count):
            print("    - " + yours_only[i])
        if len(yours_only) > 20:
            print("    ... and " + str(len(yours_only) - 20) + " more")
        print()


def check3_expression(results_path, prev_expression_path):
    """
    CHECK 3: Expression Level Comparison
    Question: Do our expression levels (like baseMean) correlate?

    WHY: Even if you don't have the exact same DE genes, the overall
    expression patterns should be similar if both experiments worked.
    A high correlation means "when a gene is highly expressed in their data,
    it's also highly expressed in yours."
    """
    print("=" * 60)
    print("  CHECK 3: Expression Level Comparison")
    print("  File: " + str(prev_expression_path))
    print("=" * 60)

    previous_data = load_previous_student_data(prev_expression_path)
    results_data = load_results_as_table(results_path)

    print("  Previous student genes: " + str(len(previous_data)))
    print("  Previous student columns: " + str(list(previous_data.columns)))
    print("  Your genes:             " + str(len(results_data)))
    print()

    # Find genes in common
    previous_gene_ids = set(previous_data.index)
    results_gene_ids = set(results_data.index)
    shared_genes = sorted(previous_gene_ids & results_gene_ids)

    print("  Genes in common: " + str(len(shared_genes)))
    print()

    if len(shared_genes) == 0:
        print("  No overlapping genes to compare.")
        print()
        return

    # Look for expression/count columns in previous data
    expression_column_candidates = [
        "baseMean", "mean_counts", "total_counts", "expression",
        "counts", "RPKM", "FPKM", "TPM", "count", "mean",
    ]
    prev_expression_column = find_column_name(previous_data, expression_column_candidates)

    your_expression_candidates = ["baseMean", "mean_counts", "total_counts"]
    your_expression_column = find_column_name(results_data, your_expression_candidates)

    # If both have expression columns, calculate correlation
    if prev_expression_column is not None and your_expression_column is not None:
        previous_values = previous_data.loc[shared_genes, prev_expression_column].astype(float)
        your_values = results_data.loc[shared_genes, your_expression_column].astype(float)

        # Only compare genes where both have valid (non-missing) data
        both_have_data = previous_values.notna() & your_values.notna()
        genes_with_data = both_have_data.sum()

        if genes_with_data >= 3:
            # Pearson correlation: 1.0 = perfect match, 0 = no relationship
            correlation = previous_values[both_have_data].corr(your_values[both_have_data])
            print("  Expression correlation (Pearson): " + str(round(correlation, 3)))
            print("    Compared: prev '" + prev_expression_column + "' vs yours '" + your_expression_column + "'")
            print("    Genes with data: " + str(genes_with_data))

            if correlation > 0.7:
                print("    Strong positive correlation — expression patterns agree")
            elif correlation > 0.3:
                print("    Moderate correlation")
            else:
                print("    Weak/no correlation — different expression patterns")
            print()

    # Show the top 10 most highly expressed shared genes
    if your_expression_column is not None:
        shared_results = results_data.loc[shared_genes]
        top_10 = shared_results.nlargest(10, your_expression_column)
        print("  Top 10 shared genes by your " + your_expression_column + ":")

        for gene_id in top_10.index:
            your_value = top_10.at[gene_id, your_expression_column]
            line = "    " + gene_id + "  yours=" + str(round(your_value, 1))

            # Also show the previous student's value if available
            if prev_expression_column is not None:
                if gene_id in previous_data.index:
                    prev_value = previous_data.at[gene_id, prev_expression_column]
                    if pd.notna(prev_value):
                        line = line + "  prev " + prev_expression_column + "=" + str(round(float(prev_value), 1))

            print(line)
        print()

    # If the previous data has many numeric columns, show a summary
    numeric_column_names = []
    for column_name in previous_data.columns:
        if previous_data[column_name].dtype in [np.float64, np.int64, float, int]:
            numeric_column_names.append(column_name)

    if len(numeric_column_names) > 2:
        print("  Previous data has " + str(len(numeric_column_names)) + " numeric columns:")
        show_count = min(10, len(numeric_column_names))
        for i in range(show_count):
            col = numeric_column_names[i]
            mean_value = previous_data[col].mean()
            print("    " + col + ": mean=" + str(round(mean_value, 1)))
        if len(numeric_column_names) > 10:
            print("    ... and " + str(len(numeric_column_names) - 10) + " more")
        print()


def check4_log2fold(results_path, prev_log2fold_path):
    """
    CHECK 4: Log2FoldChange Direction & Magnitude
    Question: When the previous student says a gene is upregulated in Root,
    do I also see it as upregulated in Root?

    WHY: This is the most important check. If most genes point the same
    direction (both say "up in root" or both say "up in leaf"), that means
    you're getting consistent biological results.
    """
    print("=" * 60)
    print("  CHECK 4: Log2FoldChange Direction & Magnitude")
    print("  File: " + str(prev_log2fold_path))
    print("=" * 60)

    previous_data = load_previous_student_data(prev_log2fold_path)
    results_data = load_results_as_table(results_path)

    print("  Previous student genes: " + str(len(previous_data)))
    print("  Previous student columns: " + str(list(previous_data.columns)))
    print("  Your genes:             " + str(len(results_data)))
    print()

    previous_gene_ids = set(previous_data.index)
    results_gene_ids = set(results_data.index)
    shared_genes = sorted(previous_gene_ids & results_gene_ids)

    print("  Genes in common: " + str(len(shared_genes)))
    print()

    if len(shared_genes) == 0:
        print("  No overlapping genes to compare.")
        print()
        return

    # Find the log2FoldChange column in the previous student's data
    prev_lfc_candidates = [
        "log2FoldChange", "log2fold", "log2fc", "lfc", "FC",
        "log2_fold_change", "foldchange", "fold_change",
    ]
    prev_lfc_column = find_column_name(previous_data, prev_lfc_candidates)

    # Find the log2FoldChange column in your results
    your_lfc_candidates = ["log2FoldChange", "log2fold", "log2fc"]
    your_lfc_column = find_column_name(results_data, your_lfc_candidates)

    # If we can't find the column by name, but there's only 1 numeric column,
    # assume that's the log2FC
    if prev_lfc_column is None:
        numeric_columns = []
        for col in previous_data.columns:
            if previous_data[col].dtype in [np.float64, np.int64, float, int]:
                numeric_columns.append(col)
        if len(numeric_columns) == 1:
            prev_lfc_column = numeric_columns[0]
            print("  (Assuming '" + prev_lfc_column + "' is the log2FC column)")
            print()

    if prev_lfc_column is None:
        print("  Cannot find log2FoldChange column in previous data.")
        print("  Available columns: " + str(list(previous_data.columns)))
        print()
        report_overlap_only(shared_genes, previous_gene_ids, results_gene_ids)
        return

    if your_lfc_column is None:
        print("  Cannot find log2FoldChange column in your results.")
        print("  Available columns: " + str(list(results_data.columns)))
        print()
        report_overlap_only(shared_genes, previous_gene_ids, results_gene_ids)
        return

    # ---- Compare direction for each shared gene ----
    agree_count = 0
    disagree_count = 0
    no_data_count = 0
    disagree_genes = []     # Genes where you and prev student disagree
    agree_genes = []        # Genes where you agree

    for gene in shared_genes:
        # Get the previous student's log2FC for this gene
        prev_lfc_value = previous_data.at[gene, prev_lfc_column]
        your_lfc_value = results_data.at[gene, your_lfc_column]

        # Skip if either value is missing (NaN)
        if pd.isna(prev_lfc_value) or pd.isna(your_lfc_value):
            no_data_count = no_data_count + 1
            continue

        # Convert to regular Python numbers
        prev_lfc_value = float(prev_lfc_value)
        your_lfc_value = float(your_lfc_value)

        # Determine direction: positive = up in root, negative = up in leaf
        if prev_lfc_value > 0:
            prev_direction = "up_root"
        else:
            prev_direction = "up_leaf"

        if your_lfc_value > 0:
            your_direction = "up_root"
        else:
            your_direction = "up_leaf"

        # Do they agree?
        if prev_direction == your_direction:
            agree_count = agree_count + 1
            agree_genes.append((gene, prev_lfc_value, your_lfc_value, your_direction))
        else:
            disagree_count = disagree_count + 1
            disagree_genes.append((gene, prev_lfc_value, your_lfc_value, prev_direction, your_direction))

    total_compared = agree_count + disagree_count
    if total_compared > 0:
        agree_percentage = agree_count / total_compared * 100
    else:
        agree_percentage = 0

    print("  Direction comparison (" + str(total_compared) + " genes with log2FC in both):")
    print("    AGREE  (same direction):   " + str(agree_count) + "  (" + str(round(agree_percentage)) + "%)")
    print("    DISAGREE (opposite):       " + str(disagree_count))
    if no_data_count > 0:
        print("    Missing data:              " + str(no_data_count))
    print()

    if agree_percentage >= 80:
        print("  Strong agreement — " + str(round(agree_percentage)) + "% of shared genes have same direction")
    elif agree_percentage >= 50:
        print("  Moderate agreement — " + str(round(agree_percentage)) + "% match, review disagreements")
    else:
        print("  Low agreement — only " + str(round(agree_percentage)) + "% match, check experimental conditions")
    print()

    # ---- Calculate magnitude correlation ----
    prev_values_list = []
    your_values_list = []

    for gene in shared_genes:
        prev_val = previous_data.at[gene, prev_lfc_column]
        your_val = results_data.at[gene, your_lfc_column]
        if pd.notna(prev_val) and pd.notna(your_val):
            prev_values_list.append(float(prev_val))
            your_values_list.append(float(your_val))

    if len(prev_values_list) >= 3:
        prev_series = pd.Series(prev_values_list)
        your_series = pd.Series(your_values_list)
        correlation = prev_series.corr(your_series)
        print("  Log2FC magnitude correlation (Pearson): " + str(round(correlation, 3)))
        if correlation > 0.7:
            print("    Strong — similar fold-change magnitudes")
        elif correlation > 0.3:
            print("    Moderate — some consistency in magnitudes")
        else:
            print("    Weak — fold-change magnitudes differ")
        print()

    # ---- Show genes that DISAGREE ----
    if len(disagree_genes) > 0:
        print("  Genes with OPPOSITE direction (" + str(len(disagree_genes)) + "):")
        show_count = min(20, len(disagree_genes))
        for i in range(show_count):
            gene = disagree_genes[i][0]
            prev_lfc = disagree_genes[i][1]
            your_lfc = disagree_genes[i][2]
            prev_dir = disagree_genes[i][3]
            your_dir = disagree_genes[i][4]
            print("    " + gene
                  + "  prev=" + format(prev_lfc, "+.2f") + " (" + prev_dir + ")"
                  + "  yours=" + format(your_lfc, "+.2f") + " (" + your_dir + ")")
        if len(disagree_genes) > 20:
            print("    ... and " + str(len(disagree_genes) - 20) + " more")
        print()

    # ---- Show top genes that AGREE ----
    if len(agree_genes) > 0:
        # Sort by the absolute value of YOUR log2FC (biggest changes first)
        agree_genes.sort(key=lambda x: abs(x[2]), reverse=True)
        top_agree = agree_genes[:10]
        print("  Top 10 AGREEING genes (by your |log2FC|):")
        for i in range(len(top_agree)):
            gene = top_agree[i][0]
            prev_lfc = top_agree[i][1]
            your_lfc = top_agree[i][2]
            direction = top_agree[i][3]
            print("    " + gene
                  + "  prev=" + format(prev_lfc, "+.2f")
                  + "  yours=" + format(your_lfc, "+.2f")
                  + "  (" + direction + ")")
        print()

    # ---- Compare statistical significance (padj) ----
    prev_padj_column = find_column_name(previous_data, ["padj", "p_adj", "adjusted_pvalue", "fdr", "qvalue"])
    your_padj_column = find_column_name(results_data, ["padj"])

    if prev_padj_column is not None and your_padj_column is not None:
        both_significant = 0
        prev_significant_only = 0
        your_significant_only = 0

        for gene in shared_genes:
            prev_pval = previous_data.at[gene, prev_padj_column]
            your_pval = results_data.at[gene, your_padj_column]

            # Treat missing p-values as not significant (1.0)
            if pd.isna(prev_pval):
                prev_pval = 1.0
            if pd.isna(your_pval):
                your_pval = 1.0

            prev_is_significant = float(prev_pval) < 0.05
            your_is_significant = float(your_pval) < 0.05

            if prev_is_significant and your_is_significant:
                both_significant = both_significant + 1
            elif prev_is_significant:
                prev_significant_only = prev_significant_only + 1
            elif your_is_significant:
                your_significant_only = your_significant_only + 1

        print("  Significance comparison (padj < 0.05):")
        print("    Both significant:          " + str(both_significant))
        print("    Previous only:             " + str(prev_significant_only))
        print("    Yours only:                " + str(your_significant_only))
        print()


def check5_upregulated(results_path, upregulated_path):
    """
    CHECK 5: Upregulated Gene Agreement
    Question: Genes the previous student found upregulated in Root —
    are they also upregulated in Root in my data?

    WHY: If both studies say the same genes are "turned on more in Root",
    that confirms the biology is real and reproducible.
    """
    print("=" * 60)
    print("  CHECK 5: Upregulated Gene Agreement")
    print("  File: " + str(upregulated_path))
    print("=" * 60)

    # Read the upregulated gene list (plain text, one gene per line)
    raw_gene_ids = []
    file_handle = open(upregulated_path)
    for line in file_handle:
        cleaned = line.strip()
        if cleaned == "":
            continue
        if cleaned.startswith("#"):
            continue
        raw_gene_ids.append(cleaned)
    file_handle.close()

    # Some gene IDs might be just numbers — add "LOC" prefix if needed
    previous_upregulated = set()
    for gene_id in raw_gene_ids:
        if gene_id.startswith("LOC"):
            previous_upregulated.add(gene_id)
        elif gene_id.isdigit():
            previous_upregulated.add("LOC" + gene_id)
        else:
            previous_upregulated.add(gene_id)

    results_data = load_results_as_table(results_path)
    results_gene_ids = set(results_data.index)
    shared_genes = sorted(previous_upregulated & results_gene_ids)
    previous_only = sorted(previous_upregulated - results_gene_ids)

    print("  Previous upregulated genes: " + str(len(previous_upregulated)))
    print("  Your DE genes:              " + str(len(results_gene_ids)))
    print("  In common:                  " + str(len(shared_genes)))
    print()

    if len(shared_genes) == 0:
        print("  No overlapping genes to compare.")
        if len(previous_only) > 0:
            print("  All " + str(len(previous_only)) + " previous upregulated genes are missing from your results")
            print("  (may not pass your DE cutoffs)")
        print()
        return

    lfc_column = find_column_name(results_data, ["log2FoldChange", "log2fold", "log2fc"])
    if lfc_column is None:
        print("  Cannot find log2FoldChange in your results.")
        print()
        return

    # Check each shared gene: is it also upregulated (positive log2FC) in your data?
    confirmed_up = []        # Genes that ARE upregulated in your data too
    flipped_to_down = []     # Genes that are downregulated in your data (opposite!)
    no_data_count = 0

    for gene in shared_genes:
        lfc_value = results_data.at[gene, lfc_column]

        if pd.isna(lfc_value):
            no_data_count = no_data_count + 1
            continue

        lfc_value = float(lfc_value)

        if lfc_value > 0:
            # Positive = upregulated in Root — matches!
            confirmed_up.append((gene, lfc_value))
        else:
            # Negative = downregulated in Root — opposite!
            flipped_to_down.append((gene, lfc_value))

    total_checked = len(confirmed_up) + len(flipped_to_down)
    if total_checked > 0:
        confirmed_percentage = len(confirmed_up) / total_checked * 100
    else:
        confirmed_percentage = 0

    print("  Direction check (" + str(total_checked) + " genes with log2FC):")
    print("    CONFIRMED upregulated:     " + str(len(confirmed_up)) + "  (" + str(round(confirmed_percentage)) + "%)")
    print("    FLIPPED to downregulated:  " + str(len(flipped_to_down)))
    if no_data_count > 0:
        print("    Missing log2FC:            " + str(no_data_count))
    print()

    if confirmed_percentage >= 80:
        print("  Strong agreement — " + str(round(confirmed_percentage)) + "% confirmed as upregulated")
    elif confirmed_percentage >= 50:
        print("  Moderate — " + str(round(confirmed_percentage)) + "% confirmed, review flipped genes")
    else:
        print("  Low agreement — only " + str(round(confirmed_percentage)) + "% confirmed, check contrast direction")
        print("  (your contrast may be R vs L while theirs was L vs R)")
    print()

    # Show top confirmed upregulated genes
    if len(confirmed_up) > 0:
        confirmed_up.sort(key=lambda x: x[1], reverse=True)
        show_count = min(10, len(confirmed_up))
        print("  Top confirmed upregulated (by your log2FC):")
        for i in range(show_count):
            gene = confirmed_up[i][0]
            lfc = confirmed_up[i][1]
            print("    " + gene + "  log2FC=" + format(lfc, "+.2f"))
        print()

    # Show flipped genes
    if len(flipped_to_down) > 0:
        flipped_to_down.sort(key=lambda x: x[1])
        print("  FLIPPED genes (prev=up, yours=down):")
        for i in range(len(flipped_to_down)):
            gene = flipped_to_down[i][0]
            lfc = flipped_to_down[i][1]
            print("    " + gene + "  log2FC=" + format(lfc, "+.2f"))
        print()

    # Show genes that were in previous but not in your results
    if len(previous_only) > 0:
        print("  Previous upregulated NOT in your results (" + str(len(previous_only)) + "):")
        print("  (may not pass your DE cutoffs or not in count matrix)")
        show_count = min(20, len(previous_only))
        for i in range(show_count):
            print("    - " + previous_only[i])
        if len(previous_only) > 20:
            print("    ... and " + str(len(previous_only) - 20) + " more")
        print()


def check6_downregulated(results_path, downregulated_path):
    """
    CHECK 6: Downregulated Gene Agreement
    Question: Genes the previous student found downregulated in Root —
    are they also downregulated in Root in my data?

    WHY: Same logic as Check 5, but for downregulated genes.
    Negative log2FC = higher in Leaf = downregulated in Root.
    """
    print("=" * 60)
    print("  CHECK 6: Downregulated Gene Agreement")
    print("  File: " + str(downregulated_path))
    print("=" * 60)

    # Read the downregulated gene list
    raw_gene_ids = []
    file_handle = open(downregulated_path)
    for line in file_handle:
        cleaned = line.strip()
        if cleaned == "":
            continue
        if cleaned.startswith("#"):
            continue
        raw_gene_ids.append(cleaned)
    file_handle.close()

    previous_downregulated = set()
    for gene_id in raw_gene_ids:
        if gene_id.startswith("LOC"):
            previous_downregulated.add(gene_id)
        elif gene_id.isdigit():
            previous_downregulated.add("LOC" + gene_id)
        else:
            previous_downregulated.add(gene_id)

    results_data = load_results_as_table(results_path)
    results_gene_ids = set(results_data.index)
    shared_genes = sorted(previous_downregulated & results_gene_ids)
    previous_only = sorted(previous_downregulated - results_gene_ids)

    print("  Previous downregulated genes: " + str(len(previous_downregulated)))
    print("  Your DE genes:                " + str(len(results_gene_ids)))
    print("  In common:                    " + str(len(shared_genes)))
    print()

    if len(shared_genes) == 0:
        print("  No overlapping genes to compare.")
        if len(previous_only) > 0:
            print("  All " + str(len(previous_only)) + " previous downregulated genes are missing from your results")
            print("  (may not pass your DE cutoffs)")
        print()
        return

    lfc_column = find_column_name(results_data, ["log2FoldChange", "log2fold", "log2fc"])
    if lfc_column is None:
        print("  Cannot find log2FoldChange in your results.")
        print()
        return

    confirmed_down = []      # Genes that ARE downregulated in your data too
    flipped_to_up = []       # Genes that are upregulated in your data (opposite!)
    no_data_count = 0

    for gene in shared_genes:
        lfc_value = results_data.at[gene, lfc_column]

        if pd.isna(lfc_value):
            no_data_count = no_data_count + 1
            continue

        lfc_value = float(lfc_value)

        if lfc_value < 0:
            # Negative = downregulated in Root — matches!
            confirmed_down.append((gene, lfc_value))
        else:
            # Positive = upregulated in Root — opposite!
            flipped_to_up.append((gene, lfc_value))

    total_checked = len(confirmed_down) + len(flipped_to_up)
    if total_checked > 0:
        confirmed_percentage = len(confirmed_down) / total_checked * 100
    else:
        confirmed_percentage = 0

    print("  Direction check (" + str(total_checked) + " genes with log2FC):")
    print("    CONFIRMED downregulated:   " + str(len(confirmed_down)) + "  (" + str(round(confirmed_percentage)) + "%)")
    print("    FLIPPED to upregulated:    " + str(len(flipped_to_up)))
    if no_data_count > 0:
        print("    Missing log2FC:            " + str(no_data_count))
    print()

    if confirmed_percentage >= 80:
        print("  Strong agreement — " + str(round(confirmed_percentage)) + "% confirmed as downregulated")
    elif confirmed_percentage >= 50:
        print("  Moderate — " + str(round(confirmed_percentage)) + "% confirmed, review flipped genes")
    else:
        print("  Low agreement — only " + str(round(confirmed_percentage)) + "% confirmed, check contrast direction")
        print("  (your contrast should be R vs L, same as previous student)")
    print()

    if len(confirmed_down) > 0:
        confirmed_down.sort(key=lambda x: x[1])
        show_count = min(10, len(confirmed_down))
        print("  Top confirmed downregulated (most negative log2FC):")
        for i in range(show_count):
            gene = confirmed_down[i][0]
            lfc = confirmed_down[i][1]
            print("    " + gene + "  log2FC=" + format(lfc, "+.2f"))
        print()

    if len(flipped_to_up) > 0:
        flipped_to_up.sort(key=lambda x: x[1], reverse=True)
        print("  FLIPPED genes (prev=down, yours=up):")
        for i in range(len(flipped_to_up)):
            gene = flipped_to_up[i][0]
            lfc = flipped_to_up[i][1]
            print("    " + gene + "  log2FC=" + format(lfc, "+.2f"))
        print()

    if len(previous_only) > 0:
        print("  Previous downregulated NOT in your results (" + str(len(previous_only)) + "):")
        print("  (may not pass your DE cutoffs or not in count matrix)")
        show_count = min(20, len(previous_only))
        for i in range(show_count):
            print("    - " + previous_only[i])
        if len(previous_only) > 20:
            print("    ... and " + str(len(previous_only) - 20) + " more")
        print()


def report_overlap_only(shared_genes, previous_gene_ids, results_gene_ids):
    """
    Fallback function: when we can't compare column values,
    just report which genes overlap between the two datasets.
    """
    previous_only = sorted(previous_gene_ids - results_gene_ids)
    yours_only = sorted(results_gene_ids - previous_gene_ids)

    if len(shared_genes) > 0:
        print("  Shared genes (" + str(len(shared_genes)) + "):")
        show_count = min(30, len(shared_genes))
        for i in range(show_count):
            print("    - " + shared_genes[i])
        if len(shared_genes) > 30:
            print("    ... and " + str(len(shared_genes) - 30) + " more")
        print()

    if len(previous_only) > 0:
        print("  In previous ONLY (" + str(len(previous_only)) + "):")
        show_count = min(15, len(previous_only))
        for i in range(show_count):
            print("    - " + previous_only[i])
        if len(previous_only) > 15:
            print("    ... and " + str(len(previous_only) - 15) + " more")
        print()

    if len(yours_only) > 0:
        print("  In yours ONLY (" + str(len(yours_only)) + "):")
        show_count = min(15, len(yours_only))
        for i in range(show_count):
            print("    - " + yours_only[i])
        if len(yours_only) > 15:
            print("    ... and " + str(len(yours_only) - 15) + " more")
        print()


# =======================================================================
# COMPARISON TABLE: Build a big side-by-side table of all genes
# =======================================================================


def build_comparison_table(results_path, list_path,
                           prev_de_path=None, prev_expr_path=None,
                           prev_lfc_path=None, prev_upreg_path=None,
                           prev_downreg_path=None):
    """
    Build one giant table where every row is a gene and columns show
    your data vs the previous student's data side by side.
    """
    results_data = load_results_as_table(results_path)
    gene_list_ids = load_gene_ids_from_list_file(list_path)

    # Combine all gene IDs from your results AND the gene list
    all_gene_ids = set(gene_list_ids)
    for gene_id in results_data.index:
        all_gene_ids.add(gene_id)
    all_gene_ids = sorted(all_gene_ids)

    # Build the table row by row
    rows = []
    for gene_id in all_gene_ids:
        row = {"gene_id": gene_id}

        # Is this gene in your results?
        is_in_results = gene_id in results_data.index
        if is_in_results:
            row["in_your_results"] = "yes"
        else:
            row["in_your_results"] = "no"

        # Is this gene in the official gene list?
        if gene_id in gene_list_ids:
            row["in_gene_list"] = "yes"
        else:
            row["in_gene_list"] = "no"

        # If the gene is in your results, copy over the important columns
        if is_in_results:
            important_columns = ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
            for col in important_columns:
                if col in results_data.columns:
                    value = results_data.at[gene_id, col]
                    if pd.notna(value):
                        row["your_" + col] = value
                    else:
                        row["your_" + col] = ""

            # Determine direction (up in root or up in leaf)
            if "log2FoldChange" in results_data.columns:
                lfc = results_data.at[gene_id, "log2FoldChange"]
                if pd.notna(lfc):
                    if float(lfc) > 0:
                        row["your_direction"] = "up_root"
                    else:
                        row["your_direction"] = "up_leaf"
                else:
                    row["your_direction"] = ""

            if "mean_counts" in results_data.columns:
                row["your_mean_counts"] = results_data.at[gene_id, "mean_counts"]
            if "total_counts" in results_data.columns:
                row["your_total_counts"] = results_data.at[gene_id, "total_counts"]

        rows.append(row)

    table = pd.DataFrame(rows)
    table = table.set_index("gene_id")

    # ---- Add previous student DE list info (Check 2) ----
    if prev_de_path is not None and Path(prev_de_path).exists():
        prev_de_data = load_previous_student_data(prev_de_path)
        column_values = []
        for gene_id in table.index:
            if gene_id in prev_de_data.index:
                column_values.append("yes")
            else:
                column_values.append("no")
        table["in_prev_DE_list"] = column_values

    # ---- Add previous student expression info (Check 3) ----
    if prev_expr_path is not None and Path(prev_expr_path).exists():
        prev_expr_data = load_previous_student_data(prev_expr_path)

        column_values = []
        for gene_id in table.index:
            if gene_id in prev_expr_data.index:
                column_values.append("yes")
            else:
                column_values.append("no")
        table["in_prev_expression"] = column_values

        for col in ["baseMean", "log2FoldChange", "padj"]:
            if col in prev_expr_data.columns:
                column_values = []
                for gene_id in table.index:
                    if gene_id in prev_expr_data.index and pd.notna(prev_expr_data.at[gene_id, col]):
                        column_values.append(prev_expr_data.at[gene_id, col])
                    else:
                        column_values.append("")
                table["prev_expr_" + col] = column_values

    # ---- Add previous student log2fold info (Check 4) ----
    if prev_lfc_path is not None and Path(prev_lfc_path).exists():
        prev_lfc_data = load_previous_student_data(prev_lfc_path)

        column_values = []
        for gene_id in table.index:
            if gene_id in prev_lfc_data.index:
                column_values.append("yes")
            else:
                column_values.append("no")
        table["in_prev_log2fold"] = column_values

        for col in ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]:
            if col in prev_lfc_data.columns:
                column_values = []
                for gene_id in table.index:
                    if gene_id in prev_lfc_data.index and pd.notna(prev_lfc_data.at[gene_id, col]):
                        column_values.append(prev_lfc_data.at[gene_id, col])
                    else:
                        column_values.append("")
                table["prev_lfc_" + col] = column_values

        lfc_col = find_column_name(prev_lfc_data, ["log2FoldChange", "log2fold", "log2fc"])
        if lfc_col is not None:
            direction_values = []
            for gene_id in table.index:
                if gene_id in prev_lfc_data.index and pd.notna(prev_lfc_data.at[gene_id, lfc_col]):
                    if float(prev_lfc_data.at[gene_id, lfc_col]) > 0:
                        direction_values.append("up_root")
                    else:
                        direction_values.append("up_leaf")
                else:
                    direction_values.append("")
            table["prev_lfc_direction"] = direction_values

    # ---- Add previous upregulated list (Check 5) ----
    if prev_upreg_path is not None and Path(prev_upreg_path).exists():
        raw_ids = []
        fh = open(prev_upreg_path)
        for line in fh:
            cleaned = line.strip()
            if cleaned != "" and not cleaned.startswith("#"):
                raw_ids.append(cleaned)
        fh.close()

        upreg_ids = set()
        for gid in raw_ids:
            if gid.isdigit():
                upreg_ids.add("LOC" + gid)
            else:
                upreg_ids.add(gid)

        column_values = []
        for gene_id in table.index:
            if gene_id in upreg_ids:
                column_values.append("yes")
            else:
                column_values.append("no")
        table["prev_upregulated"] = column_values

    # ---- Add previous downregulated list (Check 6) ----
    if prev_downreg_path is not None and Path(prev_downreg_path).exists():
        raw_ids = []
        fh = open(prev_downreg_path)
        for line in fh:
            cleaned = line.strip()
            if cleaned != "" and not cleaned.startswith("#"):
                raw_ids.append(cleaned)
        fh.close()

        downreg_ids = set()
        for gid in raw_ids:
            if gid.isdigit():
                downreg_ids.add("LOC" + gid)
            else:
                downreg_ids.add(gid)

        column_values = []
        for gene_id in table.index:
            if gene_id in downreg_ids:
                column_values.append("yes")
            else:
                column_values.append("no")
        table["prev_downregulated"] = column_values

    # ---- Add direction agreement column ----
    if "your_direction" in table.columns and "prev_lfc_direction" in table.columns:
        agreement_values = []
        for gene_id in table.index:
            your_dir = table.at[gene_id, "your_direction"]
            prev_dir = table.at[gene_id, "prev_lfc_direction"]

            if your_dir == "" or prev_dir == "":
                agreement_values.append("")
            elif your_dir == prev_dir:
                agreement_values.append("agree")
            else:
                agreement_values.append("DISAGREE")
        table["direction_agreement"] = agreement_values

    # ---- Sort the table ----
    # Genes in your results first, then by padj (most significant first)
    if "your_padj" in table.columns:
        # Create a temporary column for sorting
        sort_padj = pd.to_numeric(table["your_padj"], errors="coerce").fillna(1)
        table["_sort_padj"] = sort_padj
        table = table.sort_values(
            ["in_your_results", "_sort_padj"],
            ascending=[False, True]
        )
        table = table.drop(columns=["_sort_padj"])
    elif "in_your_results" in table.columns:
        table = table.sort_values("in_your_results", ascending=False)

    return table


# =======================================================================
# SUMMARY WRITER: Create a human-readable summary text file
# =======================================================================


def write_summary(table, summary_path, args):
    """
    Write a nicely formatted summary file that shows the results
    of all 6 checks in one place — like a final report card.
    """
    from datetime import datetime

    # We'll collect all lines in a list, then write them all at once
    lines = []

    # Count some basics
    total_genes = len(table)
    in_yours_mask = table["in_your_results"] == "yes"
    genes_in_yours = in_yours_mask.sum()
    in_list_mask = table["in_gene_list"] == "yes"
    genes_in_list = in_list_mask.sum()

    # ---- Header ----
    # CHANGED: 2026-04-01 — Use --gene-family arg in header (was hardcoded "P450")
    family_label = "CYP"
    if hasattr(args, "gene_family") and args.gene_family:
        family_label = args.gene_family.upper()

    lines.append("=" * 70)
    lines.append("  " + family_label + " VERIFICATION SUMMARY")
    lines.append("  Generated: " + datetime.now().strftime("%Y-%m-%d %H:%M"))
    lines.append("=" * 70)
    lines.append("")

    # ---- Overview ----
    lines.append("OVERVIEW")
    lines.append("-" * 70)
    lines.append("  Your results file:          " + str(args.results))
    lines.append("  Original gene list:         " + str(args.gene_list))
    lines.append("  Gene family:                " + family_label)
    lines.append("  Total unique genes:         " + str(total_genes))
    lines.append("  Genes in your results:      " + str(genes_in_yours))
    lines.append("  Genes in " + family_label + " list:         " + str(genes_in_list))
    lines.append("")

    # ---- Check 1 ----
    both_count = (in_yours_mask & in_list_mask).sum()
    yours_only_count = (in_yours_mask & ~in_list_mask).sum()
    list_only_count = (~in_yours_mask & in_list_mask).sum()

    lines.append("CHECK 1: GENE LIST MEMBERSHIP")
    lines.append("-" * 70)
    lines.append("  In both (verified):         " + str(both_count))

    unexpected_label = ""
    if yours_only_count > 0:
        unexpected_label = "  *** UNEXPECTED ***"
    lines.append("  In your results only:       " + str(yours_only_count) + unexpected_label)
    lines.append("  In gene list only:          " + str(list_only_count) + "  (filtered out by DE cutoffs)")

    if genes_in_list > 0:
        retention_pct = both_count / genes_in_list * 100
    else:
        retention_pct = 0
    lines.append("  Retention rate:             " + str(round(retention_pct, 1)) + "% of gene list passed your filters")

    if yours_only_count == 0:
        lines.append("  Result:                     PASS")
    else:
        lines.append("  Result:                     FAIL")
    lines.append("")

    # ---- Check 2 ----
    overlap_count = 0
    n_prev_de = 0
    if "in_prev_DE_list" in table.columns:
        prev_de_mask = table["in_prev_DE_list"] == "yes"
        n_prev_de = prev_de_mask.sum()
        overlap_count = (in_yours_mask & prev_de_mask).sum()
        prev_de_only = (prev_de_mask & ~in_yours_mask).sum()
        yours_new = (in_yours_mask & ~prev_de_mask).sum()

        lines.append("CHECK 2: PREVIOUS STUDENT DE OVERLAP")
        lines.append("-" * 70)
        lines.append("  Previous student DE genes:   " + str(n_prev_de))
        lines.append("  Overlap with your results:   " + str(overlap_count))
        lines.append("  Previous DE only:            " + str(prev_de_only) + "  (didn't pass your cutoffs)")
        lines.append("  Your new DE genes:           " + str(yours_new) + "  (not in previous student's DE)")
        if n_prev_de > 0:
            repro_pct = overlap_count / n_prev_de * 100
            lines.append("  Reproducibility:             " + str(round(repro_pct)) + "% of previous DE genes reproduced")
        lines.append("")

    # ---- Check 3 ----
    if "in_prev_expression" in table.columns:
        prev_expr_mask = table["in_prev_expression"] == "yes"
        n_prev_expr = prev_expr_mask.sum()
        expr_overlap = (in_yours_mask & prev_expr_mask).sum()

        lines.append("CHECK 3: EXPRESSION DATA COMPARISON")
        lines.append("-" * 70)
        lines.append("  Previous student total genes: " + str(n_prev_expr))
        lines.append("  Overlap with your DE genes:   " + str(expr_overlap))

        if "your_baseMean" in table.columns and "prev_expr_baseMean" in table.columns:
            shared_mask = in_yours_mask & prev_expr_mask
            yours_bm = pd.to_numeric(table.loc[shared_mask, "your_baseMean"], errors="coerce")
            prev_bm = pd.to_numeric(table.loc[shared_mask, "prev_expr_baseMean"], errors="coerce")
            valid_mask = yours_bm.notna() & prev_bm.notna()
            if valid_mask.sum() >= 3:
                corr = yours_bm[valid_mask].corr(prev_bm[valid_mask])
                if corr > 0.7:
                    strength = "strong"
                elif corr > 0.3:
                    strength = "moderate"
                else:
                    strength = "weak"
                lines.append("  baseMean correlation:         " + str(round(corr, 3)) + "  (" + strength + ")")
        lines.append("")

    # ---- Check 4 ----
    agree_count = 0
    disagree_count = 0
    total_compared = 0
    pct_agree = 0
    if "direction_agreement" in table.columns:
        agree_count = (table["direction_agreement"] == "agree").sum()
        disagree_count = (table["direction_agreement"] == "DISAGREE").sum()
        total_compared = agree_count + disagree_count
        if total_compared > 0:
            pct_agree = agree_count / total_compared * 100
        else:
            pct_agree = 0

        lines.append("CHECK 4: LOG2FOLDCHANGE DIRECTION")
        lines.append("-" * 70)
        lines.append("  Genes compared:              " + str(total_compared))
        lines.append("  Same direction (agree):      " + str(agree_count) + "  (" + str(round(pct_agree)) + "%)")
        lines.append("  Opposite direction:          " + str(disagree_count))
        if pct_agree >= 80:
            lines.append("  Assessment:                  STRONG agreement")
        elif pct_agree >= 50:
            lines.append("  Assessment:                  Moderate agreement")
        else:
            lines.append("  Assessment:                  Low agreement — check contrast direction")
        lines.append("")

        if "your_log2FoldChange" in table.columns and "prev_lfc_log2FoldChange" in table.columns:
            yours_lfc = pd.to_numeric(table["your_log2FoldChange"], errors="coerce")
            prev_lfc = pd.to_numeric(table["prev_lfc_log2FoldChange"], errors="coerce")
            valid_mask = yours_lfc.notna() & prev_lfc.notna()
            if valid_mask.sum() >= 3:
                corr = yours_lfc[valid_mask].corr(prev_lfc[valid_mask])
                if corr > 0.7:
                    strength = "strong"
                elif corr > 0.3:
                    strength = "moderate"
                else:
                    strength = "weak"
                lines.append("  log2FC magnitude correlation: " + str(round(corr, 3)) + "  (" + strength + ")")
                lines.append("")

        if disagree_count > 0:
            lines.append("  Genes with OPPOSITE direction:")
            disagree_mask = table["direction_agreement"] == "DISAGREE"
            disagree_gene_ids = list(table[disagree_mask].index)
            show_count = min(15, len(disagree_gene_ids))
            for i in range(show_count):
                gid = disagree_gene_ids[i]
                ylfc = ""
                plfc = ""
                if "your_log2FoldChange" in table.columns:
                    ylfc = str(table.at[gid, "your_log2FoldChange"])
                if "prev_lfc_log2FoldChange" in table.columns:
                    plfc = str(table.at[gid, "prev_lfc_log2FoldChange"])
                lines.append("    " + gid + "  yours=" + ylfc + "  prev=" + plfc)
            if len(disagree_gene_ids) > 15:
                lines.append("    ... and " + str(len(disagree_gene_ids) - 15) + " more (see comparison TSV)")
            lines.append("")

    # ---- Check 5 ----
    confirmed_up = 0
    up_in_yours = 0
    pct_confirmed_up = 0
    if "prev_upregulated" in table.columns:
        prev_up_mask = table["prev_upregulated"] == "yes"
        n_prev_up = prev_up_mask.sum()
        up_in_yours = (in_yours_mask & prev_up_mask).sum()

        lines.append("CHECK 5: UPREGULATED GENE AGREEMENT")
        lines.append("-" * 70)
        lines.append("  Previous upregulated genes:  " + str(n_prev_up))
        lines.append("  Found in your results:       " + str(up_in_yours))

        if up_in_yours > 0 and "your_log2FoldChange" in table.columns:
            shared_up_mask = in_yours_mask & prev_up_mask
            yours_lfc = pd.to_numeric(table.loc[shared_up_mask, "your_log2FoldChange"], errors="coerce")
            confirmed_up = (yours_lfc > 0).sum()
            flipped_down = (yours_lfc <= 0).sum()
            total_up = confirmed_up + flipped_down
            if total_up > 0:
                pct_confirmed_up = confirmed_up / total_up * 100
            else:
                pct_confirmed_up = 0
            lines.append("  Confirmed upregulated:       " + str(confirmed_up) + "  (" + str(round(pct_confirmed_up)) + "%)")
            lines.append("  Flipped to downregulated:    " + str(flipped_down))
            if pct_confirmed_up >= 80:
                lines.append("  Assessment:                  STRONG upregulation agreement")
            elif pct_confirmed_up >= 50:
                lines.append("  Assessment:                  Moderate — review flipped genes")
            else:
                lines.append("  Assessment:                  Low — check contrast direction (R vs L)")
        lines.append("")

    # ---- Check 6 ----
    down_in_yours = 0
    confirmed_down = 0
    pct_confirmed_dn = 0
    if "prev_downregulated" in table.columns:
        prev_down_mask = table["prev_downregulated"] == "yes"
        n_prev_down = prev_down_mask.sum()
        down_in_yours = (in_yours_mask & prev_down_mask).sum()

        lines.append("CHECK 6: DOWNREGULATED GENE AGREEMENT")
        lines.append("-" * 70)
        lines.append("  Previous downregulated genes: " + str(n_prev_down))
        lines.append("  Found in your results:       " + str(down_in_yours))

        if down_in_yours > 0 and "your_log2FoldChange" in table.columns:
            shared_down_mask = in_yours_mask & prev_down_mask
            yours_lfc_dn = pd.to_numeric(table.loc[shared_down_mask, "your_log2FoldChange"], errors="coerce")
            confirmed_down = (yours_lfc_dn < 0).sum()
            flipped_up = (yours_lfc_dn >= 0).sum()
            total_dn = confirmed_down + flipped_up
            if total_dn > 0:
                pct_confirmed_dn = confirmed_down / total_dn * 100
            else:
                pct_confirmed_dn = 0
            lines.append("  Confirmed downregulated:     " + str(confirmed_down) + "  (" + str(round(pct_confirmed_dn)) + "%)")
            lines.append("  Flipped to upregulated:      " + str(flipped_up))
            if pct_confirmed_dn >= 80:
                lines.append("  Assessment:                  STRONG downregulation agreement")
            elif pct_confirmed_dn >= 50:
                lines.append("  Assessment:                  Moderate — review flipped genes")
            else:
                lines.append("  Assessment:                  Low — check contrast direction (R vs L)")
        lines.append("")

    # ---- Contrast verification note ----
    lines.append("CONTRAST VERIFICATION")
    lines.append("-" * 70)
    lines.append("  Both analyses use: contrast = R vs L (Root vs Leaf)")
    lines.append("    Positive log2FC  = upregulated in ROOT")
    lines.append("    Negative log2FC  = upregulated in LEAF (= downregulated in ROOT)")
    lines.append("  Previous student (R DESeq2): results(dds, contrast=c('condition','R','L'))")
    lines.append("  Your pipeline (PyDESeq2):    contrast=['condition', 'R', 'L']")
    lines.append("  Directions should match across studies.")
    lines.append("")

    # ---- Top DE genes in your results ----
    if "your_padj" in table.columns:
        lines.append("TOP DE GENES IN YOUR RESULTS")
        lines.append("-" * 70)

        yours_padj = pd.to_numeric(table["your_padj"], errors="coerce")
        top_mask = in_yours_mask & yours_padj.notna()
        top_table = table[top_mask].copy()
        top_table["_padj"] = pd.to_numeric(top_table["your_padj"], errors="coerce")
        top_table = top_table.nsmallest(15, "_padj")

        # Build header line
        header = "  " + "gene_id".ljust(18) + "your_log2FC".rjust(12) + "your_padj".rjust(12) + "your_dir".rjust(10)
        if "prev_lfc_log2FoldChange" in table.columns:
            header = header + "prev_log2FC".rjust(12) + "agree?".rjust(10)
        lines.append(header)
        lines.append("  " + "-" * (len(header) - 2))

        for gene_id in top_table.index:
            # Get your values
            ylfc = ""
            ypadj = ""
            ydir = ""
            if "your_log2FoldChange" in top_table.columns:
                ylfc = top_table.at[gene_id, "your_log2FoldChange"]
            if "your_padj" in top_table.columns:
                ypadj = top_table.at[gene_id, "your_padj"]
            if "your_direction" in top_table.columns:
                ydir = top_table.at[gene_id, "your_direction"]

            # Format the numbers nicely
            try:
                if ylfc != "":
                    ylfc_str = format(float(ylfc), "+.2f")
                else:
                    ylfc_str = ""
            except (ValueError, TypeError):
                ylfc_str = str(ylfc)

            try:
                if ypadj != "":
                    ypadj_str = format(float(ypadj), ".2e")
                else:
                    ypadj_str = ""
            except (ValueError, TypeError):
                ypadj_str = str(ypadj)

            line = "  " + gene_id.ljust(18) + ylfc_str.rjust(12) + ypadj_str.rjust(12) + str(ydir).rjust(10)

            # Add previous student's values if available
            if "prev_lfc_log2FoldChange" in table.columns:
                plfc = ""
                agr = ""
                if "prev_lfc_log2FoldChange" in top_table.columns:
                    plfc = top_table.at[gene_id, "prev_lfc_log2FoldChange"]
                if "direction_agreement" in top_table.columns:
                    agr = top_table.at[gene_id, "direction_agreement"]

                try:
                    if plfc != "":
                        plfc_str = format(float(plfc), "+.2f")
                    else:
                        plfc_str = ""
                except (ValueError, TypeError):
                    plfc_str = str(plfc)

                line = line + plfc_str.rjust(12) + str(agr).rjust(10)

            lines.append(line)
        lines.append("")

    # ---- Final verdict ----
    lines.append("=" * 70)
    lines.append("FINAL VERDICT")
    lines.append("=" * 70)

    if yours_only_count == 0:
        lines.append("  Check 1 (gene list):       PASS")
    else:
        lines.append("  Check 1 (gene list):       FAIL")

    if "in_prev_DE_list" in table.columns:
        lines.append("  Check 2 (DE overlap):      " + str(overlap_count) + "/" + str(n_prev_de) + " reproduced")

    if "direction_agreement" in table.columns:
        lines.append("  Check 4 (direction):       "
                      + str(agree_count) + "/" + str(total_compared)
                      + " agree (" + str(round(pct_agree)) + "%)")

    if "prev_upregulated" in table.columns and up_in_yours > 0:
        lines.append("  Check 5 (upregulated):     "
                      + str(confirmed_up) + "/" + str(up_in_yours)
                      + " confirmed (" + str(round(pct_confirmed_up)) + "%)")

    if "prev_downregulated" in table.columns and down_in_yours > 0:
        lines.append("  Check 6 (downregulated):   "
                      + str(confirmed_down) + "/" + str(down_in_yours)
                      + " confirmed (" + str(round(pct_confirmed_dn)) + "%)")

    lines.append("  Contrast match:            CONFIRMED (both R vs L)")
    lines.append("")
    lines.append("FILES PRODUCED")
    lines.append("  Comparison table:  " + str(args.output))
    lines.append("  This summary:      " + str(summary_path))
    lines.append("=" * 70)

    # Write all lines to the file
    file_handle = open(summary_path, "w")
    file_handle.write("\n".join(lines) + "\n")
    file_handle.close()

    # Also print to the terminal
    for line in lines:
        print(line)


# =======================================================================
# MAIN: This is where the script actually starts running
# =======================================================================


def main():
    """
    The main function:
    1. Reads command-line arguments (file paths)
    2. Runs all 6 checks
    3. Builds the comparison table
    4. Writes the summary
    """

    # ---- Set up command-line argument parsing ----
    # This is what lets you run: python verify_genelist.py --results myfile.tsv --gene-list genes.txt
    parser = argparse.ArgumentParser(
        description="Verify pipeline results with 6 checks + output comparison table"
    )

    # Required arguments (you MUST provide these)
    parser.add_argument("--results", required=True,
                        help="Your results file (TSV)")
    parser.add_argument("--gene-list", required=True,
                        help="Original gene list (P450_list_RefSeq.txt)")

    # Optional arguments (you CAN provide these, but don't have to)
    parser.add_argument("--prev-de", default=None,
                        help="Check 2: Previous student's DE gene list")
    parser.add_argument("--prev-expression", default=None,
                        help="Check 3: Previous student's expression data")
    parser.add_argument("--prev-log2fold", default=None,
                        help="Check 4: Previous student's log2fold data")
    parser.add_argument("--prev-upregulated", default=None,
                        help="Check 5: Previous student's upregulated gene list")
    parser.add_argument("--prev-downregulated", default=None,
                        help="Check 6: Previous student's downregulated gene list")
    parser.add_argument("-o", "--output", default=None,
                        help="Where to save the comparison table (TSV file)")
    # CHANGED: 2026-04-01 — Added --gene-family so the summary header says
    # "CYP VERIFICATION SUMMARY" or "OMT VERIFICATION SUMMARY" instead of
    # the old hardcoded "P450 VERIFICATION SUMMARY".
    parser.add_argument("--gene-family", default="CYP",
                        help="Gene family being verified (CYP or OMT). "
                             "Used in the summary header. Default: CYP")

    # Parse (read) the arguments
    args = parser.parse_args()

    # Convert file paths to Path objects
    results_path = Path(args.results)
    list_path = Path(args.gene_list)

    # ---- Check that required files exist ----
    if not results_path.exists():
        print("ERROR: results file not found: " + str(results_path))
        sys.exit(1)

    if not list_path.exists():
        print("ERROR: gene list not found: " + str(list_path))
        sys.exit(1)

    # ---- Run Check 1 (always runs) ----
    passed = check1_gene_list(results_path, list_path)

    # ---- Run Check 2 (only if the user provided --prev-de) ----
    if args.prev_de is not None:
        prev_de_path = Path(args.prev_de)
        if prev_de_path.exists():
            check2_prev_de_overlap(results_path, prev_de_path)
        else:
            print("  WARNING: " + str(prev_de_path) + " not found, skipping Check 2")
            print()

    # ---- Run Check 3 (only if the user provided --prev-expression) ----
    if args.prev_expression is not None:
        prev_expr_path = Path(args.prev_expression)
        if prev_expr_path.exists():
            check3_expression(results_path, prev_expr_path)
        else:
            print("  WARNING: " + str(prev_expr_path) + " not found, skipping Check 3")
            print()

    # ---- Run Check 4 (only if the user provided --prev-log2fold) ----
    if args.prev_log2fold is not None:
        prev_lfc_path = Path(args.prev_log2fold)
        if prev_lfc_path.exists():
            check4_log2fold(results_path, prev_lfc_path)
        else:
            print("  WARNING: " + str(prev_lfc_path) + " not found, skipping Check 4")
            print()

    # ---- Run Check 5 (only if the user provided --prev-upregulated) ----
    if args.prev_upregulated is not None:
        prev_upreg_path = Path(args.prev_upregulated)
        if prev_upreg_path.exists():
            check5_upregulated(results_path, prev_upreg_path)
        else:
            print("  WARNING: " + str(prev_upreg_path) + " not found, skipping Check 5")
            print()

    # ---- Run Check 6 (only if the user provided --prev-downregulated) ----
    if args.prev_downregulated is not None:
        prev_downreg_path = Path(args.prev_downregulated)
        if prev_downreg_path.exists():
            check6_downregulated(results_path, prev_downreg_path)
        else:
            print("  WARNING: " + str(prev_downreg_path) + " not found, skipping Check 6")
            print()

    # ---- Build and save the comparison table ----
    if args.output is not None:
        print("=" * 60)
        print("  Building comparison table...")
        print("=" * 60)

        table = build_comparison_table(
            results_path, list_path,
            prev_de_path=args.prev_de,
            prev_expr_path=args.prev_expression,
            prev_lfc_path=args.prev_log2fold,
            prev_upreg_path=args.prev_upregulated,
            prev_downreg_path=args.prev_downregulated,
        )

        # Create the output directory if it doesn't exist
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Save the table as a TSV file
        table.to_csv(output_path, sep="\t")

        # Count how many genes are in your results
        genes_in_yours = (table["in_your_results"] == "yes").sum()
        total_genes = len(table)
        print("  Rows: " + str(total_genes) + " genes (" + str(genes_in_yours) + " in your results)")
        print()
        print("  Saved: " + str(output_path))
        print()

        # Write the summary text file
        summary_filename = output_path.stem + "_SUMMARY.txt"
        summary_path = output_path.with_name(summary_filename)
        write_summary(table, summary_path, args)
        print("  Saved: " + str(summary_path))
        print()

    # ---- Done! ----
    print("=" * 60)
    print("  VERIFICATION COMPLETE")
    print("=" * 60)

    # Exit with code 0 (success) if Check 1 passed, or 1 (failure) if it didn't
    if passed:
        sys.exit(0)
    else:
        sys.exit(1)


# This line means: "only run main() if this script is being run directly"
# (not if it's being imported by another script)
if __name__ == "__main__":
    main()
