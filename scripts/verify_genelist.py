#!/usr/bin/env python3
"""
Double-check that every gene in a results file appears in the original
gene list. Run this AFTER step 3 (or any path) to verify nothing snuck
in or went missing.

Usage:
  python scripts/verify_genelist.py \
      --results 07_NRdatabase/sukman_database/geneious_candidates.tsv \
      --gene-list 07_NRdatabase/sukman_database/P450_list_RefSeq.txt

  python scripts/verify_genelist.py \
      --results 07_NRdatabase/cyp450_database/cyp_expressed_list.tsv \
      --gene-list 07_NRdatabase/sukman_database/P450_list_RefSeq.txt

  # Also works on step 3 plot output gene lists:
  python scripts/verify_genelist.py \
      --results 06_analysis/pydeseq2_DC_step3_plots_*/cyp_gene_list.tsv \
      --gene-list 07_NRdatabase/sukman_database/P450_list_RefSeq.txt
"""

import argparse
import sys
from pathlib import Path


def load_ids_from_list(filepath):
    """Load gene IDs from .txt (one per line), .csv, or .tsv."""
    path_str = str(filepath)

    if path_str.endswith(".txt"):
        with open(filepath) as fh:
            ids = [line.strip() for line in fh
                   if line.strip() and not line.startswith("#")]
        return set(ids)

    sep = "\t" if path_str.endswith((".tsv", ".tab")) else ","
    import pandas as pd
    df = pd.read_csv(filepath, sep=sep)

    if "gene_id" in df.columns:
        return set(df["gene_id"].astype(str))
    if df.columns[0].startswith("LOC") or df.shape[1] == 1:
        return set(df.iloc[:, 0].astype(str))

    raise ValueError(f"Cannot find gene_id column in {filepath}: {list(df.columns)}")


def load_ids_from_results(filepath):
    """Load gene IDs from a results TSV (gene_id as index or first column)."""
    import pandas as pd
    df = pd.read_csv(filepath, sep="\t", comment="#")

    if "gene_id" in df.columns:
        return set(df["gene_id"].astype(str))
    if df.columns[0] == "Unnamed: 0" or df.iloc[:, 0].astype(str).str.startswith("LOC").any():
        return set(df.iloc[:, 0].astype(str))

    df2 = pd.read_csv(filepath, sep="\t", index_col=0, comment="#")
    if df2.index.astype(str).str.startswith("LOC").any():
        return set(df2.index.astype(str))

    raise ValueError(f"Cannot find gene IDs in {filepath}: {list(df.columns)}")


def main():
    parser = argparse.ArgumentParser(
        description="Verify that results genes match the original gene list",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--results", required=True,
                        help="Results file to verify (TSV from any path or step 3)")
    parser.add_argument("--gene-list", required=True,
                        help="Original gene list (.txt, .csv, .tsv)")
    args = parser.parse_args()

    results_path = Path(args.results)
    list_path = Path(args.gene_list)

    if not results_path.exists():
        print(f"ERROR: results file not found: {results_path}")
        sys.exit(1)
    if not list_path.exists():
        print(f"ERROR: gene list not found: {list_path}")
        sys.exit(1)

    results_ids = load_ids_from_results(results_path)
    list_ids = load_ids_from_list(list_path)

    in_both = results_ids & list_ids
    results_only = results_ids - list_ids
    list_only = list_ids - results_ids

    print("=" * 60)
    print("  Gene List Verification")
    print("=" * 60)
    print(f"  Results file: {results_path}")
    print(f"    Genes: {len(results_ids)}")
    print(f"  Gene list:    {list_path}")
    print(f"    Genes: {len(list_ids)}")
    print()
    print(f"  In BOTH (verified):       {len(in_both)}")
    print(f"  In results ONLY:          {len(results_only)}")
    print(f"  In gene list ONLY:        {len(list_only)}")
    print()

    if results_only:
        print("  *** UNEXPECTED: genes in results but NOT in gene list ***")
        for g in sorted(results_only):
            print(f"    - {g}")
        print()

    if list_only:
        print(f"  Genes in gene list but NOT in results ({len(list_only)}):")
        print("  (filtered out by DE cutoffs, or not in count matrix)")
        sorted_missing = sorted(list_only)
        if len(sorted_missing) <= 30:
            for g in sorted_missing:
                print(f"    - {g}")
        else:
            for g in sorted_missing[:15]:
                print(f"    - {g}")
            print(f"    ... and {len(sorted_missing) - 15} more")
        print()

    # Verdict
    if not results_only:
        pct = len(in_both) / len(list_ids) * 100 if list_ids else 0
        print(f"  PASS: All {len(results_ids)} genes in results are in the gene list")
        print(f"        ({pct:.1f}% of the gene list made it through filtering)")
    else:
        print(f"  FAIL: {len(results_only)} gene(s) in results are NOT in the gene list")
        print("        Something unexpected happened — investigate above genes")

    print("=" * 60)
    sys.exit(0 if not results_only else 1)


if __name__ == "__main__":
    main()
