#!/usr/bin/env python3
"""
Filter a count matrix to genes in a gene list, run PyDESeq2 to compute
expression stats (baseMean, log2FoldChange, padj), and optionally filter
with standard DE cutoffs. Output is ready for step 3 plots.

How PyDESeq2 computes the stats (replicated here):
  1. Load raw counts (genes x samples) + sample metadata (condition column)
  2. Estimate size factors (normalize for sequencing depth differences)
  3. Estimate gene-wise dispersions (variance model per gene)
  4. Fit a negative binomial GLM per gene: counts ~ condition
  5. Wald test → log2FoldChange, pvalue per gene
  6. Benjamini-Hochberg correction → padj (adjusted for multiple testing)

Usage — run PyDESeq2 fresh on candidate genes (recommended):
  python scripts/filter_count_by_genelist.py \
      --counts 03_count_tables/00_1_DC/gene_count_matrix.tsv \
      --gene-list 07_NRdatabase/sukman_database/P450_list_RefSeq.txt \
      --metadata 03_count_tables/00_1_DC/sample_metadata.tsv \
      --padj-cutoff 0.05 --lfc-cutoff 2.0 \
      --output 07_NRdatabase/sukman_database/geneious_candidates.tsv

Usage — use pre-computed step 1 results instead:
  python scripts/filter_count_by_genelist.py \
      --counts 03_count_tables/00_1_DC/gene_count_matrix.tsv \
      --gene-list 07_NRdatabase/sukman_database/P450_list_RefSeq.txt \
      --deseq 06_analysis/pydeseq2_DC_step1_unfiltered/pydeseq2_results_UNFILTERED.tsv \
      --padj-cutoff 0.05 --lfc-cutoff 2.0 \
      --output 07_NRdatabase/sukman_database/geneious_candidates.tsv

  Then feed into step 3:
    sbatch scripts/run_pydeseq2_step3_plots.sbatch DC \
        /path/to/geneious_candidates.tsv
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
import numpy as np


def load_gene_ids(filepath):
    """Load gene IDs from .txt (one per line), .csv, or .tsv."""
    path_str = str(filepath)

    if path_str.endswith(".txt"):
        with open(filepath) as fh:
            ids = [line.strip() for line in fh
                   if line.strip() and not line.startswith("#")]
        return list(dict.fromkeys(ids))

    sep = "\t" if path_str.endswith((".tsv", ".tab")) else ","
    df = pd.read_csv(filepath, sep=sep)

    if "gene_id" not in df.columns:
        if df.columns[0].startswith("LOC") or df.shape[1] == 1:
            df = df.rename(columns={df.columns[0]: "gene_id"})
        else:
            raise ValueError(
                f"{filepath} must have a 'gene_id' column "
                f"(found: {list(df.columns)})"
            )

    return df["gene_id"].astype(str).drop_duplicates().tolist()


def run_pydeseq2(counts_df, metadata_path, contrast_a="R", contrast_b="L"):
    """Run PyDESeq2 on a count matrix and return the results DataFrame.

    Uses the FULL count matrix for size factor and dispersion estimation
    (statistically correct), then returns stats for all genes.
    """
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    meta = pd.read_csv(metadata_path, sep="\t", index_col=0)

    shared = sorted(set(counts_df.columns) & set(meta.index))
    if not shared:
        print("  ERROR: no samples in common between count matrix and metadata")
        sys.exit(1)

    counts_aligned = counts_df[shared].T
    meta_aligned = meta.loc[shared]

    print(f"  Running PyDESeq2: {len(counts_df)} genes x {len(shared)} samples")
    print(f"  Contrast: condition {contrast_a} vs {contrast_b}")
    print()

    dds = DeseqDataSet(
        counts=counts_aligned,
        metadata=meta_aligned,
        design_factors="condition",
    )
    dds.deseq2()

    stats = DeseqStats(dds, contrast=["condition", contrast_a, contrast_b])
    stats.summary()

    results = stats.results_df
    results.index.name = "gene_id"
    return results


def merge_precomputed_deseq(result, deseq_path):
    """Merge pre-computed PyDESeq2 stats from step 1 onto the result table."""
    deseq = pd.read_csv(deseq_path, sep="\t", comment="#")
    if "gene_id" not in deseq.columns:
        deseq = deseq.rename(columns={deseq.columns[0]: "gene_id"})
    deseq["gene_id"] = deseq["gene_id"].astype(str)
    deseq = deseq.set_index("gene_id")

    stat_cols = ["baseMean", "log2FoldChange", "lfcSE", "stat",
                 "pvalue", "padj"]
    for col in stat_cols:
        if col in deseq.columns:
            result[col] = result.index.map(
                lambda g, c=col: deseq.at[g, c]
                if g in deseq.index else np.nan
            )

    matched = result["baseMean"].notna().sum() if "baseMean" in result.columns else 0
    print(f"  Genes with DESeq2 stats: {matched} / {len(result)}")
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Filter count matrix by gene list + run/merge PyDESeq2",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--counts", required=True,
                        help="Count matrix TSV (genes x samples)")
    parser.add_argument("--gene-list", required=True,
                        help="Gene list: .txt (one LOC per line), .csv, or .tsv")
    parser.add_argument("--metadata", default=None,
                        help="Sample metadata TSV — triggers a fresh PyDESeq2 "
                             "run on the full count matrix")
    parser.add_argument("--deseq", default=None,
                        help="Pre-computed PyDESeq2 results TSV (alternative "
                             "to --metadata; merges existing stats)")
    parser.add_argument("--contrast-a", default="R",
                        help="Condition A for contrast (default: R = root)")
    parser.add_argument("--contrast-b", default="L",
                        help="Condition B for contrast (default: L = leaf)")
    parser.add_argument("--padj-cutoff", type=float, default=None,
                        help="Filter: keep genes with padj < this (e.g. 0.05)")
    parser.add_argument("--lfc-cutoff", type=float, default=None,
                        help="Filter: keep genes with |log2FC| > this (e.g. 2.0)")
    parser.add_argument("-o", "--output", required=True,
                        help="Output TSV — curated candidate list")
    args = parser.parse_args()

    counts_path = Path(args.counts)
    list_path = Path(args.gene_list)
    if not counts_path.exists():
        print(f"ERROR: count matrix not found: {counts_path}")
        sys.exit(1)
    if not list_path.exists():
        print(f"ERROR: gene list not found: {list_path}")
        sys.exit(1)
    if args.metadata and not Path(args.metadata).exists():
        print(f"ERROR: metadata not found: {args.metadata}")
        sys.exit(1)
    if args.deseq and not Path(args.deseq).exists():
        print(f"ERROR: DESeq2 file not found: {args.deseq}")
        sys.exit(1)

    print("=" * 60)
    print("  Filter count matrix by gene list")
    print("=" * 60)
    print(f"  Count matrix: {counts_path}")
    print(f"  Gene list:    {list_path}")
    if args.metadata:
        print(f"  Metadata:     {args.metadata}  (will run PyDESeq2 fresh)")
    elif args.deseq:
        print(f"  DESeq2:       {args.deseq}  (pre-computed)")
    if args.padj_cutoff or args.lfc_cutoff:
        padj_str = f"padj < {args.padj_cutoff}" if args.padj_cutoff else ""
        lfc_str = f"|log2FC| > {args.lfc_cutoff}" if args.lfc_cutoff else ""
        filt = " AND ".join(x for x in [padj_str, lfc_str] if x)
        print(f"  DE filter:    {filt}")
    print(f"  Output:       {args.output}")
    print()

    # ── Load gene list ──
    gene_ids = load_gene_ids(list_path)
    print(f"  Genes in list (unique): {len(gene_ids)}")

    # ── Load count matrix ──
    counts = pd.read_csv(counts_path, sep="\t", index_col=0)
    counts = counts[~counts.index.str.startswith("N_")]
    sample_cols = list(counts.columns)
    print(f"  Genes in count matrix:  {counts.shape[0]}")
    print(f"  Samples:                {counts.shape[1]}")
    print(f"    {sample_cols}")
    print()

    # ── Match gene list against count matrix ──
    found = [g for g in gene_ids if g in counts.index]
    missing = [g for g in gene_ids if g not in counts.index]

    print(f"  Matched in count matrix:  {len(found)}")
    if missing:
        print(f"  NOT in count matrix:      {len(missing)}")
        if len(missing) <= 20:
            for g in missing:
                print(f"    - {g}")
        else:
            for g in missing[:10]:
                print(f"    - {g}")
            print(f"    ... and {len(missing) - 10} more")
    print()

    if not found:
        print("  ERROR: No genes matched. Check gene_id format (LOC...).")
        sys.exit(1)

    # ── Build result: filtered counts ──
    result = counts.loc[found].copy()
    result.index.name = "gene_id"
    result["total_counts"] = result[sample_cols].sum(axis=1)
    result["mean_counts"] = result[sample_cols].mean(axis=1).round(1)

    # ── Get DESeq2 stats ──
    if args.metadata:
        print("  ── Running PyDESeq2 ──")
        deseq_results = run_pydeseq2(
            counts, args.metadata, args.contrast_a, args.contrast_b
        )
        stat_cols = ["baseMean", "log2FoldChange", "lfcSE", "stat",
                     "pvalue", "padj"]
        for col in stat_cols:
            if col in deseq_results.columns:
                result[col] = result.index.map(
                    lambda g, c=col: deseq_results.at[g, c]
                    if g in deseq_results.index else np.nan
                )
        matched = result["baseMean"].notna().sum()
        print(f"  Genes with DESeq2 stats: {matched} / {len(result)}")
        print()

    elif args.deseq:
        print("  ── Merging pre-computed DESeq2 stats ──")
        result = merge_precomputed_deseq(result, args.deseq)
        print()

    # ── Apply DE filter if requested ──
    has_stats = "padj" in result.columns and "log2FoldChange" in result.columns
    filtered_out = 0

    if has_stats and (args.padj_cutoff or args.lfc_cutoff):
        n_before = len(result)
        mask = pd.Series(True, index=result.index)
        if args.padj_cutoff:
            mask &= result["padj"].notna() & (result["padj"] < args.padj_cutoff)
        if args.lfc_cutoff:
            mask &= result["log2FoldChange"].notna() & (result["log2FoldChange"].abs() > args.lfc_cutoff)
        result = result[mask].copy()
        filtered_out = n_before - len(result)
        print(f"  DE filter applied:")
        print(f"    Before: {n_before} genes")
        print(f"    After:  {len(result)} genes  ({filtered_out} removed)")

        if has_stats and len(result) > 0:
            n_up = (result["log2FoldChange"] > 0).sum()
            n_down = (result["log2FoldChange"] <= 0).sum()
            print(f"    Upregulated (root):  {n_up}")
            print(f"    Downregulated (leaf): {n_down}")
        print()

    if len(result) == 0:
        print("  WARNING: No genes passed the DE filter.")
        print("  Try relaxing cutoffs: --padj-cutoff 0.1 --lfc-cutoff 1.0")
        print()

    result = result.sort_values("total_counts", ascending=False)

    nonzero = (result["total_counts"] > 0).sum()
    zero = (result["total_counts"] == 0).sum()
    print(f"  Candidates with expression: {nonzero}")
    if zero:
        print(f"  Candidates with zero counts: {zero}")
    print()

    # ── Save ──
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, sep="\t")

    print(f"  Top genes by total counts:")
    top = result.head(min(10, len(result)))
    for gid in top.index:
        total = int(top.at[gid, "total_counts"])
        mean = top.at[gid, "mean_counts"]
        stats = ""
        if "padj" in top.columns and pd.notna(top.at[gid, "padj"]):
            stats = (f"  padj={top.at[gid, 'padj']:.2e}"
                     f"  log2FC={top.at[gid, 'log2FoldChange']:+.2f}")
        print(f"    {gid}  total={total:>8,}  mean={mean:>8.1f}{stats}")
    print()
    print(f"  Saved: {args.output}")
    print(f"  Rows: {len(result)} genes x {len(sample_cols)} samples")
    if has_stats:
        print()
        print("  Ready for step 3 plots:")
        print(f"    sbatch scripts/run_pydeseq2_step3_plots.sbatch DC {args.output}")
    print("=" * 60)


if __name__ == "__main__":
    main()
