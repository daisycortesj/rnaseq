#!/usr/bin/env python3
"""
Intersect a gene-family master list with PyDESeq2 (count matrix) results.

Supports any gene-family database (CYP, UGT, GST, etc.) via --database.
Optionally compares your list against a previous student's list (--prev-list)
and tags each gene as 'both', 'hmmer_gtf_only', or 'geneious_only'
(label names are configurable with --current-name / --prev-name).

Usage — CYP (default, backward-compatible):
  python scripts/cyp_intersect_pydeseq2.py \
      --gene-list 07_NRdatabase/cyp450_database/cyp_master_list.csv \
      --deseq 06_analysis/pydeseq2_DC_step1_unfiltered/pydeseq2_results_UNFILTERED.tsv \
      --gff 04_reference/dc_genomic.gtf \
      --output 07_NRdatabase/cyp450_database/cyp_expressed_list.tsv

Usage — Geneious list separately (plain .txt, one LOC per line):
  python scripts/cyp_intersect_pydeseq2.py \
      --gene-list 07_NRdatabase/sukman_database/P450_list_RefSeq.txt \
      --deseq 06_analysis/pydeseq2_DC_step1_unfiltered/pydeseq2_results_UNFILTERED.tsv \
      --gff 04_reference/dc_genomic.gtf \
      --evidence all \
      --output 07_NRdatabase/cyp450_database/geneious_expressed_list.tsv

Usage — different database (e.g. UGT):
  python scripts/cyp_intersect_pydeseq2.py \
      --database UGT \
      --gene-list 07_NRdatabase/ugt_database/ugt_master_list.csv \
      --deseq 06_analysis/pydeseq2_DC_step1_unfiltered/pydeseq2_results_UNFILTERED.tsv \
      --evidence Both,HMMER_only,keyword_only \
      --output 07_NRdatabase/ugt_database/ugt_expressed_list.tsv

Next steps:
  1. cyp_extract_proteins.py           — extract protein FASTA for these genes
  2. blastp_discoveryfilter.sbatch     — BLAST proteins against swissprot
  3. run_combine_filter.sbatch ... cyp — combine BLAST + expression, filter DE
  4. run_pydeseq2_step3_plots.sbatch   — heatmap, volcano, MA plots
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
import numpy as np

# Default paths (override with CLI)
CYP_LIST = "07_NRdatabase/cyp450_database/cyp_master_list.csv"
DESEQ_FILE = "06_analysis/pydeseq2_DC_step1_unfiltered/pydeseq2_results_UNFILTERED.tsv"
GFF_FILE = "04_reference/dc_genomic.gtf"
OUTPUT_FILE = "07_NRdatabase/cyp450_database/cyp_expressed_list.tsv"

DEFAULT_EVIDENCE = "Both,HMMER_only"


def _parse_attrs_auto(attr_string):
    """Auto-detect GFF3 vs GTF attribute format and parse to dict."""
    attrs = {}
    if '=' in attr_string and '"' not in attr_string.split('=')[0]:
        for item in attr_string.split(';'):
            item = item.strip()
            if '=' in item:
                k, v = item.split('=', 1)
                attrs[k] = v
    else:
        for item in attr_string.split(';'):
            item = item.strip()
            if item:
                parts = item.split(' ', 1)
                if len(parts) == 2:
                    attrs[parts[0]] = parts[1].strip('"')
    return attrs


def build_gene_to_protein_map(gff_path, gene_ids):
    """Build gene_id → protein_id mapping from GFF/GTF CDS features.

    Only maps gene_ids in the provided set (for speed on large files).
    """
    gff_path = Path(gff_path)
    if not gff_path.exists():
        return {}

    print(f"  Building gene→protein mapping from: {gff_path.name}")
    gene_to_prot = {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'CDS':
                continue
            attrs = _parse_attrs_auto(fields[8])
            gene_id = attrs.get('gene_id', '')
            protein_id = attrs.get('protein_id', '')
            if gene_id in gene_ids and protein_id and gene_id not in gene_to_prot:
                gene_to_prot[gene_id] = protein_id

    print(f"  Mapped {len(gene_to_prot)} gene→protein IDs")
    return gene_to_prot


def load_gene_list(csv_path, evidence_keep=None):
    """Load a gene-family list from CSV, TSV, or plain text.

    Supported formats:
      - .csv / .tsv with a 'gene_id' column (+ optional 'evidence' column)
      - .txt with one gene_id per line (no header)

    If the file has an 'evidence' column and evidence_keep is provided,
    rows are filtered to those values.  If there is no 'evidence' column
    (e.g. a simple gene list from a previous student), all rows are kept.
    """
    path_str = str(csv_path)

    if path_str.endswith(".txt"):
        with open(csv_path) as fh:
            ids = [line.strip() for line in fh if line.strip()
                   and not line.startswith("#")]
        df = pd.DataFrame({"gene_id": ids})
        df["gene_id"] = df["gene_id"].astype(str)
        return df

    sep = "\t" if path_str.endswith((".tsv", ".tab")) else ","
    df = pd.read_csv(csv_path, sep=sep)

    if "gene_id" not in df.columns:
        if df.columns[0].startswith("LOC") or df.shape[1] == 1:
            df = df.rename(columns={df.columns[0]: "gene_id"})
        else:
            raise ValueError(
                f"{csv_path} must have a 'gene_id' column "
                f"(found: {list(df.columns)})"
            )

    if "evidence" in df.columns and evidence_keep:
        df = df[df["evidence"].isin(evidence_keep)].copy()

    df["gene_id"] = df["gene_id"].astype(str)
    return df


def load_deseq(deseq_path):
    """Load PyDESeq2 or combined annotated TSV."""
    df = pd.read_csv(deseq_path, sep="\t", comment="#")
    if "gene_id" not in df.columns:
        df = df.rename(columns={df.columns[0]: "gene_id"})
    df["gene_id"] = df["gene_id"].astype(str)
    return df


def compare_gene_lists(current_ids, previous_ids, db_name,
                       current_name="hmmer_gtf", prev_name="geneious"):
    """Print a comparison between two gene-ID sets and return tagged dict.

    Labels use the caller-supplied method names so output reads e.g.:
        both, hmmer_gtf_only, geneious_only
    Returns {gene_id: source_label} for the union of both sets.
    """
    shared = current_ids & previous_ids
    cur_only = current_ids - previous_ids
    prev_only = previous_ids - current_ids
    union = current_ids | previous_ids

    tag_both = "both"
    tag_cur = f"{current_name}_only"
    tag_prev = f"{prev_name}_only"

    print(f"  ── {db_name} list comparison ──")
    print(f"    {current_name}:    {len(current_ids):>5} genes")
    print(f"    {prev_name}:      {len(previous_ids):>5} genes")
    print(f"    {tag_both}:            {len(shared):>5} genes")
    print(f"    {tag_cur}:  {len(cur_only):>5} genes  (new finds)")
    print(f"    {tag_prev}:   {len(prev_only):>5} genes  (check these)")
    print(f"    union:           {len(union):>5} genes")
    print()

    source = {}
    for g in shared:
        source[g] = tag_both
    for g in cur_only:
        source[g] = tag_cur
    for g in prev_only:
        source[g] = tag_prev
    return source


def main():
    parser = argparse.ArgumentParser(
        description="Intersect a gene-family master list with PyDESeq2 results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--gene-list", "--cyp-list", dest="gene_list",
                        default=CYP_LIST,
                        help=f"Your gene-family master list CSV/TSV (default: {CYP_LIST})")
    parser.add_argument("--prev-list", default=None,
                        help="Previous student's gene list CSV/TSV for comparison")
    parser.add_argument("--current-name", default="hmmer_gtf",
                        help="Label for YOUR list method (default: hmmer_gtf)")
    parser.add_argument("--prev-name", default="geneious",
                        help="Label for previous student's method (default: geneious)")
    parser.add_argument("--database", "--db", default="CYP",
                        help="Gene-family label for output messages (default: CYP)")
    parser.add_argument("--evidence", default=DEFAULT_EVIDENCE,
                        help=f"Comma-separated evidence values to keep "
                             f"(default: {DEFAULT_EVIDENCE}). "
                             f"Use 'all' to skip evidence filtering.")
    parser.add_argument("--deseq", default=DESEQ_FILE,
                        help=f"PyDESeq2 results TSV (default: {DESEQ_FILE})")
    parser.add_argument("--gff", default=GFF_FILE,
                        help=f"GFF/GTF for gene→protein mapping (default: {GFF_FILE})")
    parser.add_argument("-o", "--output", default=None,
                        help=f"Output TSV (default: auto-named from --database)")
    args = parser.parse_args()

    db = args.database.upper()
    if args.output is None:
        args.output = OUTPUT_FILE if db == "CYP" else \
            f"07_NRdatabase/{db.lower()}_database/{db.lower()}_expressed_list.tsv"

    evidence_keep = None if args.evidence.lower() == "all" \
        else set(args.evidence.split(","))

    gene_path = Path(args.gene_list)
    deseq_path = Path(args.deseq)
    if not gene_path.exists():
        print(f"ERROR: gene list not found: {gene_path}")
        sys.exit(1)
    if not deseq_path.exists():
        print(f"ERROR: DESeq file not found: {deseq_path}")
        sys.exit(1)

    print("=" * 60)
    print(f"STEP 1: {db} master list ∩ PyDESeq2")
    print("=" * 60)
    print(f"  Database:    {db}")
    print(f"  Gene list:   {gene_path}")
    if args.prev_list:
        print(f"  Prev list:   {args.prev_list}")
    print(f"  Evidence:    {args.evidence}")
    print(f"  DESeq:       {deseq_path}")
    print(f"  Output:      {args.output}")
    print()

    # ── Load your gene list ──
    gene_df = load_gene_list(gene_path, evidence_keep)
    current_ids = set(gene_df["gene_id"])
    print(f"  {db} master list: {len(current_ids)} genes "
          f"(evidence filter: {args.evidence})")

    # ── Optionally load & compare previous student's list ──
    source_map = {g: args.current_name for g in current_ids}
    combined_ids = set(current_ids)

    if args.prev_list:
        prev_path = Path(args.prev_list)
        if not prev_path.exists():
            print(f"ERROR: previous list not found: {prev_path}")
            sys.exit(1)
        prev_df = load_gene_list(prev_path, evidence_keep=None)
        prev_ids = set(prev_df["gene_id"])
        print(f"  Previous student list: {len(prev_ids)} genes")
        print()

        source_map = compare_gene_lists(current_ids, prev_ids, db,
                                        args.current_name, args.prev_name)
        combined_ids = current_ids | prev_ids

        # Build a merged DataFrame: your rows + previous-only rows
        prev_only_ids = prev_ids - current_ids
        if prev_only_ids:
            prev_only_df = prev_df[prev_df["gene_id"].isin(prev_only_ids)].copy()
            gene_df = pd.concat([gene_df, prev_only_df], ignore_index=True)

    # ── Intersect with PyDESeq2 ──
    deseq_df = load_deseq(deseq_path)
    deseq_ids = set(deseq_df["gene_id"])
    print(f"  Genes in DESeq2/count matrix:  {len(deseq_ids)} genes")

    overlap = combined_ids & deseq_ids
    print(f"  Intersection (keep these):     {len(overlap)} genes")
    missing = combined_ids - deseq_ids
    if missing:
        print(f"  {db} genes NOT in count matrix: {len(missing)} (dropped)")
    print()

    if not overlap:
        print("  ERROR: No overlap. Check gene_id format (LOC...).")
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(columns=["gene_id"]).to_csv(args.output, sep="\t", index=False)
        sys.exit(1)

    out = gene_df[gene_df["gene_id"].isin(overlap)].drop_duplicates(
        subset="gene_id"
    ).copy()
    deseq_lookup = deseq_df.set_index("gene_id")

    # Tag source (current / previous / both)
    out["list_source"] = out["gene_id"].map(source_map).fillna(args.current_name)

    # Attach expression columns from DESeq2
    expr_cols = ["baseMean", "log2FoldChange", "pvalue", "padj"]
    blast_cols = ["blast_description", "blast_species", "sseqid",
                  "pident", "qcovhsp", "evalue", "bitscore"]
    for col in expr_cols + blast_cols:
        if col in deseq_lookup.columns:
            out[col] = out["gene_id"].map(
                lambda g, c=col: deseq_lookup.at[g, c] if g in deseq_lookup.index else np.nan
            )

    # Direction label
    lfc = out.get("log2FoldChange", pd.Series(dtype=float))
    padj = out.get("padj", pd.Series(dtype=float))
    direction = pd.Series("no_data", index=out.index)
    valid = pd.notna(padj) & pd.notna(lfc)
    direction.loc[valid & (padj <= 0.05) & (lfc.abs() >= 1.0) & (lfc > 0)] = "root_up"
    direction.loc[valid & (padj <= 0.05) & (lfc.abs() >= 1.0) & (lfc <= 0)] = "leaf_up"
    direction.loc[valid & ~((padj <= 0.05) & (lfc.abs() >= 1.0))] = "ns"
    out["direction"] = direction

    # Ensure protein_id column is present and filled
    if "protein_id" not in out.columns:
        out["protein_id"] = np.nan

    n_missing_prot = out["protein_id"].isna().sum()
    if n_missing_prot > 0:
        print(f"  {n_missing_prot} genes missing protein_id — filling from GFF...")
        genes_need = set(out.loc[out["protein_id"].isna(), "gene_id"])
        gff_map = build_gene_to_protein_map(args.gff, genes_need)
        for idx, row in out.iterrows():
            if pd.isna(row["protein_id"]) and row["gene_id"] in gff_map:
                out.at[idx, "protein_id"] = gff_map[row["gene_id"]]
        filled = n_missing_prot - out["protein_id"].isna().sum()
        print(f"  Filled {filled}, still missing: {out['protein_id'].isna().sum()}")

    out = out.sort_values("gene_id").reset_index(drop=True)

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False)

    # ── Summary ──
    print()
    print("  Summary:")
    print(f"    Total {db} genes with expression: {len(out)}")
    print(f"    With protein_id:                 {out['protein_id'].notna().sum()}")
    if args.prev_list:
        source_tags = ["both",
                       f"{args.current_name}_only",
                       f"{args.prev_name}_only"]
        for src in source_tags:
            n = (out["list_source"] == src).sum()
            if n:
                print(f"    source {src:20s}: {n}")
    for d in ["root_up", "leaf_up", "ns", "no_data"]:
        n = (out["direction"] == d).sum()
        if n:
            print(f"    direction {d:10s}: {n}")
    print()
    print(f"  Saved: {args.output}")
    print("=" * 60)


if __name__ == "__main__":
    main()
