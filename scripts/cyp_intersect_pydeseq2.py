#!/usr/bin/env python3
"""
Intersect CYP master list with PyDESeq2 (count matrix) results.

Keeps only genes that are in BOTH:
  - cyp_master_list.csv (Both + HMMER_only)
  - PyDESeq2 results (genes that have expression stats from the count matrix)

Output is your "CYP expressed list": CYPs in the count matrix with expression
stats and protein_id (needed for the BLAST step).

Usage:
  python scripts/cyp_intersect_pydeseq2.py \
      --cyp-list 07_NRdatabase/cyp450_database/cyp_master_list.csv \
      --deseq 06_analysis/pydeseq2_DC_step1_unfiltered/pydeseq2_results_UNFILTERED.tsv \
      --gff 04_reference/dc_genomic.gtf \
      --output 07_NRdatabase/cyp450_database/cyp_expressed_list.tsv

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

EVIDENCE_KEEP = {"Both", "HMMER_only"}


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


def load_cyp_master_list(csv_path, evidence_keep=None):
    """Load CYP master list, return DataFrame restricted to evidence_keep."""
    if evidence_keep is None:
        evidence_keep = EVIDENCE_KEEP
    df = pd.read_csv(csv_path)
    if "gene_id" not in df.columns or "evidence" not in df.columns:
        raise ValueError("cyp_master_list must have columns gene_id, evidence")
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


def main():
    parser = argparse.ArgumentParser(
        description="Intersect CYP master list with PyDESeq2 results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--cyp-list", default=CYP_LIST,
                        help=f"CYP master list CSV (default: {CYP_LIST})")
    parser.add_argument("--deseq", default=DESEQ_FILE,
                        help=f"PyDESeq2 results TSV (default: {DESEQ_FILE})")
    parser.add_argument("--gff", default=GFF_FILE,
                        help=f"GFF/GTF for gene→protein mapping fallback (default: {GFF_FILE})")
    parser.add_argument("-o", "--output", default=OUTPUT_FILE,
                        help=f"Output TSV (default: {OUTPUT_FILE})")
    args = parser.parse_args()

    cyp_path = Path(args.cyp_list)
    deseq_path = Path(args.deseq)
    if not cyp_path.exists():
        print(f"ERROR: CYP list not found: {cyp_path}")
        sys.exit(1)
    if not deseq_path.exists():
        print(f"ERROR: DESeq file not found: {deseq_path}")
        sys.exit(1)

    print("=" * 60)
    print("STEP 1: CYP master list ∩ PyDESeq2")
    print("=" * 60)
    print(f"  CYP list:  {cyp_path}")
    print(f"  DESeq:     {deseq_path}")
    print(f"  Output:    {args.output}")
    print()

    cyp_df = load_cyp_master_list(cyp_path)
    cyp_ids = set(cyp_df["gene_id"])
    print(f"  CYP master list (Both + HMMER_only): {len(cyp_ids)} genes")

    deseq_df = load_deseq(deseq_path)
    deseq_ids = set(deseq_df["gene_id"])
    print(f"  Genes in DESeq2/count matrix:         {len(deseq_ids)} genes")

    overlap = cyp_ids & deseq_ids
    print(f"  Intersection (keep these):            {len(overlap)} genes")
    missing = cyp_ids - deseq_ids
    if missing:
        print(f"  CYP genes NOT in count matrix:        {len(missing)} (dropped)")
    print()

    if not overlap:
        print("  ERROR: No overlap. Check gene_id format (LOC...).")
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(columns=["gene_id"]).to_csv(args.output, sep="\t", index=False)
        sys.exit(1)

    out = cyp_df[cyp_df["gene_id"].isin(overlap)].copy()
    deseq_lookup = deseq_df.set_index("gene_id")

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

    print()
    print("  Summary:")
    print(f"    Total CYP genes with expression: {len(out)}")
    print(f"    With protein_id:                 {out['protein_id'].notna().sum()}")
    for d in ["root_up", "leaf_up", "ns", "no_data"]:
        n = (out["direction"] == d).sum()
        if n:
            print(f"    direction {d:10s}: {n}")
    print()
    print(f"  Saved: {args.output}")
    print("=" * 60)


if __name__ == "__main__":
    main()
