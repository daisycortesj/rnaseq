#!/usr/bin/env python3
# CHANGED: 2026-03-31 — Renamed from cyp_extract_proteins.py → gene_extract_proteins.py.
#   Generalized docstring and banner to support any gene family (CYP, OMT, etc.).
#   Auto-detects family name from input path for the step banner.
"""
Extract protein sequences for a gene family from the full protein FASTA.

Works for any gene family (CYP, OMT, UGT, etc.).
Reads an expressed list TSV (output of gene_intersect_pydeseq2.py) to get
the protein_id (XP_*) and gene_id (LOC*) values, then pulls matching
sequences from the full protein FASTA.

Output FASTA headers use the gene_id (e.g. >LOC108194654) so that
BLAST results join directly to gene_ids in combine_blast_deseq.py.

Usage:
  # CYP proteins:
  python scripts/08_gene_families/gene_extract_proteins.py \
      --expressed 07_NRdatabase/cyp450_database/cyp_expressed_list.tsv \
      --fasta 04_reference/GCF_001625215.2_DH1_v3.0_protein.faa \
      --output 07_NRdatabase/cyp450_database/cyp_proteins.fasta

  # OMT proteins:
  python scripts/08_gene_families/gene_extract_proteins.py \
      --expressed 07_NRdatabase/omt_database/omt_expressed_list.tsv \
      --fasta 04_reference/GCF_001625215.2_DH1_v3.0_protein.faa \
      --output 07_NRdatabase/omt_database/omt_proteins.fasta
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

EXPRESSED_FILE = "07_NRdatabase/cyp450_database/cyp_expressed_list.tsv"
PROTEIN_FASTA = "04_reference/GCF_001625215.2_DH1_v3.0_protein.faa"
OUTPUT_FASTA = "07_NRdatabase/cyp450_database/cyp_proteins.fasta"


def extract_proteins(expressed_path, fasta_path, output_path):
    """Pull protein sequences whose accession is in the expressed list.

    Rewrites FASTA headers to use gene_id (LOC*) so that BLAST results
    join directly to gene_ids in combine_blast_deseq.py and step 2/3.
    """

    df = pd.read_csv(expressed_path, sep="\t")
    if "protein_id" not in df.columns or "gene_id" not in df.columns:
        print("ERROR: expressed list must have gene_id and protein_id columns",
              file=sys.stderr)
        sys.exit(1)

    sub = df[["gene_id", "protein_id"]].dropna(subset=["protein_id"]).copy()
    sub["gene_id"] = sub["gene_id"].astype(str)
    sub["protein_id"] = sub["protein_id"].astype(str)

    # Build protein_id → gene_id map (also index by base accession without version)
    prot_to_gene = {}
    for _, row in sub.iterrows():
        pid = row["protein_id"]
        gid = row["gene_id"]
        prot_to_gene[pid] = gid
        base = pid.rsplit(".", 1)[0] if "." in pid else pid
        prot_to_gene.setdefault(base, gid)

    wanted = set(sub["protein_id"])
    wanted_bases = {p.rsplit(".", 1)[0] for p in wanted if "." in p}
    print(f"  Protein IDs to extract: {len(wanted)}")

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    found = set()
    writing = False
    written = 0

    with open(fasta_path) as fin, open(output_path, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                acc = line[1:].split()[0]
                acc_base = acc.rsplit(".", 1)[0] if "." in acc else acc
                if acc in wanted or acc_base in wanted_bases:
                    gene_id = prot_to_gene.get(acc) or prot_to_gene.get(acc_base, acc)
                    rest = line[1:].split(None, 1)
                    desc = rest[1] if len(rest) > 1 else "\n"
                    fout.write(f">{gene_id} {desc}")
                    writing = True
                    found.add(acc)
                    written += 1
                else:
                    writing = False
            elif writing:
                fout.write(line)

    still_missing = wanted - found - {a.rsplit(".", 1)[0] for a in found}
    print(f"  Sequences written: {written}")
    if still_missing:
        print(f"  WARNING: {len(still_missing)} protein_ids not found in FASTA:")
        for pid in sorted(still_missing)[:10]:
            print(f"    {pid}")

    return written


def main():
    parser = argparse.ArgumentParser(
        description="Extract gene-family protein sequences from full FASTA",
    )
    parser.add_argument("--expressed", default=EXPRESSED_FILE,
                        help=f"Expressed list TSV (default: {EXPRESSED_FILE})")
    parser.add_argument("--fasta", default=PROTEIN_FASTA,
                        help=f"Full protein FASTA (default: {PROTEIN_FASTA})")
    parser.add_argument("-o", "--output", default=OUTPUT_FASTA,
                        help=f"Output subset FASTA (default: {OUTPUT_FASTA})")
    args = parser.parse_args()

    if not Path(args.expressed).exists():
        print(f"ERROR: Expressed list not found: {args.expressed}")
        sys.exit(1)
    if not Path(args.fasta).exists():
        print(f"ERROR: Protein FASTA not found: {args.fasta}")
        sys.exit(1)

    # Auto-detect family name from the expressed list path for display
    family_label = Path(args.expressed).stem.split("_")[0].upper()

    print("=" * 60)
    print(f"STEP 2: Extract {family_label} protein sequences")
    print("=" * 60)
    print(f"  Expressed list: {args.expressed}")
    print(f"  Full FASTA:     {args.fasta}")
    print(f"  Output FASTA:   {args.output}")
    print()

    n = extract_proteins(args.expressed, args.fasta, args.output)
    print()
    print(f"  Saved: {args.output} ({n} sequences)")
    print("  FASTA headers use gene_id (LOC*) for direct matching.")
    print("=" * 60)


if __name__ == "__main__":
    main()
