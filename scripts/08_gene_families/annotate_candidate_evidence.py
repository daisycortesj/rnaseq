#!/usr/bin/env python3
"""
Add evidence columns to an existing CYP/OMT candidate TSV.

WHAT THIS DOES:
  Reads cyp_candidates_MF.tsv (or OMT) WITHOUT changing it.
  Looks up each gene in HMMER, BLAST genefinder, InterProScan, and BLASTX-nr.
  Writes a NEW file: cyp_candidates_MF_annotated.tsv  (original is left alone)

WHY:
  run_filter_genelist.sbatch only tags gene_family=CYP/OMT and DE stats.
  This script answers: which tool put this gene on the family list, and
  what name did BLASTX vs nr give it (if any)?

HOW TO RUN (on ARC, from the repo):
  conda activate rnaseq
  python scripts/08_gene_families/annotate_candidate_evidence.py --family both

OUTPUT (same folder as the originals; originals are NOT overwritten):
  /projects/tholl_lab_1/daisy_analysis/07_NRdatabase/00_2_MF/cyp_candidates_MF_annotated.tsv
  /projects/tholl_lab_1/daisy_analysis/07_NRdatabase/00_2_MF/omt_candidates_MF_annotated.tsv
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd


# ---------------------------------------------------------------------------
# Default HPC folders (override with flags if your paths differ)
# ---------------------------------------------------------------------------
DEFAULT_BASE = Path("/projects/tholl_lab_1/daisy_analysis")
# Annotated tables go here, next to cyp_candidates_MF.tsv / omt_candidates_MF.tsv
DEFAULT_OUT_DIR = Path(
    "/projects/tholl_lab_1/daisy_analysis/07_NRdatabase/00_2_MF"
)


def to_gene_id(raw_id):
    """Turn a protein/transcript ID into a Trinity GENE id.

    Examples:
      TRINITY_DN123_c0_g1_i2.p1  ->  TRINITY_DN123_c0_g1
      TRINITY_DN123_c0_g1_i2     ->  TRINITY_DN123_c0_g1
      TRINITY_DN123_c0_g1        ->  TRINITY_DN123_c0_g1
    """
    text = str(raw_id).strip().split()[0]
    # Strip TransDecoder ORF suffix (.p1, .p2, ...)
    text = re.sub(r"\.p[0-9]+$", "", text)
    # Strip Trinity isoform suffix (_i1, _i2, ...)
    text = re.sub(r"_i[0-9]+$", "", text)
    return text


def load_id_set(path):
    """Load a one-ID-per-line file and collapse every ID to gene level."""
    path = Path(path)
    if not path.is_file():
        print(f"  WARNING: missing {path}  (that evidence column will be no)")
        return set()
    gene_ids = set()
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line == "" or line.startswith("#"):
                continue
            gene_ids.add(to_gene_id(line))
    return gene_ids


def load_blastx_best_per_gene(path):
    """Keep one BLASTX hit per GENE: the isoform with the highest bitscore.

    BLASTX outfmt 6 columns (no header):
      qseqid sacc pident length evalue bitscore qcovs staxids stitle
    """
    path = Path(path)
    empty = {}
    if not path.is_file():
        print(f"  WARNING: missing {path}  (BLASTX name columns will be empty)")
        print("           BLASTX may still be running, or only the old 79-seq file exists.")
        return empty

    table = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=[
            "qseqid", "sacc", "pident", "length", "evalue",
            "bitscore", "qcovs", "staxids", "stitle",
        ],
    )
    table["gene_id"] = table["qseqid"].map(to_gene_id)
    table["bitscore"] = pd.to_numeric(table["bitscore"], errors="coerce")
    table = table.sort_values("bitscore", ascending=False)
    # One row per gene: first row after sorting by bitscore is the best
    table = table.drop_duplicates(subset=["gene_id"], keep="first")
    best = {}
    for _, row in table.iterrows():
        best[row["gene_id"]] = row
    print(f"  BLASTX names for {len(best)} genes from {path.name}")
    return best


def load_interpro_for_genes(summary_path, wanted_genes):
    """Collect InterProScan rows whose protein maps to a wanted gene.

    Summary columns:
      Protein, Database, Signature, Description,
      InterPro_ID, InterPro_description, GO_terms
    """
    summary_path = Path(summary_path)
    by_gene = {}
    if not summary_path.is_file():
        print(f"  WARNING: missing {summary_path}")
        return by_gene

    print(f"  Reading InterProScan summary (only genes in your candidate list)...")
    table = pd.read_csv(summary_path, sep="\t")
    if "Protein" not in table.columns:
        print(f"  WARNING: unexpected InterProScan header: {list(table.columns)}")
        return by_gene

    table["gene_id"] = table["Protein"].map(to_gene_id)
    table = table[table["gene_id"].isin(wanted_genes)]

    for gene_id, group in table.groupby("gene_id"):
        ipr_ids = []
        ipr_desc = []
        signatures = []
        for _, row in group.iterrows():
            ipr = row.get("InterPro_ID", "")
            desc = row.get("InterPro_description", "")
            db = row.get("Database", "")
            sig = row.get("Signature", "")
            dshort = row.get("Description", "")
            if pd.notna(ipr) and str(ipr) not in ("", "-", "NaN"):
                if str(ipr) not in ipr_ids:
                    ipr_ids.append(str(ipr))
                    if pd.notna(desc) and str(desc) not in ("", "-"):
                        ipr_desc.append(str(desc))
            if pd.notna(sig) and str(sig) not in ("", "-"):
                piece = f"{db}:{sig}"
                if pd.notna(dshort) and str(dshort) not in ("", "-"):
                    piece = piece + f" ({dshort})"
                if piece not in signatures:
                    signatures.append(piece)
        by_gene[gene_id] = {
            "interpro_ids": "|".join(ipr_ids),
            "interpro_descriptions": " | ".join(ipr_desc),
            "interpro_signatures": " | ".join(signatures[:12]),
        }
    print(f"  InterProScan hits for {len(by_gene)} candidate genes")
    return by_gene


def annotate_one(candidates_path, out_path, hmmer_ids, blast_ids, ips_by_gene, blastx_by_gene):
    """Copy every original column, then append evidence columns."""
    candidates_path = Path(candidates_path)
    if not candidates_path.is_file():
        print(f"ERROR: candidate table not found: {candidates_path}")
        sys.exit(1)

    result = pd.read_csv(candidates_path, sep="\t")
    if "gene_id" not in result.columns:
        result = result.rename(columns={result.columns[0]: "gene_id"})
    result["gene_id"] = result["gene_id"].astype(str)

    print(f"  Candidates: {len(result)} genes in {candidates_path.name}")

    in_hmmer = []
    in_blast = []
    in_ips = []
    sources = []
    ipr_ids_col = []
    ipr_desc_col = []
    ipr_sig_col = []
    bx_acc = []
    bx_pident = []
    bx_title = []

    for gene_id in result["gene_id"]:
        hmmer_hit = gene_id in hmmer_ids
        blast_hit = gene_id in blast_ids
        ips_hit = gene_id in ips_by_gene
        in_hmmer.append("yes" if hmmer_hit else "no")
        in_blast.append("yes" if blast_hit else "no")
        in_ips.append("yes" if ips_hit else "no")

        tags = []
        if hmmer_hit:
            tags.append("HMMER")
        if blast_hit:
            tags.append("BLAST_genefinder")
        if ips_hit:
            tags.append("InterProScan")
        if gene_id in blastx_by_gene:
            tags.append("BLASTX_nr")
        sources.append(";".join(tags) if tags else "none")

        if ips_hit:
            info = ips_by_gene[gene_id]
            ipr_ids_col.append(info["interpro_ids"])
            ipr_desc_col.append(info["interpro_descriptions"])
            ipr_sig_col.append(info["interpro_signatures"])
        else:
            ipr_ids_col.append("")
            ipr_desc_col.append("")
            ipr_sig_col.append("")

        if gene_id in blastx_by_gene:
            hit = blastx_by_gene[gene_id]
            bx_acc.append(hit["sacc"])
            bx_pident.append(hit["pident"])
            bx_title.append(hit["stitle"])
        else:
            bx_acc.append("")
            bx_pident.append("")
            bx_title.append("")

    result["evidence_hmmer"] = in_hmmer
    result["evidence_blast_genefinder"] = in_blast
    result["evidence_interproscan"] = in_ips
    result["evidence_sources"] = sources
    result["interpro_ids"] = ipr_ids_col
    result["interpro_descriptions"] = ipr_desc_col
    result["interpro_signatures"] = ipr_sig_col
    result["blastx_accession"] = bx_acc
    result["blastx_pident"] = bx_pident
    result["blastx_title"] = bx_title

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(out_path, sep="\t", index=False)
    print(f"  Saved: {out_path}")
    print(f"  Original NOT changed: {candidates_path}")
    n_named = sum(1 for t in bx_title if str(t).strip() != "")
    print(f"  Rows with a BLASTX name: {n_named} / {len(result)}")
    print()


def main():
    parser = argparse.ArgumentParser(
        description="Annotate CYP/OMT candidate TSVs with HMMER/BLAST/InterProScan/BLASTX evidence",
    )
    parser.add_argument(
        "--family",
        choices=["CYP", "OMT", "both"],
        default="both",
        help="Which candidate table to annotate (default: both)",
    )
    parser.add_argument(
        "--base-dir",
        default=str(DEFAULT_BASE),
        help="HPC project root (HMMER / BLAST / InterProScan / BLASTX inputs)",
    )
    parser.add_argument(
        "--out-dir",
        default=str(DEFAULT_OUT_DIR),
        help="Where to write *_annotated.tsv (default: 07_NRdatabase/00_2_MF)",
    )
    args = parser.parse_args()

    base = Path(args.base_dir)
    cand_dir = Path(args.out_dir)
    cand_dir.mkdir(parents=True, exist_ok=True)
    print(f"  Output folder: {cand_dir}")
    hmmer_dir = base / "06_analysis" / "hmmer_genefinder_MF"
    blast_dir = base / "06_analysis" / "blast_genefinder_MF"
    ips_summary = base / "06_analysis" / "interproscan_MF" / "MF_interproscan_summary.tsv"
    blastx_dir = base / "06_analysis" / "blastx_nr_MF_candidates"

    families = ["CYP", "OMT"] if args.family == "both" else [args.family]

    print("=" * 60)
    print("  Annotate candidate evidence (original TSVs are not modified)")
    print("=" * 60)

    for family in families:
        print()
        print(f"### {family}")
        if family == "CYP":
            candidates = cand_dir / "cyp_candidates_MF.tsv"
            out_file = cand_dir / "cyp_candidates_MF_annotated.tsv"
            hmmer_file = hmmer_dir / "cyp450_hmmer_gene_ids.txt"
            blast_file = blast_dir / "cyp450_blast_ids.txt"
            blastx_file = blastx_dir / "cyp450_vs_nr.besthit.tsv"
        else:
            candidates = cand_dir / "omt_candidates_MF.tsv"
            out_file = cand_dir / "omt_candidates_MF_annotated.tsv"
            hmmer_file = hmmer_dir / "omt_hmmer_gene_ids.txt"
            blast_file = blast_dir / "omt_blast_ids.txt"
            blastx_file = blastx_dir / "omt_vs_nr.besthit.tsv"

        hmmer_ids = load_id_set(hmmer_file)
        blast_ids = load_id_set(blast_file)
        print(f"  HMMER genes: {len(hmmer_ids)} from {hmmer_file.name}")
        print(f"  BLAST genefinder genes: {len(blast_ids)} from {blast_file.name}")

        cand_table = pd.read_csv(candidates, sep="\t")
        if "gene_id" not in cand_table.columns:
            cand_table = cand_table.rename(columns={cand_table.columns[0]: "gene_id"})
        wanted = set(cand_table["gene_id"].astype(str))

        ips_by_gene = load_interpro_for_genes(ips_summary, wanted)
        blastx_by_gene = load_blastx_best_per_gene(blastx_file)
        annotate_one(
            candidates, out_file,
            hmmer_ids, blast_ids, ips_by_gene, blastx_by_gene,
        )

    print("Done. Open the *_annotated.tsv files; the originals are unchanged.")
    print("=" * 60)


if __name__ == "__main__":
    main()
