#!/usr/bin/env python3
"""
Fetch protein sequences from NCBI using edirect (efetch).

Reads a TSV that contains NCBI protein accessions (e.g. PyDESeq2
with sseqid, or any TSV with an accession column) and downloads protein FASTA
via NCBI E-utilities.

Requires: conda env with entrez-direct (edirect) and pandas.
  conda activate rnaseq   # or your env with entrez-direct

Usage:
  # From combined PyDESeq2 TSV (uses sseqid column)
  python scripts/fetch_proteins_edirect.py 06_analysis/combined_candidates_with_blast.tsv -o proteins.fasta

  # Specify accession column and optional gene ID column (for FASTA headers)
  python scripts/fetch_proteins_edirect.py results.tsv --acc-column sseqid --gene-column gene_id -o proteins.fasta

  # Rate limit (default 0.34 s between requests to respect NCBI)
  python scripts/fetch_proteins_edirect.py results.tsv -o proteins.fasta --delay 0.5
"""

import argparse
import subprocess
import sys
import time
from pathlib import Path
from shutil import which

import pandas as pd


def normalize_accession(sseqid):
    """
    Extract clean NCBI accession from BLAST sseqid.
    Examples: 'ref|NP_123.1|' -> 'NP_123.1'; 'XP_012345.1' -> 'XP_012345.1'
    """
    if pd.isna(sseqid) or not str(sseqid).strip():
        return None
    s = str(sseqid).strip()
    if "|" in s:
        parts = s.split("|")
        # ref|ACCESSION| or gb|ACCESSION|
        if len(parts) >= 2:
            return parts[1]
        return parts[0]
    return s


def efetch_protein(accession, email=None):
    """Fetch one protein in FASTA format using edirect efetch. Returns bytes or None."""
    cmd = ["efetch", "-db", "protein", "-id", accession, "-format", "fasta"]
    env = None
    if email:
        env = {**__import__("os").environ, "EMAIL": email}
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=False,
            timeout=30,
            env=env,
        )
        if result.returncode == 0 and result.stdout.strip():
            return result.stdout
        return None
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return None


def main():
    parser = argparse.ArgumentParser(
        description="Fetch protein sequences from NCBI using edirect (TSV with accessions)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "tsv",
        help="TSV file (e.g. combined BLAST+PyDESeq2 with sseqid, or any TSV with accession column)",
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output FASTA file",
    )
    parser.add_argument(
        "--acc-column",
        default="sseqid",
        help="Column containing NCBI protein accessions (default: sseqid)",
    )
    parser.add_argument(
        "--gene-column",
        default="gene_id",
        help="Optional column for gene ID to use in FASTA header (default: gene_id)",
    )
    parser.add_argument(
        "--delay",
        type=float,
        default=0.34,
        help="Seconds between NCBI requests (default: 0.34)",
    )
    parser.add_argument(
        "--email",
        default="",
        help="Email for NCBI (optional; set EMAIL env or use this)",
    )
    parser.add_argument(
        "--max",
        type=int,
        default=5,
        help="Max number of accessions to fetch (0 = all; use for testing)",
    )
    args = parser.parse_args()

    if not which("efetch"):
        print("ERROR: 'efetch' not found. Activate conda env with entrez-direct:", file=sys.stderr)
        print("  conda activate rnaseq", file=sys.stderr)
        return 1

    tsv_path = Path(args.tsv)
    if not tsv_path.exists():
        print(f"ERROR: TSV not found: {tsv_path}", file=sys.stderr)
        return 1

    # Load TSV
    try:
        df = pd.read_csv(tsv_path, sep="\t")
    except Exception as e:
        print(f"ERROR: Could not read TSV: {e}", file=sys.stderr)
        return 1

    if args.acc_column not in df.columns:
        print(f"ERROR: Column '{args.acc_column}' not in TSV. Columns: {list(df.columns)}", file=sys.stderr)
        return 1

    # Unique accessions, keep first occurrence for optional gene_id
    acc_to_gene = {}
    if args.gene_column and args.gene_column in df.columns:
        for _, row in df[[args.acc_column, args.gene_column]].drop_duplicates(args.acc_column).iterrows():
            acc = normalize_accession(row[args.acc_column])
            if acc:
                acc_to_gene[acc] = row[args.gene_column]

    accessions = [normalize_accession(x) for x in df[args.acc_column].dropna().unique()]
    accessions = [a for a in accessions if a]
    accessions = list(dict.fromkeys(accessions))  # preserve order, unique

    if not accessions:
        print("ERROR: No valid accessions found in column", args.acc_column, file=sys.stderr)
        return 1

    if args.max > 0:
        accessions = accessions[: args.max]
        print(f"Limiting to first {args.max} accessions (--max)")

    print(f"TSV: {tsv_path}")
    print(f"Accession column: {args.acc_column}")
    print(f"Unique accessions: {len(accessions)}")
    print(f"Output: {args.output}")
    print()

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    n_ok = 0
    n_fail = 0
    email = args.email or __import__("os").environ.get("EMAIL", "")

    with open(out_path, "wb") as out:
        for i, acc in enumerate(accessions):
            if (i + 1) % 50 == 0 or i == 0:
                print(f"  Fetching {i+1}/{len(accessions)}: {acc} ...")
            fasta = efetch_protein(acc, email=email or None)
            if fasta:
                # Optionally prepend gene_id to header line for traceability
                if acc in acc_to_gene and args.gene_column and args.gene_column in df.columns:
                    gene = acc_to_gene[acc]
                    lines = fasta.split(b"\n")
                    if lines and lines[0].startswith(b">"):
                        lines[0] = lines[0] + f" | gene_id={gene}".encode()
                    fasta = b"\n".join(lines)
                out.write(fasta)
                if not fasta.endswith(b"\n"):
                    out.write(b"\n")
                n_ok += 1
            else:
                n_fail += 1
                print(f"  WARNING: No sequence for {acc}", file=sys.stderr)
            time.sleep(args.delay)

    print()
    print(f"Done. Fetched: {n_ok}, failed: {n_fail}. Output: {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
