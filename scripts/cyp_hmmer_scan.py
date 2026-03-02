#!/usr/bin/env python3
"""
===============================================================================
STEP 2: CYP450 HMMER DOMAIN SCAN
===============================================================================

Identify CYP450 genes by scanning a protein FASTA file for the P450 Pfam
domain (PF00067) using HMMER's hmmsearch.

Workflow:
  1. Download the PF00067 HMM profile from Pfam/InterPro (if not cached)
  2. Prepare the HMM with hmmpress
  3. Run hmmsearch against the carrot protein FASTA
  4. Parse --tblout results and filter by E-value
  5. Save confirmed P450 gene IDs to CSV

Output:  cyp_hmmer_confirmed.csv
Columns: gene_id, protein_id, hmmer_evalue, hmmer_score, hmmer_bias,
         domain_name, domain_accession

REQUIREMENTS:
  - HMMER3 installed (hmmsearch, hmmpress in PATH)
  - Python packages: pandas

USAGE:
  python scripts/cyp_hmmer_scan.py

  # Override defaults:
  python scripts/cyp_hmmer_scan.py \
      --proteins /path/to/protein.faa \
      --output cyp_hmmer_confirmed.csv \
      --evalue 1e-5

  # Skip download (if you already have the HMM):
  python scripts/cyp_hmmer_scan.py --hmm /path/to/PF00067.hmm

NOTE:
  - Uses hmmsearch (one profile vs many sequences), NOT hmmscan
  - The E-value threshold of 1e-5 is deliberately permissive for discovery;
    tighten to 1e-10 or 1e-20 for higher stringency
  - Gene IDs are extracted from protein FASTA headers by splitting on '|'
    and common delimiters — adjust extract_gene_id() if your headers differ

Author: Daisy Cortes
===============================================================================
"""

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

# ============================================================================
# USER CONFIGURATION — change these to match your file locations
# ============================================================================

PROTEIN_FASTA = "04_reference/GCF_001625215.2_DH1_v3.0_protein.faa"
OUTPUT_FILE = "07_NRdatabase/cyp450_database/cyp_hmmer_confirmed.csv"
EVALUE_CUTOFF = 1e-5

# Where to cache the downloaded HMM file
HMM_CACHE_DIR = "07_NRdatabase/cyp450_database/hmm_profiles"
PF00067_ACCESSION = "PF00067"

# Pfam HMM download URLs (try in order; mirrors change over time)
HMM_URLS = [
    f"https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/{PF00067_ACCESSION}?annotation=hmm",
    f"https://pfam.xfam.org/family/{PF00067_ACCESSION}/hmm",
]


# ============================================================================
# DEPENDENCY CHECKS
# ============================================================================

def check_hmmer():
    """Verify that HMMER3 tools are available on PATH."""
    for tool in ['hmmsearch', 'hmmpress']:
        if shutil.which(tool) is None:
            print(f"ERROR: '{tool}' not found in PATH.")
            print("  Install HMMER3:")
            print("    conda install -c bioconda hmmer")
            print("    # or: brew install hmmer")
            print("    # or: sudo apt install hmmer")
            sys.exit(1)

    result = subprocess.run(['hmmsearch', '-h'], capture_output=True, text=True)
    version_line = result.stdout.split('\n')[1] if result.stdout else 'unknown'
    print(f"  HMMER found: {version_line.strip()}")


# ============================================================================
# DOWNLOAD PF00067 HMM
# ============================================================================

def download_hmm(cache_dir):
    """
    Download the P450 Pfam HMM profile (PF00067).

    Tries multiple URLs since Pfam mirrors change. Falls back to a manual
    download message if all URLs fail.
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    hmm_path = cache_dir / f"{PF00067_ACCESSION}.hmm"

    if hmm_path.exists() and hmm_path.stat().st_size > 100:
        print(f"  Using cached HMM: {hmm_path}")
        return str(hmm_path)

    print(f"  Downloading {PF00067_ACCESSION} HMM profile...")

    # Try wget first, then curl
    for url in HMM_URLS:
        print(f"  Trying: {url}")

        if shutil.which('wget'):
            cmd = ['wget', '-q', '-O', str(hmm_path), url]
        elif shutil.which('curl'):
            cmd = ['curl', '-sL', '-o', str(hmm_path), url]
        else:
            print("  WARNING: Neither wget nor curl found. Cannot download.")
            continue

        try:
            result = subprocess.run(cmd, capture_output=True, text=True,
                                    timeout=120)
            if result.returncode == 0 and hmm_path.exists() and hmm_path.stat().st_size > 100:
                # Validate it looks like an HMM file
                with open(hmm_path) as f:
                    first_line = f.readline()
                if first_line.startswith('HMMER'):
                    print(f"  Downloaded successfully: {hmm_path}")
                    return str(hmm_path)
                else:
                    print(f"  WARNING: Downloaded file does not look like an HMM")
                    hmm_path.unlink(missing_ok=True)
        except subprocess.TimeoutExpired:
            print(f"  Timed out")
            hmm_path.unlink(missing_ok=True)

    # All URLs failed — provide manual instructions
    print()
    print("  COULD NOT AUTO-DOWNLOAD the HMM. Please download manually:")
    print(f"    1. Go to: https://www.ebi.ac.uk/interpro/entry/pfam/{PF00067_ACCESSION}/")
    print(f"    2. Click 'Curation' tab → download the HMM")
    print(f"    3. Save as: {hmm_path}")
    print(f"    4. Re-run this script, or use: --hmm {hmm_path}")
    sys.exit(1)


# ============================================================================
# PREPARE AND RUN HMMER
# ============================================================================

def prepare_hmm(hmm_path):
    """Run hmmpress to prepare the HMM for searching."""
    pressed_files = [hmm_path + ext for ext in ['.h3m', '.h3i', '.h3f', '.h3p']]
    if all(Path(f).exists() for f in pressed_files):
        print(f"  HMM already pressed (index files exist)")
        return

    print(f"  Running hmmpress on {hmm_path}...")
    result = subprocess.run(
        ['hmmpress', '-f', hmm_path],
        capture_output=True, text=True, timeout=60,
    )
    if result.returncode != 0:
        print(f"  ERROR: hmmpress failed:\n{result.stderr.strip()}")
        sys.exit(1)
    print(f"  hmmpress complete")


def run_hmmsearch(hmm_path, protein_fasta, evalue, cpus=4):
    """
    Run hmmsearch: one HMM profile against many protein sequences.

    Returns path to the --tblout (per-sequence table) output file.
    """
    protein_fasta = Path(protein_fasta)
    if not protein_fasta.exists():
        print(f"ERROR: Protein FASTA not found: {protein_fasta}")
        sys.exit(1)

    tblout = tempfile.NamedTemporaryFile(
        mode='w', suffix='_hmmsearch_tbl.txt', delete=False,
    )
    tblout.close()

    domtblout = tempfile.NamedTemporaryFile(
        mode='w', suffix='_hmmsearch_domtbl.txt', delete=False,
    )
    domtblout.close()

    print(f"  Running hmmsearch...")
    print(f"    HMM:      {hmm_path}")
    print(f"    Proteins: {protein_fasta}")
    print(f"    E-value:  {evalue}")

    cmd = [
        'hmmsearch',
        '--tblout', tblout.name,
        '--domtblout', domtblout.name,
        '-E', str(evalue),
        '--cpu', str(cpus),
        '--noali',
        hmm_path,
        str(protein_fasta),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)

    if result.returncode != 0:
        print(f"  ERROR: hmmsearch failed:\n{result.stderr.strip()}")
        Path(tblout.name).unlink(missing_ok=True)
        Path(domtblout.name).unlink(missing_ok=True)
        sys.exit(1)

    print(f"  hmmsearch complete")
    return tblout.name, domtblout.name


# ============================================================================
# PARSE HMMER OUTPUT
# ============================================================================

def extract_gene_id(protein_id):
    """
    Extract a gene-level ID from a protein accession.

    NCBI RefSeq protein headers look like:
      >XP_017313308.1 ...     → gene ID is the LOC number in the GFF
      >NP_001316995.1 ...

    For GCF_001625215.2, protein IDs (XP_*) map to gene IDs (LOC*) via the
    GFF. Since we don't have a direct mapping here, we keep the protein ID
    and let the combine step handle the join.

    If your FASTA headers include gene IDs (e.g., ">LOC108192212|XP_..."),
    split on '|' and take the first part.
    """
    # If header has pipe-delimited fields, first field is often the gene ID
    if '|' in protein_id:
        return protein_id.split('|')[0]
    return protein_id


def parse_hmmer_tblout(tblout_path, evalue_cutoff):
    """
    Parse hmmsearch --tblout file.

    Tblout format (space-delimited, # lines are comments):
      target_name  target_acc  query_name  query_acc  E-value  score  bias  ...
    """
    print(f"  Parsing hmmsearch results...")
    rows = []

    with open(tblout_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 9:
                continue

            target_name = parts[0]   # protein ID (e.g., XP_017313308.1)
            target_acc  = parts[1]   # "-" or accession
            query_name  = parts[2]   # HMM name (e.g., p450)
            query_acc   = parts[3]   # PF00067.25
            evalue      = float(parts[4])
            score       = float(parts[5])
            bias        = float(parts[6])

            if evalue > evalue_cutoff:
                continue

            gene_id = extract_gene_id(target_name)

            rows.append({
                'gene_id':            gene_id,
                'protein_id':         target_name,
                'hmmer_evalue':       evalue,
                'hmmer_score':        score,
                'hmmer_bias':         bias,
                'domain_name':        query_name,
                'domain_accession':   query_acc,
            })

    df = pd.DataFrame(rows)

    if df.empty:
        print(f"  WARNING: No hits passed the E-value cutoff ({evalue_cutoff})")
        return df

    # Keep best hit per gene (lowest E-value)
    df = df.sort_values('hmmer_evalue')
    df = df.drop_duplicates(subset='gene_id', keep='first')
    df = df.sort_values('hmmer_score', ascending=False).reset_index(drop=True)

    print(f"  Hits passing E-value < {evalue_cutoff}: {len(df)}")
    return df


# ============================================================================
# OUTPUT AND SUMMARY
# ============================================================================

def save_results(df, output_path):
    """Save HMMER-confirmed P450 genes to CSV."""
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    if df.empty:
        print("\nWARNING: No HMMER-confirmed P450 genes found!")
        pd.DataFrame(columns=[
            'gene_id', 'protein_id', 'hmmer_evalue', 'hmmer_score',
            'hmmer_bias', 'domain_name', 'domain_accession',
        ]).to_csv(output_path, index=False)
        return

    df.to_csv(output_path, index=False)
    print(f"\nSaved: {output_path}")
    print(f"  Rows: {len(df)}")
    print(f"  Columns: {', '.join(df.columns)}")


def print_summary(df):
    """Print HMMER scan summary."""
    if df.empty:
        return

    print()
    print("=" * 60)
    print("STEP 2 SUMMARY — HMMER DOMAIN SCAN")
    print("=" * 60)
    print(f"  Confirmed P450 domain proteins: {len(df)}")
    print(f"  Unique gene IDs:                {df['gene_id'].nunique()}")

    print(f"\n  Score distribution:")
    print(f"    Min score:  {df['hmmer_score'].min():.1f}")
    print(f"    Max score:  {df['hmmer_score'].max():.1f}")
    print(f"    Mean score: {df['hmmer_score'].mean():.1f}")

    print(f"\n  E-value distribution:")
    print(f"    Best (lowest):  {df['hmmer_evalue'].min():.2e}")
    print(f"    Worst (highest): {df['hmmer_evalue'].max():.2e}")

    print(f"\n  Top 5 hits (by score):")
    for _, row in df.head(5).iterrows():
        print(f"    {row['gene_id']:25s}  score={row['hmmer_score']:.1f}  "
              f"E={row['hmmer_evalue']:.2e}")

    print("=" * 60)


# ============================================================================
# CLI
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Step 2: HMMER domain scan for P450 (PF00067)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        '--proteins', default=PROTEIN_FASTA,
        help=f'Protein FASTA file (default: {PROTEIN_FASTA})',
    )
    parser.add_argument(
        '--hmm', default=None,
        help='Path to pre-downloaded PF00067 HMM (skips download)',
    )
    parser.add_argument(
        '--output', '-o', default=OUTPUT_FILE,
        help=f'Output CSV path (default: {OUTPUT_FILE})',
    )
    parser.add_argument(
        '--evalue', type=float, default=EVALUE_CUTOFF,
        help=f'E-value cutoff (default: {EVALUE_CUTOFF})',
    )
    parser.add_argument(
        '--cpus', type=int, default=4,
        help='Number of CPUs for hmmsearch (default: 4)',
    )
    args = parser.parse_args()

    print("=" * 60)
    print("CYP450 DATABASE — STEP 2: HMMER DOMAIN SCAN")
    print("=" * 60)
    print(f"  Proteins: {args.proteins}")
    print(f"  E-value:  {args.evalue}")
    print(f"  Output:   {args.output}")
    print()

    # Check dependencies
    print("Checking dependencies...")
    check_hmmer()
    print()

    # Get HMM file
    if args.hmm:
        hmm_path = args.hmm
        if not Path(hmm_path).exists():
            print(f"ERROR: HMM file not found: {hmm_path}")
            sys.exit(1)
        print(f"Using provided HMM: {hmm_path}")
    else:
        print("Obtaining PF00067 HMM profile...")
        hmm_path = download_hmm(HMM_CACHE_DIR)

    print()

    # Prepare HMM
    prepare_hmm(hmm_path)
    print()

    # Run hmmsearch
    tblout_path, domtblout_path = run_hmmsearch(
        hmm_path, args.proteins, args.evalue, args.cpus,
    )

    # Parse results
    df = parse_hmmer_tblout(tblout_path, args.evalue)

    # Clean up temp files
    Path(tblout_path).unlink(missing_ok=True)
    Path(domtblout_path).unlink(missing_ok=True)

    # Save and summarize
    save_results(df, args.output)
    print_summary(df)

    return 0


if __name__ == '__main__':
    sys.exit(main())
