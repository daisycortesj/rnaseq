#!/bin/bash
################################################################################
# Run BLASTp against BLAST database (protein query from .fasta protein file)
################################################################################
#
# Uses the BLAST database (swissprot, nr, refseq_protein). 
# Database is downloaded in 07_NRdatabase directory and set as the BLASTDB environment variable.
# Database created with download_blastdb.sbatch script.
# Uses database for discovery mode: more hits and permissive e-value.
#
# BLASTp (protein vs protein) is 3-6x FASTER than BLASTX (nucleotide vs protein)
# Use this AFTER translating CDS to protein with run_translate_cds.sbatch
#
# RESUME SUPPORT:
#   If a previous run timed out, just resubmit the SAME command.
#   The script detects the existing partial output TSV and only BLASTs
#   sequences that don't have results yet, then appends new results.
#   You can resubmit as many times as needed until all sequences finish.
#
# How to run (or resume):
#   sbatch scripts/blastp_discoveryfilter.sbatch DC nr
#   sbatch scripts/blastp_discoveryfilter.sbatch DC swissprot
#   sbatch scripts/blastp_discoveryfilter.sbatch DC refseq_protein
#   sbatch scripts/blastp_discoveryfilter.sbatch DC DCdb
#   sbatch scripts/blastp_discoveryfilter.sbatch DG nr
#   sbatch scripts/blastp_discoveryfilter.sbatch DG swissprot
#   sbatch scripts/blastp_discoveryfilter.sbatch DG refseq_protein
#   sbatch scripts/blastp_discoveryfilter.sbatch DG DCdb
#   sbatch scripts/blastp_discoveryfilter.sbatch MF nr
#   sbatch scripts/blastp_discoveryfilter.sbatch MF swissprot
#   sbatch scripts/blastp_discoveryfilter.sbatch MF refseq_protein
#   sbatch scripts/blastp_discoveryfilter.sbatch MF DCdb
#
#   DCdb is a small custom DB (~50-100K seqs) extracted from nr.
#   Build it first: sbatch scripts/extract_DCdb_db.sbatch
#   TIP: Override resources for the tiny DB:
#     sbatch --mem=8G --cpus-per-task=4 --time=1:00:00 scripts/blastp_discoveryfilter.sbatch DC DCdb
################################################################################

#SBATCH --job-name=blastp_discovery
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=3-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/06_analysis
#SBATCH -o blastp_discoveryfilter_%j.out
#SBATCH -e blastp_discoveryfilter_%j.err


# Strict error handling (pipefail but not -e, so resume logic works on timeout)
set -uo pipefail

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq

# Set BLAST database location
export BLASTDB=/projects/tholl_lab_1/daisy_analysis/07_NRdatabase/blastdb

# ============================================================================
# USER-EDITABLE VARIABLES - DISCOVERY FILTER MODE
# ============================================================================
# Command-line arguments with defaults
# Usage: sbatch blastp_discoveryfilter.sbatch [SPECIES] [DATABASE]
#   Examples:
#     sbatch blastp_discoveryfilter.sbatch                    # DC + nr (defaults)
#     sbatch blastp_discoveryfilter.sbatch DC                 # DC + nr
#     sbatch blastp_discoveryfilter.sbatch DG swissprot       # DG + swissprot
#     sbatch blastp_discoveryfilter.sbatch MF refseq_protein  # MF + refseq_protein
#     sbatch blastp_discoveryfilter.sbatch DC DCdb   # DC + DCdb (fast)
SPECIES="${1:-DC}"     # DC (Daucus carota), DG (Daucus glochidiatus), MF (Melaleuca)
DBNAME="${2:-nr}"      # nr, swissprot, refseq_protein, DCdb

# Validate species
if [[ ! "$SPECIES" =~ ^(DC|DG|MF)$ ]]; then
    echo "ERROR: Invalid species '$SPECIES'. Must be DC, DG, or MF"
    echo "Usage: sbatch blastp_discoveryfilter.sbatch [SPECIES] [DATABASE]"
    echo "  SPECIES:  DC, DG, or MF"
    echo "  DATABASE: nr, swissprot, refseq_protein, or DCdb"
    exit 1
fi

# Validate database
if [[ ! "$DBNAME" =~ ^(nr|swissprot|refseq_protein|DCdb)$ ]]; then
    echo "ERROR: Invalid database '$DBNAME'. Must be nr, swissprot, refseq_protein, or DCdb"
    echo "Usage: sbatch blastp_discoveryfilter.sbatch [SPECIES] [DATABASE]"
    echo "  SPECIES:  DC, DG, or MF"
    echo "  DATABASE: nr, swissprot, refseq_protein, or DCdb"
    exit 1
fi

# Auto-configure paths based on species
BASE_DIR=/projects/tholl_lab_1/daisy_analysis
QUERY_FASTA="${BASE_DIR}/06_analysis/blast_input_${SPECIES}/all_genes_protein.fasta"
OUTDIR="${BASE_DIR}/06_analysis/blastp_${SPECIES}"

# Validate input file exists
if [ ! -f "$QUERY_FASTA" ]; then
    echo "ERROR: Query file not found: $QUERY_FASTA"
    echo "Make sure to run: sbatch scripts/run_translate_cds.sbatch $SPECIES"
    exit 1
fi

# BLAST filter parameters
EVALUE=1e-4              # DISCOVERY MODE: permissive e-value
QCOV_HSP_PERC=40         # DISCOVERY MODE: low coverage threshold
MAX_TARGET_SEQS=25      # DISCOVERY MODE: more hits
# ============================================================================

# Create output directory and move into it
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# Print reproducibility information
echo "=========================================="
echo "BLASTp Analysis Started (DISCOVERY MODE)"
echo "=========================================="
echo "Species:  $SPECIES"
echo "Query:    $QUERY_FASTA"
echo "Database: $DBNAME"
echo "Output:   $OUTDIR"
echo "------------------------------------------"
date
echo "Hostname: $(hostname)"
echo "BLASTDB: $BLASTDB"
echo "PATH: $PATH"
echo "------------------------------------------"
which blastp
blastp -version
which blastdbcmd
echo "=========================================="

# Export variables for Python to read
export QUERY_FASTA
export DBNAME
export OUTDIR
export EVALUE
export QCOV_HSP_PERC
export MAX_TARGET_SEQS
export SPECIES

# Trap SIGTERM (sent by SLURM on timeout) to merge any partial temp results
OUTPUT_FILE="${OUTDIR}/blastp_${SPECIES}_${DBNAME}_discovery.tsv"
TMP_OUTPUT="${OUTDIR}/_tmp_${SPECIES}_${DBNAME}_discovery.tsv"
trap 'echo ">>> SIGTERM received (likely SLURM timeout)"; \
      if [ -f "$TMP_OUTPUT" ] && [ -s "$TMP_OUTPUT" ]; then \
          echo ">>> Merging partial temp results into main output for next resume..."; \
          cat "$TMP_OUTPUT" >> "$OUTPUT_FILE"; \
          rm -f "$TMP_OUTPUT"; \
          echo ">>> Done. Resubmit the same command to continue."; \
      fi; \
      echo ">>> Partial results preserved. Resubmit to resume."; \
      exit 0' SIGTERM

# Run BLASTp using Python with RESUME support
python - <<'PY'
import os
import subprocess
import sys
import tempfile
from pathlib import Path

# Read configuration from environment variables
query_fasta = os.environ['QUERY_FASTA']
dbname = os.environ['DBNAME']
outdir = os.environ['OUTDIR']
evalue = os.environ['EVALUE']
qcov_hsp_perc = os.environ['QCOV_HSP_PERC']
max_target_seqs = os.environ['MAX_TARGET_SEQS']
species = os.environ['SPECIES']
num_threads = os.environ.get('SLURM_CPUS_PER_TASK', '1')

# Validate input file exists
if not os.path.exists(query_fasta):
    print(f"ERROR: Query file not found: {query_fasta}", file=sys.stderr)
    sys.exit(1)

# Define output file
output_file = f"blastp_{species}_{dbname}_discovery.tsv"

# ── RESUME LOGIC ──────────────────────────────────────────────────────────
# If a partial output file exists from a previous timed-out run,
# figure out which query sequences already have results and skip them.
completed_ids = set()
actual_query = query_fasta

if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
    with open(output_file) as fh:
        for line in fh:
            if line.strip():
                completed_ids.add(line.split('\t')[0])

    print(f"\n>>> RESUME MODE: Found existing {output_file}")
    print(f">>> {len(completed_ids)} query sequences already have results")

    # Count total sequences in the original FASTA
    total_seqs = 0
    with open(query_fasta) as fh:
        for line in fh:
            if line.startswith('>'):
                total_seqs += 1
    print(f">>> Total sequences in query FASTA: {total_seqs}")

    remaining = total_seqs - len(completed_ids)
    if remaining <= 0:
        print(">>> All sequences already processed! Nothing to do.")
        print(f"\nTotal results in {output_file}: {len(completed_ids)} queries")
        sys.exit(0)

    print(f">>> Sequences remaining: {remaining}")

    # Write a filtered FASTA with only unfinished sequences
    filtered_fasta = os.path.join(outdir, f"_remaining_{species}_{dbname}.fasta")
    written = 0
    skip = False
    with open(query_fasta) as fin, open(filtered_fasta, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                seq_id = line[1:].split()[0]
                skip = seq_id in completed_ids
                if not skip:
                    fout.write(line)
                    written += 1
            else:
                if not skip:
                    fout.write(line)

    print(f">>> Wrote {written} remaining sequences to temporary FASTA")
    print(f">>> New results will be APPENDED to {output_file}")
    actual_query = filtered_fasta
else:
    print(f"\n>>> Starting fresh (no existing {output_file} found)")

# ── RUN BLASTp ────────────────────────────────────────────────────────────
# When resuming, write to a temp file and then append
if completed_ids:
    tmp_output = os.path.join(outdir, f"_tmp_{species}_{dbname}_discovery.tsv")
else:
    tmp_output = output_file

blastp_cmd = [
    'blastp',
    '-query', actual_query,
    '-db', dbname,
    '-evalue', evalue,
    '-max_target_seqs', max_target_seqs,
    '-qcov_hsp_perc', qcov_hsp_perc,
    '-max_hsps', '1',
    '-num_threads', num_threads,
    '-outfmt', '6 qseqid sseqid stitle sscinames pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp',
    '-out', tmp_output
]

print("\n========================================")
print("Running BLASTp command:")
print("========================================")
print(' '.join(blastp_cmd))
print("========================================\n")

try:
    subprocess.run(blastp_cmd, check=True)
except subprocess.CalledProcessError as e:
    print(f"ERROR: BLASTp failed with exit code {e.returncode}", file=sys.stderr)
    # Even on failure, append any partial results so next resume picks them up
    if completed_ids and os.path.exists(tmp_output) and os.path.getsize(tmp_output) > 0:
        print("Appending partial new results to main output before exiting...")
        with open(tmp_output) as fin, open(output_file, 'a') as fout:
            fout.write(fin.read())
        os.remove(tmp_output)
    sys.exit(1)

# ── MERGE RESULTS ─────────────────────────────────────────────────────────
if completed_ids:
    if os.path.exists(tmp_output) and os.path.getsize(tmp_output) > 0:
        with open(tmp_output) as fin, open(output_file, 'a') as fout:
            fout.write(fin.read())
        os.remove(tmp_output)
        print(f">>> Appended new results to {output_file}")
    # Clean up temporary filtered FASTA
    if os.path.exists(actual_query):
        os.remove(actual_query)
        print(">>> Cleaned up temporary FASTA file")

# ── SUMMARY ───────────────────────────────────────────────────────────────
print("\n========================================")
print(f"Done: BLASTp results written to {output_file}")
print("========================================")

output_path = Path(output_file)
if output_path.exists():
    total_lines = sum(1 for _ in open(output_file))
    unique_queries = len(set(line.split('\t')[0] for line in open(output_file) if line.strip()))
    print(f"\nTotal result lines: {total_lines}")
    print(f"Unique query sequences with hits: {unique_queries}")
    print("\nFirst 5 lines of output:")
    print("------------------------------------------")
    with open(output_path, 'r') as f:
        for i, line in enumerate(f):
            if i >= 5:
                break
            print(line.rstrip())
    print("------------------------------------------")
else:
    print("WARNING: Output file was not created", file=sys.stderr)

print("\nBLASTp analysis complete!")
PY

echo "=========================================="
echo "BLASTp Analysis Completed (DISCOVERY MODE)"
date
echo "=========================================="
