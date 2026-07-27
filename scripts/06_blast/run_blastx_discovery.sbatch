#!/bin/bash
################################################################################
# Run BLASTX against a protein database (ALL nucleotide genes / transcripts)
################################################################################
#
# WHAT THIS DOES:
#   Annotates EVERY sequence in your query FASTA (not just CYP/OMT).
#   BLASTX = nucleotide query vs protein database.
#
# FOR TRINITY / CD-HIT (nutmeg MF — recommended after run_cdhit.sbatch):
#   Default query = 01_processed/00_6_cdhit/MF_trinity_cdhit95.fasta
#   That is ALL CD-HIT-collapsed Trinity transcripts.
#
# FOR REFERENCE-GENOME CDS (carrot DC/DG):
#   Default query = 06_analysis/blast_input_<SP>/all_genes_cds.fasta
#
# DATABASES: swissprot (start here — smaller/faster), nr, refseq_protein, DCdb
#   Prefer swissprot first. nr can take weeks even with 14 days + resume.
#
# WALLTIME: 14 days on ARC normal_q long QoS (maximum non-preemptable).
# RESUME: if the job times out, resubmit the SAME command — unfinished
#   sequences are detected and continued automatically.
#
# How to run (or resume):
#   # ALL Trinity/CD-HIT transcripts vs swissprot (recommended start)
#   sbatch scripts/06_blast/run_blastx_discovery.sbatch MF swissprot
#
#   # Same, but against the huge nr database (very slow)
#   sbatch scripts/06_blast/run_blastx_discovery.sbatch MF nr
#
#   # Optional: pass an explicit FASTA path as arg 3
#   sbatch scripts/06_blast/run_blastx_discovery.sbatch MF swissprot \
#       /projects/tholl_lab_1/daisy_analysis/01_processed/00_6_cdhit/MF_trinity_cdhit95.fasta
#
#   # Carrot reference CDS (DC/DG)
#   sbatch scripts/06_blast/run_blastx_discovery.sbatch DC swissprot
#
# FASTER ALTERNATIVE (if TransDecoder .pep already exists):
#   BLASTp is 3–6x faster than BLASTX. Use:
#     sbatch scripts/06_blast/run_blastp_discovery.sbatch MF swissprot \
#         /projects/tholl_lab_1/daisy_analysis/01_processed/00_6_cdhit/MF_trinity_cdhit95.fasta.transdecoder.pep
################################################################################

#SBATCH --job-name=blastx_discoveryfilter
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
# Required for 14-day walltime on ARC (default QoS max is only 7 days)
#SBATCH --qos=tc_normal_long
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=14-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu

#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/06_analysis
#SBATCH -o blastx_discoveryfilter_%j.out
#SBATCH -e blastx_discoveryfilter_%j.err


# Strict error handling (pipefail but not -e, so resume logic works on timeout)
set -uo pipefail

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq

# ===== Load shared config =====
source "${SLURM_SUBMIT_DIR:-.}/scripts/config.sh"

export BLASTDB="${DB_DIR}/blastdb"

# ============================================================================
# USER-EDITABLE VARIABLES - DISCOVERY FILTER MODE
# ============================================================================
# Usage: sbatch run_blastx_discovery.sbatch [SPECIES] [DATABASE] [QUERY_FASTA]
#   Examples:
#     sbatch ... MF swissprot          # ALL CD-HIT Trinity transcripts
#     sbatch ... MF nr                 # same, vs nr (very slow)
#     sbatch ... DC swissprot          # carrot CDS from GTF extract
#     sbatch ... MF swissprot /path/to/custom.fasta
SPECIES="${1:-MF}"          # MF (Trinity), DC, or DG
DBNAME="${2:-swissprot}"    # swissprot (recommended), nr, refseq_protein, DCdb
CUSTOM_QUERY="${3:-}"       # optional: explicit nucleotide FASTA

# Validate species
if [[ ! "$SPECIES" =~ ^(DC|DG|MF)$ ]]; then
    echo "ERROR: Invalid species '$SPECIES'. Must be DC, DG, or MF"
    echo "Usage: sbatch scripts/06_blast/run_blastx_discovery.sbatch [SPECIES] [DATABASE] [QUERY_FASTA]"
    echo "  SPECIES:  DC, DG, or MF"
    echo "  DATABASE: swissprot, nr, refseq_protein, or DCdb"
    exit 1
fi

# Validate database
if [[ ! "$DBNAME" =~ ^(nr|swissprot|refseq_protein|DCdb)$ ]]; then
    echo "ERROR: Invalid database '$DBNAME'. Must be nr, swissprot, refseq_protein, or DCdb"
    echo "Usage: sbatch scripts/06_blast/run_blastx_discovery.sbatch [SPECIES] [DATABASE] [QUERY_FASTA]"
    echo "  SPECIES:  DC, DG, or MF"
    echo "  DATABASE: swissprot, nr, refseq_protein, or DCdb"
    exit 1
fi

# ---- Pick the query FASTA (ALL genes / transcripts) ----
# Priority:
#   1) explicit path as arg 3
#   2) MF (Trinity): CD-HIT assembly
#   3) DC/DG (reference): extracted CDS from GTF
#   4) fallback to whichever of those files exists
CDHIT_FASTA="${PROCESSED_DIR}/00_6_cdhit/${SPECIES}_trinity_cdhit95.fasta"
CDS_FASTA="${ANALYSIS_DIR}/blast_input_${SPECIES}/all_genes_cds.fasta"

if [ -n "$CUSTOM_QUERY" ]; then
    QUERY_FASTA="$CUSTOM_QUERY"
    QUERY_BASENAME=$(basename "$QUERY_FASTA")
    QUERY_BASENAME="${QUERY_BASENAME%.fasta}"
    QUERY_BASENAME="${QUERY_BASENAME%.fa}"
    OUTDIR="${ANALYSIS_DIR}/blastx_${SPECIES}_${QUERY_BASENAME}"
elif [[ "$SPECIES" == "MF" ]] && [ -f "$CDHIT_FASTA" ]; then
    # ALL CD-HIT Trinity transcripts (nucleotide)
    QUERY_FASTA="$CDHIT_FASTA"
    OUTDIR="${ANALYSIS_DIR}/blastx_${SPECIES}_cdhit"
elif [ -f "$CDS_FASTA" ]; then
    QUERY_FASTA="$CDS_FASTA"
    OUTDIR="${ANALYSIS_DIR}/blastx_${SPECIES}"
elif [ -f "$CDHIT_FASTA" ]; then
    QUERY_FASTA="$CDHIT_FASTA"
    OUTDIR="${ANALYSIS_DIR}/blastx_${SPECIES}_cdhit"
else
    echo "ERROR: No query FASTA found for species ${SPECIES}."
    echo "Looked for:"
    echo "  CD-HIT Trinity: ${CDHIT_FASTA}"
    echo "  Extracted CDS:  ${CDS_FASTA}"
    echo ""
    echo "For Trinity/nutmeg, run CD-HIT first:"
    echo "  sbatch scripts/04_busco/run_cdhit.sbatch ${SPECIES}"
    echo "For carrot reference CDS, run:"
    echo "  sbatch scripts/06_blast/run_extract_cds.sbatch ${SPECIES}"
    exit 1
fi

# Validate input file exists
if [ ! -f "$QUERY_FASTA" ]; then
    echo "ERROR: Query file not found: $QUERY_FASTA"
    exit 1
fi

# BLAST filter parameters
EVALUE=1e-4              # DISCOVERY MODE: permissive e-value
QCOV_HSP_PERC=40         # DISCOVERY MODE: low coverage threshold
MAX_TARGET_SEQS=25       # DISCOVERY MODE: more hits
# ============================================================================

# Create output directory and move into it
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# Print reproducibility information
echo "=========================================="
echo "BLASTX Analysis Started (DISCOVERY MODE)"
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
python3 -c "from rna_pipeline.tools.blast import check_blast; assert check_blast(), 'blastx not found — run: conda install -c bioconda blast'"
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
OUTPUT_FILE="${OUTDIR}/blastx_${SPECIES}_${DBNAME}_discovery.tsv"
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

# Run BLASTX using Python with RESUME support
python - <<'PY'
import os
import subprocess
import sys
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
output_file = f"blastx_{species}_{dbname}_discovery.tsv"

# ── RESUME LOGIC ──────────────────────────────────────────────────────────
completed_ids = set()
actual_query = query_fasta

if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
    with open(output_file) as fh:
        for line in fh:
            if line.strip():
                completed_ids.add(line.split('\t')[0])

    print(f"\n>>> RESUME MODE: Found existing {output_file}")
    print(f">>> {len(completed_ids)} query sequences already have results")

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

# ── RUN BLASTX ────────────────────────────────────────────────────────────
if completed_ids:
    tmp_output = os.path.join(outdir, f"_tmp_{species}_{dbname}_discovery.tsv")
else:
    tmp_output = output_file

blastx_cmd = [
    'blastx',
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
print("Running BLASTX command:")
print("========================================")
print(' '.join(blastx_cmd))
print("========================================\n")

try:
    subprocess.run(blastx_cmd, check=True)
except subprocess.CalledProcessError as e:
    print(f"ERROR: BLASTX failed with exit code {e.returncode}", file=sys.stderr)
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
    if os.path.exists(actual_query):
        os.remove(actual_query)
        print(">>> Cleaned up temporary FASTA file")

# ── SUMMARY ───────────────────────────────────────────────────────────────
print("\n========================================")
print(f"Done: BLASTX results written to {output_file}")
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

print("\nBLASTX analysis complete!")
PY

echo "=========================================="
echo "BLASTX Analysis Completed (DISCOVERY MODE)"
date
echo "=========================================="