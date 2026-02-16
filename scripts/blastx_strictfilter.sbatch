#!/bin/bash
################################################################################
# Run BLASTX against NCBI database (nucleotide query from .fasta CDS file)
# STRICT FILTER MODE: High-confidence hits only
################################################################################
#
# 
# Uses the Blast database (swissprot, nr, refseq_protein). 
# Database is downloaded in 07_NRdatabase directory and set as the BLASTDB environment variable.
# Database created with download_blastdb.sbatch script.
# Uses database for discovery mode: more hits and permissive e-value.
#   - Default: runs against NCBI database via -remote (no local database).
#   - Optional: set DBNAME to a local database path if you have one
#
# How to run:
#   sbatch scripts/blastx_strictfilter.sbatch DC nr
#   sbatch scripts/blastx_strictfilter.sbatch DC swissprot
#   sbatch scripts/blastx_strictfilter.sbatch DC refseq_protein
#   sbatch scripts/blastx_strictfilter.sbatch DG nr
#   sbatch scripts/blastx_strictfilter.sbatch DG swissprot
#   sbatch scripts/blastx_strictfilter.sbatch DG refseq_protein
#   sbatch scripts/blastx_strictfilter.sbatch MF nr
#   sbatch scripts/blastx_strictfilter.sbatch MF swissprot
#   sbatch scripts/blastx_strictfilter.sbatch MF refseq_protein
#
################################################################################

#SBATCH --job-name=blastx_strict
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu

#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/06_analysis/blast_input_<species>
#SBATCH -o blastx_discoveryfilter_<species>_%j.out
#SBATCH -e blastx_discoveryfilter_<species>_%j.err

# Strict error handling
set -euo pipefail

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq

# Set BLAST database location
export BLASTDB=/projects/tholl_lab_1/daisy_analysis/07_NRdatabase/blastdb

# ============================================================================
# USER-EDITABLE VARIABLES - STRICT FILTER MODE
# ============================================================================
# Command-line arguments with defaults
# Usage: sbatch blastx_strictfilter.sbatch [SPECIES] [DATABASE]
#   Examples:
#     sbatch blastx_strictfilter.sbatch                    # DC + nr (defaults)
#     sbatch blastx_strictfilter.sbatch DC                 # DC + nr
#     sbatch blastx_strictfilter.sbatch DG swissprot       # DG + swissprot
#     sbatch blastx_strictfilter.sbatch MF refseq_protein  # MF + refseq_protein
SPECIES="${1:-DC}"     # DC (Daucus carota), DG (Daucus glochidiatus), MF (Melaleuca)
DBNAME="${2:-nr}"      # nr, swissprot, refseq_protein

# Validate species
if [[ ! "$SPECIES" =~ ^(DC|DG|MF)$ ]]; then
    echo "ERROR: Invalid species '$SPECIES'. Must be DC, DG, or MF"
    echo "Usage: sbatch blastx_strictfilter.sbatch [SPECIES] [DATABASE]"
    echo "  SPECIES:  DC, DG, or MF"
    echo "  DATABASE: nr, swissprot, or refseq_protein"
    exit 1
fi

# Validate database
if [[ ! "$DBNAME" =~ ^(nr|swissprot|refseq_protein)$ ]]; then
    echo "ERROR: Invalid database '$DBNAME'. Must be nr, swissprot, or refseq_protein"
    echo "Usage: sbatch blastx_strictfilter.sbatch [SPECIES] [DATABASE]"
    echo "  SPECIES:  DC, DG, or MF"
    echo "  DATABASE: nr, swissprot, or refseq_protein"
    exit 1
fi

# Auto-configure paths based on species
BASE_DIR=/projects/tholl_lab_1/daisy_analysis
QUERY_FASTA="${BASE_DIR}/06_analysis/blast_input_${SPECIES}/all_genes_cds.fasta"
OUTDIR="${BASE_DIR}/06_analysis/blastx_${SPECIES}"

# Validate input file exists
if [ ! -f "$QUERY_FASTA" ]; then
    echo "ERROR: Query file not found: $QUERY_FASTA"
    echo "Make sure to run: python scripts/extract_cds_from_gtf.py $SPECIES"
    exit 1
fi

# BLAST filter parameters
EVALUE=1e-8             # STRICT MODE: very stringent e-value (high confidence)
QCOV_HSP_PERC=70         # STRICT MODE: high coverage requirement
MAX_TARGET_SEQS=10        # STRICT MODE: only top hits
# ============================================================================

# Create output directory and move into it
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# Print reproducibility information
echo "=========================================="
echo "BLASTX Analysis Started (STRICT FILTER MODE)"
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
which blastx
blastx -version
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

# Run BLASTX using Python for easy modification
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
output_file = f"blastx_{species}_{dbname}_strict.tsv"

# Build BLASTX command as a list (easy to modify)
blastx_cmd = [
    'blastx',
    '-query', query_fasta,
    '-db', dbname,
    '-evalue', evalue,
    '-max_target_seqs', max_target_seqs,
    '-qcov_hsp_perc', qcov_hsp_perc,
    '-max_hsps', '1',
    '-num_threads', num_threads,
    '-outfmt', '6 qseqid sseqid stitle sscinames pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp',
    '-out', output_file
]

# Print the command being executed
print("\n========================================")
print("Running BLASTX command (STRICT FILTER):")
print("========================================")
print(' '.join(blastx_cmd))
print("========================================\n")

# Execute BLASTX
try:
    subprocess.run(blastx_cmd, check=True)
except subprocess.CalledProcessError as e:
    print(f"ERROR: BLASTX failed with exit code {e.returncode}", file=sys.stderr)
    sys.exit(1)

# Print completion message
print("\n========================================")
print(f"Done: BLASTX results written to {output_file}")
print("========================================")

# Show first 5 lines of output
output_path = Path(output_file)
if output_path.exists():
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
echo "BLASTX Analysis Completed (STRICT FILTER MODE)"
date
echo "=========================================="
