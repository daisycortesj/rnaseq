#!/bin/bash
################################################################################
# Run HMMER hmmscan for protein domain identification
################################################################################
#
# Uses HMMER to scan protein sequences against Pfam database for domain
# identification and protein family classification.
#
# HMMER is more sensitive and specific than BLAST for domain detection.
# Use this AFTER:
#   1. Translating CDS to protein (run_translate_cds.sbatch)
#   2. Optionally filtering to candidate genes (filter_combined_results.py)
#
# How to run:
#   # Scan all proteins:
#   sbatch scripts/run_hmmer.sbatch DC all_genes_protein.fasta
#
#   # Scan filtered proteins:
#   sbatch scripts/run_hmmer.sbatch DC filtered_proteins.fasta
#
# Output:
#   - {filename}_pfam_domains.txt (domain table output)
#   - {filename}_pfam_domains.out (human-readable output)
#
################################################################################

#SBATCH --job-name=hmmer_scan
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu

#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/06_analysis
#SBATCH -o hmmer_%j.out
#SBATCH -e hmmer_%j.err

# Strict error handling
set -euo pipefail

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq

# ============================================================================
# USER-EDITABLE VARIABLES
# ============================================================================
SPECIES="${1:-DC}"
QUERY_FILE="${2:-all_genes_protein.fasta}"

# Validate species
if [[ ! "$SPECIES" =~ ^(DC|DG|MF)$ ]]; then
    echo "ERROR: Invalid species '$SPECIES'. Must be DC, DG, or MF"
    echo "Usage: sbatch run_hmmer.sbatch [SPECIES] [protein.fasta]"
    exit 1
fi

# Set paths
BASE_DIR=/projects/tholl_lab_1/daisy_analysis
INPUT_DIR="${BASE_DIR}/06_analysis/blast_input_${SPECIES}"
OUTPUT_DIR="${BASE_DIR}/06_analysis/hmmer_${SPECIES}"
PFAM_DB="${BASE_DIR}/07_NRdatabase/Pfam-A.hmm"

# Determine input file path
if [[ "$QUERY_FILE" == /* ]]; then
    # Absolute path provided
    QUERY_FASTA="$QUERY_FILE"
else
    # Relative path - look in blast_input directory
    QUERY_FASTA="${INPUT_DIR}/${QUERY_FILE}"
fi

# Validate input file exists
if [ ! -f "$QUERY_FASTA" ]; then
    echo "ERROR: Query file not found: $QUERY_FASTA"
    echo "Make sure to run: sbatch scripts/run_translate_cds.sbatch $SPECIES"
    exit 1
fi

# Validate Pfam database exists
if [ ! -f "$PFAM_DB" ]; then
    echo "ERROR: Pfam database not found: $PFAM_DB"
    echo ""
    echo "Download and prepare Pfam:"
    echo "  cd ${BASE_DIR}/07_NRdatabase"
    echo "  wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
    echo "  gunzip Pfam-A.hmm.gz"
    echo "  hmmpress Pfam-A.hmm"
    exit 1
fi

# E-value threshold (lower = more stringent)
EVALUE=1e-5  # Standard cutoff for domains

# ============================================================================

# Create output directory
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

# Get base filename for output
BASENAME=$(basename "$QUERY_FASTA" .fasta)

# Print info
echo "=========================================="
echo "HMMER Domain Scan Started"
echo "=========================================="
echo "Species:    $SPECIES"
echo "Query:      $QUERY_FASTA"
echo "Pfam DB:    $PFAM_DB"
echo "Output dir: $OUTPUT_DIR"
echo "E-value:    $EVALUE"
echo "CPUs:       ${SLURM_CPUS_PER_TASK:-1}"
echo "------------------------------------------"
date
echo "Hostname: $(hostname)"
echo "------------------------------------------"
which hmmscan
hmmscan -h | head -3
echo "=========================================="

# Run HMMER hmmscan
echo ""
echo "Running HMMER hmmscan..."
echo "This may take 30 min - 6 hours depending on input size"
echo ""

hmmscan \
    --cpu ${SLURM_CPUS_PER_TASK:-1} \
    --domtblout "${BASENAME}_pfam_domains.txt" \
    --noali \
    -E ${EVALUE} \
    "$PFAM_DB" \
    "$QUERY_FASTA" \
    > "${BASENAME}_pfam_domains.out"

# Check if successful
if [ $? -ne 0 ]; then
    echo "ERROR: HMMER failed"
    exit 1
fi

echo ""
echo "=========================================="
echo "HMMER Scan Completed"
echo "=========================================="
echo "Output files:"
echo "  Domain table: ${BASENAME}_pfam_domains.txt"
echo "  Full output:  ${BASENAME}_pfam_domains.out"
echo ""

# Quick summary
if [ -f "${BASENAME}_pfam_domains.txt" ]; then
    echo "Quick Summary:"
    echo "------------------------------------------"
    
    # Count proteins with domains
    n_with_domains=$(grep -v "^#" "${BASENAME}_pfam_domains.txt" | cut -f4 | sort -u | wc -l)
    echo "Proteins with domains: $n_with_domains"
    
    # Count total domain hits
    n_domains=$(grep -v "^#" "${BASENAME}_pfam_domains.txt" | wc -l)
    echo "Total domain hits: $n_domains"
    
    # Top 10 most common domains
    echo ""
    echo "Top 10 most common domains:"
    grep -v "^#" "${BASENAME}_pfam_domains.txt" | awk '{print $4,$5}' | sort | uniq -c | sort -rn | head -10
    echo "------------------------------------------"
fi

echo ""
date
echo "=========================================="
