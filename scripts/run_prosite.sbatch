#!/bin/bash
################################################################################
# Run PROSITE motif scanning for protein pattern identification
################################################################################
#
# Uses PROSITE patterns to identify functional motifs in protein sequences.
# PROSITE is complementary to HMMER (Pfam):
#   - PROSITE: Short, highly specific motifs/signatures (e.g., active sites)
#   - Pfam: Longer protein domains and families
#
# Use this AFTER:
#   1. Translating CDS to protein (run_translate_cds.sbatch)
#   2. Running HMMER (run_hmmer.sbatch) - for domain information
#   3. Optionally filtering to candidate genes
#
# How to run:
#   # Scan all proteins:
#   sbatch scripts/run_prosite.sbatch DC all_genes_protein.fasta
#
#   # Scan filtered proteins:
#   sbatch scripts/run_prosite.sbatch DC filtered_proteins.fasta
#
# Output:
#   - {filename}_prosite_motifs.txt (motif table)
#
################################################################################

#SBATCH --job-name=prosite_scan
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu

#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/06_analysis
#SBATCH -o prosite_%j.out
#SBATCH -e prosite_%j.err

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
    echo "Usage: sbatch run_prosite.sbatch [SPECIES] [protein.fasta]"
    exit 1
fi

# Set paths
BASE_DIR=/projects/tholl_lab_1/daisy_analysis
INPUT_DIR="${BASE_DIR}/06_analysis/blast_input_${SPECIES}"
OUTPUT_DIR="${BASE_DIR}/06_analysis/prosite_${SPECIES}"
PROSITE_DIR="${BASE_DIR}/07_NRdatabase/prosite"

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

# Validate PROSITE database exists
if [ ! -f "${PROSITE_DIR}/prosite.dat" ]; then
    echo "ERROR: PROSITE database not found: ${PROSITE_DIR}/prosite.dat"
    echo ""
    echo "Download PROSITE database:"
    echo "  sbatch scripts/download_protein_databases.sbatch"
    echo ""
    echo "Or manually:"
    echo "  mkdir -p ${PROSITE_DIR}"
    echo "  cd ${PROSITE_DIR}"
    echo "  wget ftp://ftp.expasy.org/databases/prosite/prosite.dat"
    echo "  wget ftp://ftp.expasy.org/databases/prosite/prosite.psa"
    exit 1
fi

# ============================================================================

# Create output directory
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

# Get base filename for output
BASENAME=$(basename "$QUERY_FASTA" .fasta)

# Print info
echo "=========================================="
echo "PROSITE Motif Scan Started"
echo "=========================================="
echo "Species:      $SPECIES"
echo "Query:        $QUERY_FASTA"
echo "PROSITE DB:   $PROSITE_DIR"
echo "Output dir:   $OUTPUT_DIR"
echo "CPUs:         ${SLURM_CPUS_PER_TASK:-1}"
echo "------------------------------------------"
date
echo "Hostname: $(hostname)"
echo "------------------------------------------"

# Check for EMBOSS PROSITE tools
if command -v patmatmotifs &> /dev/null && command -v prosextract &> /dev/null; then
    echo "Using EMBOSS tools:"
    echo "  prosextract: $(which prosextract)"
    echo "  patmatmotifs: $(which patmatmotifs)"
else
    echo "ERROR: EMBOSS PROSITE tools not found"
    echo ""
    echo "Required: prosextract, patmatmotifs"
    echo "Install with: conda install -c bioconda emboss"
    exit 1
fi

echo "=========================================="

# Step 1: Prepare PROSITE database for EMBOSS (if not already done)
echo ""
echo "Step 1: Preparing PROSITE database for EMBOSS..."

PROSITE_EMBOSS_DIR="${PROSITE_DIR}/emboss_data"
mkdir -p "$PROSITE_EMBOSS_DIR"

if [ ! -f "${PROSITE_EMBOSS_DIR}/PROSITE/prosite.lines" ]; then
    echo "Running prosextract to index PROSITE for EMBOSS..."
    
    # Set EMBOSS data directory
    export EMBOSS_DATA="$PROSITE_EMBOSS_DIR"
    
    prosextract \
        -prositedir "$PROSITE_DIR" \
        -infname prosite.dat
    
    if [ $? -ne 0 ]; then
        echo "ERROR: prosextract failed"
        exit 1
    fi
    echo "✓ PROSITE database indexed for EMBOSS"
else
    echo "✓ PROSITE database already indexed"
    export EMBOSS_DATA="$PROSITE_EMBOSS_DIR"
fi

# Step 2: Run patmatmotifs on each sequence
echo ""
echo "Step 2: Running PROSITE motif scan..."
echo "This may take 30 min - 4 hours depending on input size"
echo ""

# patmatmotifs works on individual sequences, so we process the FASTA
# and collect results into a single output file

echo -e "sequence_id\tprosite_id\tmotif_name\tstart\tend\tscore" > "${BASENAME}_prosite_motifs.txt"

patmatmotifs \
    -sequence "$QUERY_FASTA" \
    -outfile "${BASENAME}_prosite_raw.patmatmotifs" \
    -full \
    -auto

if [ $? -ne 0 ]; then
    echo "ERROR: patmatmotifs failed"
    exit 1
fi

# Parse patmatmotifs output into table format
if [ -f "${BASENAME}_prosite_raw.patmatmotifs" ]; then
    awk '
    /^Sequence:/ { seq = $2 }
    /^Motif =/ { motif = $3 }
    /^  Start =/ { start = $3 }
    /^  End =/ { end = $3 }
    /^  Score =/ { 
        score = $3
        if (seq != "" && motif != "") {
            print seq "\t" motif "\t" motif "\t" start "\t" end "\t" score
        }
    }
    ' "${BASENAME}_prosite_raw.patmatmotifs" >> "${BASENAME}_prosite_motifs.txt"
fi

echo ""
echo "=========================================="
echo "PROSITE Scan Completed"
echo "=========================================="
echo "Output files:"
echo "  Table: ${BASENAME}_prosite_motifs.txt"
echo "  Raw:   ${BASENAME}_prosite_raw.patmatmotifs"
echo ""

# Quick summary
if [ -f "${BASENAME}_prosite_motifs.txt" ]; then
    echo "Quick Summary:"
    echo "------------------------------------------"
    
    n_with_motifs=$(tail -n +2 "${BASENAME}_prosite_motifs.txt" | cut -f1 | sort -u | wc -l)
    echo "Proteins with motifs: $n_with_motifs"
    
    n_motifs=$(tail -n +2 "${BASENAME}_prosite_motifs.txt" | wc -l)
    echo "Total motif hits: $n_motifs"
    
    echo ""
    echo "Top 10 most common motifs:"
    tail -n +2 "${BASENAME}_prosite_motifs.txt" | cut -f2,3 | sort | uniq -c | sort -rn | head -10
    echo "------------------------------------------"
fi

echo ""
date
echo "=========================================="
