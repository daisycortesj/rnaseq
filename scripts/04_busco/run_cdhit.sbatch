#!/bin/bash
# =============================================================================
# SCRIPT:  run_cdhit.sbatch
# WHAT IT DOES:
#   Runs CD-HIT-EST on a pooled Trinity assembly to collapse near-duplicate
#   transcripts (same gene, different isoforms or assembly artifacts).
#
# WHY RUN CD-HIT BEFORE RE-RUNNING BUSCO?
#   Trinity often produces many similar transcripts per gene. BUSCO counts
#   each match separately, which inflates the "Duplicated (D)" percentage.
#   CD-HIT keeps one representative sequence per cluster at 95% identity,
#   so BUSCO D% usually drops while Complete (C) stays similar.
#
# INPUT:
#   Pooled Trinity assembly (same paths as run_busco.sbatch):
#     MF   → 01_processed/00_4_MF_trinity/MF_trinity_pooled.Trinity.fasta
#     DC   → 01_processed/00_1_DC_trinity/DC_trinity_pooled.Trinity.fasta
#     DG   → 01_processed/00_2_DG_trinity/DG_trinity_pooled.Trinity.fasta
#     SK   → 01_processed/00_4_Sukman_trinity/SK_trinity_pooled.Trinity.fasta
#     DCDG → 01_processed/00_4_DC_DG_trinity/DCDG_trinity_pooled.Trinity.fasta
#
# OUTPUT:
#   01_processed/00_6_cdhit/<CODE>_trinity_cdhit95.fasta
#   01_processed/00_6_cdhit/<CODE>_trinity_cdhit95.fasta.clstr  (cluster info)
#
# HOW TO SUBMIT:
#   sbatch scripts/04_busco/run_cdhit.sbatch MF
#   sbatch scripts/04_busco/run_cdhit.sbatch DC
#
# AFTER CD-HIT FINISHES:
#   sbatch scripts/04_busco/run_busco.sbatch MF
#
# PATHS: All folder paths come from scripts/config.sh
# =============================================================================


# ---- SLURM cluster resource requests ----
# CD-HIT-EST is CPU-heavy; 16 cores and 64 GB RAM is a safe default.
#SBATCH --job-name=cdhit
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/01_processed/00_6_cdhit
#SBATCH -o cdhit_%j.out
#SBATCH -e cdhit_%j.err

set -euo pipefail


# =============================================================================
# LOAD CONFIG
# =============================================================================
CONFIG_FILE=""
for candidate in \
    "${SLURM_SUBMIT_DIR:-.}/scripts/config.sh" \
    "${SLURM_SUBMIT_DIR:-.}/config.sh" \
    "${SLURM_SUBMIT_DIR:-.}/../config.sh" \
    "${SLURM_SUBMIT_DIR:-.}/../../config.sh" \
    "/projects/tholl_lab_1/daisy_analysis/05_rnaseq-code/scripts/config.sh"
do
    if [[ -f "${candidate}" ]]; then
        CONFIG_FILE="${candidate}"
        break
    fi
done

if [[ -z "${CONFIG_FILE}" ]]; then
    echo "ERROR: Could not locate scripts/config.sh."
    echo "Submit from inside the rnaseq repo, e.g.:"
    echo "  cd /projects/tholl_lab_1/daisy_analysis/05_rnaseq-code"
    echo "  sbatch scripts/04_busco/run_cdhit.sbatch MF"
    exit 1
fi

echo "Loading config from: ${CONFIG_FILE}"
source "${CONFIG_FILE}"


# =============================================================================
# READ SPECIES CODE FROM COMMAND LINE
# =============================================================================
SPECIES_CODE="${1:-}"

if [[ -z "${SPECIES_CODE}" ]]; then
    echo "ERROR: No species code provided."
    echo ""
    echo "Usage:"
    echo "  sbatch scripts/04_busco/run_cdhit.sbatch <SPECIES_CODE>"
    echo ""
    echo "Examples:"
    echo "  sbatch scripts/04_busco/run_cdhit.sbatch MF"
    echo "  sbatch scripts/04_busco/run_cdhit.sbatch DC"
    exit 1
fi

get_sample_info "${SPECIES_CODE}"
SPECIES_CODE_UPPER=$(echo "${SPECIES_CODE}" | tr '[:lower:]' '[:upper:]')


# =============================================================================
# PICK TRINITY INPUT PATH (same mapping as run_busco.sbatch)
# =============================================================================
case "${SPECIES_CODE_UPPER}" in
    MF)
        TRINITY_DIR="${PROCESSED_DIR}/00_4_MF_trinity"
        INPUT="${TRINITY_DIR}/MF_trinity_pooled.Trinity.fasta"
        ;;
    DC)
        TRINITY_DIR="${PROCESSED_DIR}/00_1_DC_trinity"
        INPUT="${TRINITY_DIR}/DC_trinity_pooled.Trinity.fasta"
        ;;
    DG)
        TRINITY_DIR="${PROCESSED_DIR}/00_2_DG_trinity"
        INPUT="${TRINITY_DIR}/DG_trinity_pooled.Trinity.fasta"
        ;;
    SK)
        TRINITY_DIR="${PROCESSED_DIR}/00_4_Sukman_trinity"
        INPUT="${TRINITY_DIR}/SK_trinity_pooled.Trinity.fasta"
        ;;
    DCDG)
        TRINITY_DIR="${PROCESSED_DIR}/00_4_DC_DG_trinity"
        INPUT="${TRINITY_DIR}/DCDG_trinity_pooled.Trinity.fasta"
        ;;
    *)
        echo "ERROR: No Trinity input path configured for species '${SPECIES_CODE_UPPER}'."
        echo "Add a new case block to run_cdhit.sbatch (same as run_busco.sbatch)."
        exit 1
        ;;
esac


# =============================================================================
# OUTPUT PATHS
# =============================================================================
OUTPUT_DIR="${PROCESSED_DIR}/00_6_cdhit"
OUTPUT="${OUTPUT_DIR}/${SPECIES_CODE_UPPER}_trinity_cdhit95.fasta"

THREADS="${SLURM_CPUS_PER_TASK:-16}"
# CD-HIT -M flag is memory limit in MB; leave headroom below SLURM --mem.
CDHIT_MEM_MB=60000

mkdir -p "${OUTPUT_DIR}"


# =============================================================================
# JOB HEADER
# =============================================================================
echo "========================================================"
echo "CD-HIT-EST — collapse near-duplicate Trinity transcripts"
echo "Species code:  ${SPECIES_CODE_UPPER}"
echo "Species:       ${SPECIES_NAME}"
echo "Job ID:        ${SLURM_JOB_ID:-N/A}"
echo "Host:          $(hostname)"
echo "Started:       $(date)"
echo "========================================================"


# =============================================================================
# LOAD CD-HIT MODULE
# =============================================================================
# On ARC login node, check the exact module name first:
#   module spider cd-hit
echo ""
echo "Loading CD-HIT module..."
module load CD-HIT/4.8.1-GCC-12.3.0

if ! command -v cd-hit-est &> /dev/null; then
    echo "ERROR: cd-hit-est not found after module load."
    echo "Run: module spider cd-hit"
    echo "Then update the module load line in run_cdhit.sbatch."
    exit 1
fi

echo "CD-HIT version:"
cd-hit-est -h 2>&1 | head -3
echo "------------------------------------------------------------"


# =============================================================================
# CHECK INPUT
# =============================================================================
echo ""
echo "Checking input file..."
if [[ ! -f "${INPUT}" ]]; then
    echo "ERROR: Trinity assembly not found at:"
    echo "  ${INPUT}"
    echo ""
    echo "Run trinity_pooled.sbatch for ${SPECIES_CODE_UPPER} first."
    exit 1
fi

INPUT_COUNT=$(grep -c ">" "${INPUT}" || echo "0")
echo "Input transcripts: ${INPUT_COUNT}"
echo "  ${INPUT}"
echo "------------------------------------------------------------"


# =============================================================================
# RUN CD-HIT-EST
# =============================================================================
echo ""
echo "Running CD-HIT-EST at 95% identity..."
echo "Output: ${OUTPUT}"
echo "------------------------------------------------------------"

cd-hit-est \
    -i "${INPUT}" \
    -o "${OUTPUT}" \
    -c 0.95 \
    -n 8 \
    -T "${THREADS}" \
    -M "${CDHIT_MEM_MB}" \
    -d 0 \
    -aS 0.90

echo "------------------------------------------------------------"


# =============================================================================
# REPORT RESULTS
# =============================================================================
OUTPUT_COUNT=$(grep -c ">" "${OUTPUT}" || echo "0")

echo ""
echo "========================================================"
echo "CD-HIT complete"
echo "========================================================"
echo "Species:         ${SPECIES_CODE_UPPER} (${SPECIES_NAME})"
echo "Input:           ${INPUT_COUNT} transcripts"
echo "Output:          ${OUTPUT_COUNT} transcripts"
echo "Reduction:       ${INPUT_COUNT} → ${OUTPUT_COUNT}"
echo "FASTA:           ${OUTPUT}"
echo "Cluster file:    ${OUTPUT}.clstr"
echo ""
echo "Next step — rerun BUSCO on the CD-HIT assembly:"
echo "  sbatch scripts/04_busco/run_busco.sbatch ${SPECIES_CODE_UPPER} cdhit"
echo "Finished: $(date)"
echo "========================================================"

exit 0
