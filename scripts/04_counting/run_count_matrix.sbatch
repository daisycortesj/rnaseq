#!/bin/bash
################################################################################
# BUILD COUNT MATRIX + SAMPLE METADATA
################################################################################
#
# What this does:
#   Takes count files and combines them into PyDESeq2-ready files:
#     gene_count_matrix.tsv  — genes as rows, samples as columns (integers)
#     sample_metadata.tsv    — which sample belongs to which group
#     count_summary.txt      — quick sanity-check stats
#
# Supports three counting paths:
#   PATH A (default for DC/DG) — STAR GeneCounts:
#     Uses *_ReadsPerGene.out.tab from STAR --quantMode GeneCounts
#     Input folder: 03_count_tables/{species}/
#
#   PATH B — featureCounts (reference-genome pipeline):
#     Uses featurecounts.txt from run_featurecounts.sbatch
#     Input folder: 03_count_tables/{species}/
#
#   PATH C — RSEM (Trinity de novo pipeline):
#     Uses RSEM.gene.counts.matrix from run_rsem.sbatch
#     Input folder:  01_processed/00_7_RSEM/
#     Output folder: 03_count_tables/00_5_MF_trinity/  (MF Trinity — see config.sh)
#
# How to run:
#   Reference genome (STAR):
#     sbatch scripts/04_counting/run_count_matrix.sbatch DC
#
#   Trinity + RSEM (after run_rsem.sbatch finishes):
#     sbatch scripts/04_counting/run_count_matrix.sbatch MF --type rsem
#
# What you get (all paths):
#   03_count_tables/{species}/gene_count_matrix.tsv
#   03_count_tables/{species}/sample_metadata.tsv
#   03_count_tables/{species}/count_summary.txt
#
################################################################################


# ══════════════════════════════════════════════════════════════════════════════
# SLURM SETTINGS
# ══════════════════════════════════════════════════════════════════════════════

#SBATCH --job-name=count_matrix
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daisycortesj@vt.edu

#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis
#SBATCH -o 03_count_tables/count_matrix_%j.out
#SBATCH -e 03_count_tables/count_matrix_%j.err


# ══════════════════════════════════════════════════════════════════════════════
# LOAD CONFIG
# ══════════════════════════════════════════════════════════════════════════════

source "${SLURM_SUBMIT_DIR:-.}/scripts/config.sh"


# ══════════════════════════════════════════════════════════════════════════════
# GET SPECIES AND COUNT TYPE FROM COMMAND LINE
# ══════════════════════════════════════════════════════════════════════════════

SPECIES=$1
COUNT_TYPE="auto"

# Parse optional --type argument (e.g. --type rsem)
shift || true
while [[ $# -gt 0 ]]; do
    case "$1" in
        --type)
            COUNT_TYPE="$2"
            shift 2
            ;;
        *)
            echo "WARNING: Unknown argument '$1', ignoring"
            shift
            ;;
    esac
done

if [ -z "${SPECIES}" ]; then
    echo "ERROR: You need to tell me which species!"
    echo ""
    echo "Run like this:"
    echo "  sbatch scripts/04_counting/run_count_matrix.sbatch DC"
    echo "  sbatch scripts/04_counting/run_count_matrix.sbatch MF --type rsem"
    exit 1
fi

SPECIES=$(echo "${SPECIES}" | tr '[:lower:]' '[:upper:]')
get_sample_info "${SPECIES}" || exit 1


# ══════════════════════════════════════════════════════════════════════════════
# SET UP INPUT AND OUTPUT PATHS
# ══════════════════════════════════════════════════════════════════════════════

# RSEM matrix from Trinity pipeline (written by run_rsem.sbatch)
RSEM_MATRIX="${PROCESSED_DIR}/00_7_RSEM/RSEM.gene.counts.matrix"

# Default input folder for STAR / featureCounts
STAR_INPUT_DIR="${COUNT_DIR}/${SPECIES_DIR}"

# Pick the input folder based on count type
if [ "${COUNT_TYPE}" = "rsem" ]; then
    INPUT_DIR="${PROCESSED_DIR}/00_7_RSEM"
elif [ "${COUNT_TYPE}" = "auto" ]; then
    # Auto: use RSEM if the Trinity matrix exists and there are no STAR files
    NUM_STAR=$(ls "${STAR_INPUT_DIR}"/*_ReadsPerGene.out.tab 2>/dev/null | wc -l)
    if [ -f "${RSEM_MATRIX}" ] && [ "${NUM_STAR}" -eq 0 ]; then
        COUNT_TYPE="rsem"
        INPUT_DIR="${PROCESSED_DIR}/00_7_RSEM"
        echo "Auto-detected: Trinity/RSEM path (RSEM.gene.counts.matrix found)"
    else
        INPUT_DIR="${STAR_INPUT_DIR}"
    fi
else
    INPUT_DIR="${STAR_INPUT_DIR}"
fi

# Where to save PyDESeq2-ready files
# Trinity/RSEM uses TRINITY_COUNT_DIR from config.sh (MF → 00_5_MF_trinity)
if [ "${COUNT_TYPE}" = "rsem" ] && [ -n "${TRINITY_COUNT_DIR:-}" ]; then
    OUTPUT_DIR="${COUNT_DIR}/${TRINITY_COUNT_DIR}"
else
    OUTPUT_DIR="${COUNT_DIR}/${SPECIES_DIR}"
fi


# ══════════════════════════════════════════════════════════════════════════════
# SHOW WHAT WE'RE DOING
# ══════════════════════════════════════════════════════════════════════════════

echo "============================================================"
echo "Building Count Matrix for ${SPECIES_NAME}"
echo "============================================================"
echo ""
echo "Species code:  ${SPECIES}"
echo "Count type:    ${COUNT_TYPE}"
echo "Start time:    $(date)"
echo ""
echo "Input folder:  ${INPUT_DIR}"
echo "Output folder: ${OUTPUT_DIR}"
echo ""


# ══════════════════════════════════════════════════════════════════════════════
# ACTIVATE CONDA ENVIRONMENT
# ══════════════════════════════════════════════════════════════════════════════

source ~/.bashrc
conda activate rnaseq

if [ -z "${CONDA_DEFAULT_ENV}" ]; then
    echo "ERROR: Conda environment didn't activate"
    exit 1
fi

echo "Using conda environment: ${CONDA_DEFAULT_ENV}"
echo ""


# ══════════════════════════════════════════════════════════════════════════════
# CHECK INPUT FILES EXIST
# ══════════════════════════════════════════════════════════════════════════════

if [ ! -d "${INPUT_DIR}" ]; then
    echo "ERROR: Input directory not found!"
    echo "${INPUT_DIR}"
    exit 1
fi

if [ "${COUNT_TYPE}" = "rsem" ]; then
    if [ ! -f "${INPUT_DIR}/RSEM.gene.counts.matrix" ]; then
        echo "ERROR: RSEM count matrix not found!"
        echo "${INPUT_DIR}/RSEM.gene.counts.matrix"
        echo ""
        echo "Run RSEM first:"
        echo "  sbatch scripts/04_rsem/run_rsem.sbatch ${SPECIES}"
        exit 1
    fi
    echo "Found RSEM.gene.counts.matrix"
elif [ "${COUNT_TYPE}" = "featurecounts" ]; then
    NUM_FILES=$(ls "${INPUT_DIR}"/*featurecounts* 2>/dev/null | grep -v '.summary' | wc -l)
    FILE_DESC="featureCounts"
elif [ "${COUNT_TYPE}" = "star" ]; then
    NUM_FILES=$(ls "${INPUT_DIR}"/*_ReadsPerGene.out.tab 2>/dev/null | wc -l)
    FILE_DESC="STAR ReadsPerGene"
else
    NUM_STAR=$(ls "${INPUT_DIR}"/*_ReadsPerGene.out.tab 2>/dev/null | wc -l)
    NUM_FC=$(ls "${INPUT_DIR}"/*featurecounts* 2>/dev/null | grep -v '.summary' | wc -l)
    NUM_FILES=$((NUM_STAR + NUM_FC))
    FILE_DESC="count (auto-detect)"
fi

if [ "${COUNT_TYPE}" != "rsem" ] && [ "${NUM_FILES}" -eq 0 ]; then
    echo "ERROR: No ${FILE_DESC} files found in ${INPUT_DIR}"
    echo ""
    echo "For Trinity pipeline, use:"
    echo "  sbatch scripts/04_counting/run_count_matrix.sbatch ${SPECIES} --type rsem"
    exit 1
fi

if [ "${COUNT_TYPE}" != "rsem" ]; then
    echo "Found ${NUM_FILES} ${FILE_DESC} file(s) for ${SPECIES}"
fi
echo ""


# ══════════════════════════════════════════════════════════════════════════════
# BUILD THE COUNT MATRIX
# ══════════════════════════════════════════════════════════════════════════════

mkdir -p "${OUTPUT_DIR}"

echo "Building count matrix (${COUNT_TYPE})..."
echo ""

# Base command
BUILD_CMD=(python3 "${CODE_DIR}/scripts/04_counting/build_count_matrix.py"
    "${INPUT_DIR}"
    -o "${OUTPUT_DIR}"
    --type "${COUNT_TYPE}")

# Pass SAMPLE_CONDITIONS from config.sh when set (needed for MF: MFF=F MFL=L)
if [[ -n "${SAMPLE_CONDITIONS:-}" ]]; then
    echo "Using condition map from config.sh: ${SAMPLE_CONDITIONS}"
    BUILD_CMD+=(--condition-map "${SAMPLE_CONDITIONS}")
fi

"${BUILD_CMD[@]}"

if [ $? -eq 0 ]; then
    echo ""
    echo "============================================================"
    echo "SUCCESS! Count matrix created"
    echo "============================================================"
    echo ""
    echo "Output files:"
    echo "  - gene_count_matrix.tsv"
    echo "  - sample_metadata.tsv"
    echo "  - count_summary.txt"
    echo ""
    echo "Location: ${OUTPUT_DIR}"
    echo ""
    echo "Samples included:"
    head -1 "${OUTPUT_DIR}/gene_count_matrix.tsv" | cut -f2- | tr '\t' '\n' | sed 's/^/  - /'
    echo ""
    echo "Next step — run PyDESeq2:"
    if [ "${SPECIES}" = "MF" ]; then
        echo "  CONTRAST_A=F CONTRAST_B=L sbatch scripts/05_pydeseq2/run_step1_analysis.sbatch MF"
    else
        echo "  sbatch scripts/05_pydeseq2/run_step1_analysis.sbatch ${SPECIES}"
    fi
else
    echo "ERROR: Count matrix creation failed"
    exit 1
fi

echo ""
echo "End time: $(date)"

exit 0
