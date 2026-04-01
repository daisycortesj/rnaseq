#!/bin/bash
################################################################################
# BUILD COUNT MATRIX
################################################################################
#
# What this does:
#   Takes count files and combines them into one big table
#   (genes as rows, samples as columns)
#
# Supports two counting paths:
#   PATH A (default) — STAR GeneCounts:
#     Uses *_ReadsPerGene.out.tab from STAR --quantMode GeneCounts
#
#   PATH B — featureCounts (previous student's workflow):
#     Uses featurecounts.txt from featurecounts.sbatch
#
# How to run:
#   sbatch scripts/build_count_matrix.sbatch DC                       # auto-detect
#   sbatch scripts/build_count_matrix.sbatch DC --type star           # force STAR path
#   sbatch scripts/build_count_matrix.sbatch DC --type featurecounts  # force featureCounts path
#
# What you get:
#   gene_count_matrix.tsv - Big table with all your samples
#   sample_metadata.tsv   - Info about each sample
#   count_summary.txt     - Quick stats
#
################################################################################


# ══════════════════════════════════════════════════════════════════════════════
# SLURM SETTINGS - Change these if needed
# ══════════════════════════════════════════════════════════════════════════════

#SBATCH --job-name=count_matrix           # Job name
#SBATCH --account=tholl_lab_1             # Your lab account
#SBATCH --partition=normal_q              # Which queue to use
#SBATCH --cpus-per-task=2                 # Number of CPUs (2 is enough)
#SBATCH --mem=8G                          # Memory (8GB is enough)
#SBATCH --time=30:00                      # Max time (30 minutes)
#SBATCH --mail-type=ALL                   # Email when done
#SBATCH --mail-user=daisycortesj@vt.edu   # Your email

# Where to save output and error messages
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis
#SBATCH -o 03_count_tables/count_matrix_%j.out
#SBATCH -e 03_count_tables/count_matrix_%j.err


# ══════════════════════════════════════════════════════════════════════════════
# PATHS - Change these if your folder structure is different
# ══════════════════════════════════════════════════════════════════════════════

# ===== Load shared config =====
source "${SLURM_SUBMIT_DIR:-.}/scripts/config.sh"


# ══════════════════════════════════════════════════════════════════════════════
# GET SPECIES FROM COMMAND LINE
# ══════════════════════════════════════════════════════════════════════════════

# Get the species code and optional --type flag
SPECIES=$1
COUNT_TYPE="auto"

# Parse optional --type argument
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

# Check if you forgot to type a species
if [ -z "${SPECIES}" ]; then
    echo "ERROR: You need to tell me which species!"
    echo ""
    echo "Run like this:"
    echo "  sbatch scripts/build_count_matrix.sbatch DC"
    echo "  sbatch scripts/build_count_matrix.sbatch DC --type star"
    echo "  sbatch scripts/build_count_matrix.sbatch DC --type featurecounts"
    exit 1
fi

# Convert to uppercase (so 'dc' and 'DC' both work)
SPECIES=$(echo "${SPECIES}" | tr '[:lower:]' '[:upper:]')

# Map species code → directory name and display name (uses shared config)
get_sample_info "${SPECIES}" || exit 1


# ══════════════════════════════════════════════════════════════════════════════
# SET UP FILE PATHS
# ══════════════════════════════════════════════════════════════════════════════

# Where to find the organized count files
INPUT_DIR="${BASE_DIR}/03_count_tables/${SPECIES_DIR}"

# Where to save the count matrix
# same folder as the organized count files
OUTPUT_DIR="${INPUT_DIR}"


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

# Load conda
source ~/.bashrc

# Activate your RNA-seq environment
conda activate rnaseq

# Check it worked
if [ -z "${CONDA_DEFAULT_ENV}" ]; then
    echo "ERROR: Conda environment didn't activate"
    exit 1
fi

echo "Using conda environment: ${CONDA_DEFAULT_ENV}"
echo ""


# ══════════════════════════════════════════════════════════════════════════════
# CHECK FILES EXIST
# ══════════════════════════════════════════════════════════════════════════════

# Check if input directory exists
if [ ! -d "${INPUT_DIR}" ]; then
    echo "ERROR: Input directory not found!"
    echo "${INPUT_DIR}"
    echo ""
    echo "Did you run the organize script first?"
    echo "  bash scripts/organize_counts_by_species.sh"
    exit 1
fi

# Check that count files exist for the chosen type
if [ "${COUNT_TYPE}" = "featurecounts" ]; then
    NUM_FILES=$(ls ${INPUT_DIR}/*featurecounts* 2>/dev/null | grep -v '.summary' | wc -l)
    FILE_DESC="featureCounts"
elif [ "${COUNT_TYPE}" = "star" ]; then
    NUM_FILES=$(ls ${INPUT_DIR}/*_ReadsPerGene.out.tab 2>/dev/null | wc -l)
    FILE_DESC="STAR ReadsPerGene"
else
    # auto — check both
    NUM_STAR=$(ls ${INPUT_DIR}/*_ReadsPerGene.out.tab 2>/dev/null | wc -l)
    NUM_FC=$(ls ${INPUT_DIR}/*featurecounts* 2>/dev/null | grep -v '.summary' | wc -l)
    NUM_FILES=$((NUM_STAR + NUM_FC))
    FILE_DESC="count (auto-detect)"
fi

if [ ${NUM_FILES} -eq 0 ]; then
    echo "ERROR: No ${FILE_DESC} files found in ${INPUT_DIR}"
    exit 1
fi

echo "Found ${NUM_FILES} ${FILE_DESC} file(s) for ${SPECIES}"
echo ""


# ══════════════════════════════════════════════════════════════════════════════
# BUILD THE COUNT MATRIX
# ══════════════════════════════════════════════════════════════════════════════

mkdir -p "${OUTPUT_DIR}"

echo "Building count matrix (${COUNT_TYPE})..."
echo ""

python3 "${CODE_DIR}/scripts/04_counting/build_count_matrix.py" \
    "${INPUT_DIR}" \
    -o "${OUTPUT_DIR}" \
    --type "${COUNT_TYPE}"

# Check if it worked
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
    echo "  sbatch scripts/run_pydeseq2_step1_analysis.sbatch ${SPECIES}"
else
    echo "ERROR: Count matrix creation failed"
    exit 1
fi

echo ""
echo "End time: $(date)"

exit 0
