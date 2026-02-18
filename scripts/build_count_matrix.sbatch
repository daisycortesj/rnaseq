#!/bin/bash
################################################################################
# BUILD COUNT MATRIX - STEP 1
################################################################################
#
# What this does:
#   Takes individual count files and combines them into one big table
#   (genes as rows, samples as columns)
#
# How to run:
#   sbatch build_count_matrix.sbatch DC    (for D. carota)
#   sbatch build_count_matrix.sbatch DG    (for D. glaber)
#   sbatch build_count_matrix.sbatch MF    (for M. fragrans)
#
# What you get:
#   gene_count_matrix.tsv - Big table with all your samples
#   sample_metadata.tsv   - Info about each sample
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
#SBATCH -o 06_analysis/count_matrix_%j.out
#SBATCH -e 06_analysis/count_matrix_%j.err


# ══════════════════════════════════════════════════════════════════════════════
# PATHS - Change these if your folder structure is different
# ══════════════════════════════════════════════════════════════════════════════

BASE_DIR="/projects/tholl_lab_1/daisy_analysis"
CODE_DIR="${BASE_DIR}/05_rnaseq-code"


# ══════════════════════════════════════════════════════════════════════════════
# GET SPECIES FROM COMMAND LINE
# ══════════════════════════════════════════════════════════════════════════════

# Get the species code you typed (DC, DG, or MF)
SPECIES=$1

# Check if you forgot to type a species
if [ -z "${SPECIES}" ]; then
    echo "ERROR: You need to tell me which species!"
    echo ""
    echo "Run like this:"
    echo "  sbatch build_count_matrix.sbatch DC"
    echo "  sbatch build_count_matrix.sbatch DG"
    echo "  sbatch build_count_matrix.sbatch MF"
    exit 1
fi

# Convert to uppercase (so 'dc' and 'DC' both work)
SPECIES=$(echo "${SPECIES}" | tr '[:lower:]' '[:upper:]')

# Set the species name based on the code
if [ "${SPECIES}" = "DC" ]; then
    SPECIES_NAME="D. carota ssp. maximus"
elif [ "${SPECIES}" = "DG" ]; then
    SPECIES_NAME="Daucus glaber"
elif [ "${SPECIES}" = "MF" ]; then
    SPECIES_NAME="Myristica fragrans"
else
    echo "ERROR: Unknown species code: ${SPECIES}"
    echo "Valid options: DC, DG, MF"
    exit 1
fi


# ══════════════════════════════════════════════════════════════════════════════
# SET UP FILE PATHS
# ══════════════════════════════════════════════════════════════════════════════

# Where to find the organized count files
if [ "${SPECIES}" = "DC" ]; then
    INPUT_DIR="${BASE_DIR}/03_count_tables/00_1_DC"
elif [ "${SPECIES}" = "DG" ]; then
    INPUT_DIR="${BASE_DIR}/03_count_tables/00_2_DG"
elif [ "${SPECIES}" = "MF" ]; then
    INPUT_DIR="${BASE_DIR}/03_count_tables/00_3_MF"
fi

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
echo "Species code: ${SPECIES}"
echo "Start time: $(date)"
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

# Count how many files we have
NUM_FILES=$(ls ${INPUT_DIR}/*_ReadsPerGene.out.tab 2>/dev/null | wc -l)

if [ ${NUM_FILES} -eq 0 ]; then
    echo "ERROR: No count files found in ${INPUT_DIR}"
    exit 1
fi

echo "Found ${NUM_FILES} count files for ${SPECIES}"
echo ""


# ══════════════════════════════════════════════════════════════════════════════
# BUILD THE COUNT MATRIX
# ══════════════════════════════════════════════════════════════════════════════

# Make output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

echo "Building count matrix..."
echo ""

# Run the Python script to combine count files
python3 "${CODE_DIR}/build_count_matrix.py" \
    "${INPUT_DIR}" \
    -o "${OUTPUT_DIR}"

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
    echo "Next step:"
    echo "  sbatch scripts/run_pydeseq2_step1.sbatch ${SPECIES}"
else
    echo "ERROR: Count matrix creation failed"
    exit 1
fi

echo ""
echo "End time: $(date)"

exit 0
