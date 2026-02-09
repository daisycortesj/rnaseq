#!/bin/bash
################################################################################
# PYDESEQ2 ANALYSIS - STEP 2: Filter Results (Get significant genes only)
################################################################################
#
# What this does:
#   Takes the UNFILTERED results and keeps only significant genes
#   Filters by padj (adjusted p-value) and log2FoldChange
#
# How to run:
#   sbatch run_pydeseq2_step2.sbatch DC    (standard filtering)
#
# How to customize filters:
#   PADJ=0.01 LFC=3.0 sbatch run_pydeseq2_step2.sbatch DC  (very strict)
#   PADJ=0.1 LFC=1.0 sbatch run_pydeseq2_step2.sbatch DC   (very lenient)
#   ROOT_UP_ONLY=true sbatch run_pydeseq2_step2.sbatch DC  (root genes only)
#   LEAF_UP_ONLY=true sbatch run_pydeseq2_step2.sbatch DC  (leaf genes only)
#
# What you get:
#   pydeseq2_results_FILTERED.tsv - Only significant genes
#
################################################################################


# ══════════════════════════════════════════════════════════════════════════════
# SLURM SETTINGS - Change these if needed
# ══════════════════════════════════════════════════════════════════════════════

#SBATCH --job-name=pydeseq2_step2         # Job name
#SBATCH --account=tholl_lab_1             # Your lab account
#SBATCH --partition=normal_q              # Which queue to use
#SBATCH --cpus-per-task=2                 # Number of CPUs (2 is enough)
#SBATCH --mem=8G                          # Memory (8GB is enough)
#SBATCH --time=30:00                      # Max time (30 minutes)
#SBATCH --mail-type=ALL                   # Email when done
#SBATCH --mail-user=daisycortesj@vt.edu   # Your email

# Where to save output and error messages
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis
#SBATCH -o 06_analysis/pydeseq2_step2_%j.out
#SBATCH -e 06_analysis/pydeseq2_step2_%j.err


# ══════════════════════════════════════════════════════════════════════════════
# PATHS - Change these if your folder structure is different
# ══════════════════════════════════════════════════════════════════════════════

BASE_DIR="/projects/tholl_lab_1/daisy_analysis"
CODE_DIR="${BASE_DIR}/05_rnaseq-code"


# ══════════════════════════════════════════════════════════════════════════════
# FILTER SETTINGS - CHANGE THESE TO ADJUST YOUR FILTERS!
# ══════════════════════════════════════════════════════════════════════════════

# Adjusted p-value cutoff (padj)
# Standard: 0.05 (5% false discovery rate)
# Strict:   0.01 (1% false discovery rate)
# Lenient:  0.1  (10% false discovery rate)
PADJ_CUTOFF=0.05

# Log2 fold change cutoff
# Standard: 2.0 (means 4-fold change)
# Strict:   3.0 (means 8-fold change)
# Lenient:  1.0 (means 2-fold change)
LFC_CUTOFF=2.0

# Direction filter (set to "true" if you only want one direction)
# ROOT_UP_ONLY: Keep only genes higher in root (positive log2FC)
# LEAF_UP_ONLY: Keep only genes higher in leaf (negative log2FC)
# Leave both as "false" to keep both directions
ROOT_UP_ONLY=false
LEAF_UP_ONLY=false


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
    echo "  sbatch run_pydeseq2_step2.sbatch DC"
    echo "  sbatch run_pydeseq2_step2.sbatch DG"
    exit 1
fi

# Convert to uppercase
SPECIES=$(echo "${SPECIES}" | tr '[:lower:]' '[:upper:]')

# Set the species name
if [ "${SPECIES}" = "DC" ]; then
    SPECIES_NAME="D. carota ssp. maximus"
elif [ "${SPECIES}" = "DG" ]; then
    SPECIES_NAME="Daucus glaber"
elif [ "${SPECIES}" = "MF" ]; then
    SPECIES_NAME="Myristica fragrans"
else
    echo "ERROR: Unknown species: ${SPECIES}"
    exit 1
fi


# ══════════════════════════════════════════════════════════════════════════════
# SET UP FILE PATHS
# ══════════════════════════════════════════════════════════════════════════════

# Input file (from Step 1 - unfiltered results)
UNFILTERED_FILE="${BASE_DIR}/06_analysis/pydeseq2_${SPECIES}_step1_unfiltered/pydeseq2_results_UNFILTERED.tsv"

# Output directory
OUTPUT_DIR="${BASE_DIR}/06_analysis/pydeseq2_${SPECIES}_step2_filtered"
OUTPUT_FILE="${OUTPUT_DIR}/pydeseq2_results_FILTERED.tsv"


# ══════════════════════════════════════════════════════════════════════════════
# SHOW WHAT WE'RE DOING
# ══════════════════════════════════════════════════════════════════════════════

echo "============================================================"
echo "PyDESeq2 STEP 2: Filter Results"
echo "============================================================"
echo ""
echo "Species: ${SPECIES_NAME} (${SPECIES})"
echo "Start time: $(date)"
echo ""
echo "------------------------------------------------------------"
echo "FILTER SETTINGS"
echo "------------------------------------------------------------"
echo "  padj cutoff:        ${PADJ_CUTOFF}"
echo "  |log2FC| cutoff:    ${LFC_CUTOFF}"
echo "  Root-up only:       ${ROOT_UP_ONLY}"
echo "  Leaf-up only:       ${LEAF_UP_ONLY}"
echo ""
echo "This means we'll keep genes where:"
echo "  - padj < ${PADJ_CUTOFF}  (statistically significant)"
echo "  - |log2FC| > ${LFC_CUTOFF}  (big enough change)"


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

echo ""
echo "Using conda environment: ${CONDA_DEFAULT_ENV}"
echo ""


# ══════════════════════════════════════════════════════════════════════════════
# CHECK FILES EXIST
# ══════════════════════════════════════════════════════════════════════════════

echo "------------------------------------------------------------"
echo "CHECKING FILES"
echo "------------------------------------------------------------"

# Check unfiltered results exist
if [ ! -f "${UNFILTERED_FILE}" ]; then
    echo "ERROR: Unfiltered results not found!"
    echo "${UNFILTERED_FILE}"
    echo ""
    echo "Did you run Step 1 first?"
    echo "  sbatch scripts/run_pydeseq2_step1.sbatch ${SPECIES}"
    exit 1
fi
echo "✓ Unfiltered results found"
echo ""


# ══════════════════════════════════════════════════════════════════════════════
# RUN FILTERING
# ══════════════════════════════════════════════════════════════════════════════

# Make output directory
mkdir -p "${OUTPUT_DIR}"

echo "============================================================"
echo "FILTERING RESULTS"
echo "============================================================"
echo ""
echo "Reading unfiltered results..."
echo "Applying filters..."
echo ""

# Build the command
FILTER_CMD="python3 \"${CODE_DIR}/pydeseq2_filter_results.py\""
FILTER_CMD="${FILTER_CMD} \"${UNFILTERED_FILE}\""
FILTER_CMD="${FILTER_CMD} -o \"${OUTPUT_FILE}\""
FILTER_CMD="${FILTER_CMD} --padj ${PADJ_CUTOFF}"
FILTER_CMD="${FILTER_CMD} --lfc ${LFC_CUTOFF}"

# Add direction filters if set
if [ "${ROOT_UP_ONLY}" = "true" ]; then
    FILTER_CMD="${FILTER_CMD} --root-up-only"
fi

if [ "${LEAF_UP_ONLY}" = "true" ]; then
    FILTER_CMD="${FILTER_CMD} --leaf-up-only"
fi

# Run the filter
eval ${FILTER_CMD}

# Check if it worked
if [ $? -ne 0 ]; then
    echo "ERROR: Filtering failed"
    exit 1
fi


# ══════════════════════════════════════════════════════════════════════════════
# DONE!
# ══════════════════════════════════════════════════════════════════════════════

echo ""
echo "============================================================"
echo "STEP 2 COMPLETE!"
echo "============================================================"
echo ""
echo "Species: ${SPECIES_NAME}"
echo "Output: ${OUTPUT_FILE}"
echo ""

# Count how many genes passed the filter
if [ -f "${OUTPUT_FILE}" ]; then
    NUM_GENES=$(tail -n +2 "${OUTPUT_FILE}" | wc -l)
    echo "Significant genes found: ${NUM_GENES}"
fi

echo ""
echo "------------------------------------------------------------"
echo "NEXT STEPS"
echo "------------------------------------------------------------"
echo ""
echo "Generate plots:"
echo "  sbatch scripts/run_pydeseq2_step3.sbatch ${SPECIES}"
echo ""
echo "End time: $(date)"
echo "============================================================"

exit 0
