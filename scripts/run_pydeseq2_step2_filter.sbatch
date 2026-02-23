#!/bin/bash
################################################################################
# PYDESEQ2 ANALYSIS - STEP 2: Filter Results (Get significant genes only)
################################################################################
#
# What this does:
#   Filters any DESeq2-like TSV by padj and log2FoldChange.
#   Accepts different input sources via optional 2nd argument.
#
# How to run:
#   # Filter raw unfiltered DESeq2 results (default):
#   sbatch scripts/run_pydeseq2_step2_filter.sbatch DC
#
#   # Filter combined BLAST+DESeq2 annotated file:
#   sbatch scripts/run_pydeseq2_step2_filter.sbatch DC combined_annotated
#
#   # Filter gene-family verified file (CYP+OMT):
#   sbatch scripts/run_pydeseq2_step2_filter.sbatch DC gene_families
#
#   # Filter a specific file by path:
#   sbatch scripts/run_pydeseq2_step2_filter.sbatch DC /path/to/any_results.tsv
#
# How to customize filters:
#   PADJ=0.01 LFC=3.0 sbatch scripts/run_pydeseq2_step2_filter.sbatch DC
#   PADJ=0.1 LFC=1.0 sbatch scripts/run_pydeseq2_step2_filter.sbatch DC
#   ROOT_UP_ONLY=true sbatch scripts/run_pydeseq2_step2_filter.sbatch DC
#   LEAF_UP_ONLY=true sbatch scripts/run_pydeseq2_step2_filter.sbatch DC
#
# What you get:
#   pydeseq2_{SP}_step2_filtered/{source}_FILTERED.tsv
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
PADJ_CUTOFF="${PADJ:-0.05}"

# Log2 fold change cutoff
# Standard: 2.0 (means 4-fold change)
# Strict:   3.0 (means 8-fold change)
# Lenient:  1.0 (means 2-fold change)
LFC_CUTOFF="${LFC:-2.0}"

# Direction filter (set to "true" if you only want one direction)
# ROOT_UP_ONLY: Keep only genes higher in root (positive log2FC)
# LEAF_UP_ONLY: Keep only genes higher in leaf (negative log2FC)
# Leave both as "false" to keep both directions
ROOT_UP_ONLY=false
LEAF_UP_ONLY=false


# ══════════════════════════════════════════════════════════════════════════════
# GET SPECIES AND INPUT SOURCE FROM COMMAND LINE
# ══════════════════════════════════════════════════════════════════════════════

SPECIES=$1
INPUT_SOURCE="${2:-unfiltered}"

if [ -z "${SPECIES}" ]; then
    echo "ERROR: You need to tell me which species!"
    echo ""
    echo "Usage:"
    echo "  sbatch scripts/run_pydeseq2_step2_filter.sbatch SPECIES [INPUT_SOURCE]"
    echo ""
    echo "  SPECIES:      DC, DG, or MF"
    echo "  INPUT_SOURCE: unfiltered (default), combined_annotated, gene_families,"
    echo "                or a full file path"
    echo ""
    echo "Examples:"
    echo "  sbatch scripts/run_pydeseq2_step2_filter.sbatch DC"
    echo "  sbatch scripts/run_pydeseq2_step2_filter.sbatch DC combined_annotated"
    echo "  sbatch scripts/run_pydeseq2_step2_filter.sbatch DC gene_families"
    exit 1
fi

SPECIES=$(echo "${SPECIES}" | tr '[:lower:]' '[:upper:]')

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
# SET UP FILE PATHS (based on INPUT_SOURCE)
# ══════════════════════════════════════════════════════════════════════════════

OUTPUT_DIR="${BASE_DIR}/06_analysis/pydeseq2_${SPECIES}_step2_filtered"

case "${INPUT_SOURCE}" in
    unfiltered)
        INPUT_FILE="${BASE_DIR}/06_analysis/pydeseq2_${SPECIES}_step1_unfiltered/pydeseq2_results_UNFILTERED.tsv"
        OUTPUT_FILE="${OUTPUT_DIR}/pydeseq2_results_FILTERED.tsv"
        SOURCE_LABEL="Step 1 unfiltered results"
        ;;
    combined_annotated)
        INPUT_FILE="${BASE_DIR}/06_analysis/combined_${SPECIES}/${SPECIES}_swissprot_discovery_annotated.tsv"
        OUTPUT_FILE="${OUTPUT_DIR}/combined_annotated_FILTERED.tsv"
        SOURCE_LABEL="Combined BLAST+DESeq2 annotated"
        ;;
    gene_families)
        INPUT_FILE="${BASE_DIR}/06_analysis/gene_families_${SPECIES}/${SPECIES}_CYP_OMT_combined.tsv"
        OUTPUT_FILE="${OUTPUT_DIR}/gene_families_FILTERED.tsv"
        SOURCE_LABEL="Gene family verified (CYP+OMT)"
        ;;
    *)
        INPUT_FILE="${INPUT_SOURCE}"
        BASENAME=$(basename "${INPUT_SOURCE}" .tsv)
        OUTPUT_FILE="${OUTPUT_DIR}/${BASENAME}_FILTERED.tsv"
        SOURCE_LABEL="Custom file"
        ;;
esac


# ══════════════════════════════════════════════════════════════════════════════
# SHOW WHAT WE'RE DOING
# ══════════════════════════════════════════════════════════════════════════════

echo "============================================================"
echo "PyDESeq2 STEP 2: Filter Results"
echo "============================================================"
echo ""
echo "Species:      ${SPECIES_NAME} (${SPECIES})"
echo "Input source: ${SOURCE_LABEL}"
echo "Input file:   ${INPUT_FILE}"
echo "Output file:  ${OUTPUT_FILE}"
echo "Start time:   $(date)"
echo ""
echo "------------------------------------------------------------"
echo "FILTER SETTINGS"
echo "------------------------------------------------------------"
echo "  padj cutoff:        ${PADJ_CUTOFF}"
echo "  |log2FC| cutoff:    ${LFC_CUTOFF}"
echo "  Root-up only:       ${ROOT_UP_ONLY}"
echo "  Leaf-up only:       ${LEAF_UP_ONLY}"
echo ""
echo "Keeping genes where:"
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

if [ ! -f "${INPUT_FILE}" ]; then
    echo "ERROR: Input file not found!"
    echo "  ${INPUT_FILE}"
    echo ""
    case "${INPUT_SOURCE}" in
        unfiltered)
            echo "Run Step 1 first:"
            echo "  sbatch scripts/run_pydeseq2_step1_analysis.sbatch ${SPECIES}"
            ;;
        combined_annotated)
            echo "Run combine+filter first:"
            echo "  sbatch scripts/run_combine_filter.sbatch ${SPECIES} swissprot discovery"
            ;;
        gene_families)
            echo "Run gene family extraction first:"
            echo "  sbatch scripts/run_gene_family_heatmap.sbatch ${SPECIES} swissprot discovery full"
            ;;
        *)
            echo "Check the file path and try again."
            ;;
    esac
    exit 1
fi
echo "Input file found"
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
echo "Reading input: ${SOURCE_LABEL}"
echo "Applying filters..."
echo ""

# Build the command
FILTER_CMD="python3 \"${CODE_DIR}/pydeseq2_filter_results.py\""
FILTER_CMD="${FILTER_CMD} \"${INPUT_FILE}\""
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
echo "Generate plots (MA, volcano, CYP/OMT heatmaps):"
echo "  sbatch scripts/run_pydeseq2_step3_plots.sbatch ${SPECIES} ${INPUT_SOURCE}"
echo ""
echo "End time: $(date)"
echo "============================================================"

exit 0
