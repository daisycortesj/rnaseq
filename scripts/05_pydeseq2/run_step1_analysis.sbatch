#!/bin/bash
################################################################################
# PYDESEQ2 ANALYSIS - STEP 1: Get ALL genes (for BLAST)
################################################################################
#
# What this does:
#   Runs statistical analysis on ALL genes (no filtering)
#   Calculates log2FoldChange and padj for every gene
#   Output includes genes with bad p-values (you need them for BLAST!)
#
# How to run:
#   sbatch run_pydeseq2_step1.sbatch DC    (for D. carota)
#   sbatch run_pydeseq2_step1.sbatch DG    (for D. glaber)
#
# What you get:
#   pydeseq2_results_UNFILTERED.tsv - Every gene with statistics
#   (Use this file for BLAST annotation!)
#
################################################################################


# ══════════════════════════════════════════════════════════════════════════════
# SLURM SETTINGS - Change these if needed
# ══════════════════════════════════════════════════════════════════════════════

#SBATCH --job-name=pydeseq2_step1         # Job name
#SBATCH --account=tholl_lab_1             # Your lab account
#SBATCH --partition=normal_q              # Which queue to use
#SBATCH --cpus-per-task=8                 # Number of CPUs (8 for speed)
#SBATCH --mem=16G                         # Memory (16GB for PyDESeq2)
#SBATCH --time=2:00:00                    # Max time (2 hours)
#SBATCH --mail-type=ALL                   # Email when done
#SBATCH --mail-user=daisycortesj@vt.edu   # Your email

# Where to save output and error messages
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis
#SBATCH -o 06_analysis/pydeseq2_step1_%j.out
#SBATCH -e 06_analysis/pydeseq2_step1_%j.err


# ══════════════════════════════════════════════════════════════════════════════
# PATHS - Load shared config
# ══════════════════════════════════════════════════════════════════════════════

# ===== Load shared config =====
source "${SLURM_SUBMIT_DIR:-.}/scripts/config.sh"


# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS SETTINGS - Change these to modify your comparison
# ══════════════════════════════════════════════════════════════════════════════

# ("condition" for root vs leaf)
CONTRAST_FACTOR="condition"

# Which groups are you comparing?
#   CONTRAST_A = the condition you're interested in (numerator)
#   CONTRAST_B = the baseline (denominator)
#
#   log2FC = log2( A / B )
#     Positive log2FC --> gene is HIGHER in A
#     Negative log2FC --> gene is HIGHER in B
#
# Default: R vs L (positive log2FC = up in root)
# Override for DG (leaf focus):
#   CONTRAST_A=L CONTRAST_B=R sbatch scripts/05_pydeseq2/run_step1_analysis.sbatch DG
CONTRAST_A="${CONTRAST_A:-R}"
CONTRAST_B="${CONTRAST_B:-L}"

# Minimum total read counts to keep a gene (default: 20, matching previous student)
# Override: MIN_COUNTS=10 sbatch scripts/run_pydeseq2_step1_analysis.sbatch DC
MIN_COUNTS="${MIN_COUNTS:-20}" # minimum counts to keep a gene (default: 20)

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
    echo "  sbatch run_pydeseq2_step1.sbatch DC"
    echo "  sbatch run_pydeseq2_step1.sbatch DG"
    exit 1
fi

# Convert to uppercase
SPECIES=$(echo "${SPECIES}" | tr '[:lower:]' '[:upper:]')

# Set species name and directory from shared config
get_sample_info "${SPECIES}" || exit 1
COUNT_FOLDER="${SPECIES_DIR}"


# ══════════════════════════════════════════════════════════════════════════════
# SET UP FILE PATHS
# ══════════════════════════════════════════════════════════════════════════════

# Input files (from count matrix step)
COUNT_MATRIX="${BASE_DIR}/03_count_tables/${COUNT_FOLDER}/gene_count_matrix.tsv"
METADATA="${BASE_DIR}/03_count_tables/${COUNT_FOLDER}/sample_metadata.tsv"

# Output directory
OUTPUT_DIR="${BASE_DIR}/06_analysis/pydeseq2_${SPECIES}_step1_unfiltered"


# ══════════════════════════════════════════════════════════════════════════════
# SHOW WHAT WE'RE DOING
# ══════════════════════════════════════════════════════════════════════════════

echo "============================================================"
echo "PyDESeq2 STEP 1: Unfiltered Analysis"
echo "============================================================"
echo ""
echo "Species: ${SPECIES_NAME} (${SPECIES})"
echo "Start time: $(date)"
echo ""
echo "------------------------------------------------------------"
echo "INPUT FILES"
echo "------------------------------------------------------------"
echo "  Count matrix: ${COUNT_MATRIX}"
echo "  Metadata:     ${METADATA}"
echo ""
echo "------------------------------------------------------------"
echo "COMPARISON"
echo "------------------------------------------------------------"
LABEL_A=$( [ "${CONTRAST_A}" = "R" ] && echo "root" || echo "leaf" )
LABEL_B=$( [ "${CONTRAST_B}" = "L" ] && echo "leaf" || echo "root" )
echo "  Comparing: ${CONTRAST_A} (${LABEL_A}) vs ${CONTRAST_B} (${LABEL_B})"
echo "  Positive log2FC = higher in ${LABEL_A}"
echo "  Negative log2FC = higher in ${LABEL_B}"
echo "  Using column: ${CONTRAST_FACTOR}"
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

echo "------------------------------------------------------------"
echo "CHECKING FILES"
echo "------------------------------------------------------------"

# Check count matrix exists
if [ ! -f "${COUNT_MATRIX}" ]; then
    echo "ERROR: Count matrix not found!"
    echo "${COUNT_MATRIX}"
    echo ""
    echo "Did you run the count matrix step first?"
    echo "  sbatch scripts/build_count_matrix.sbatch ${SPECIES}"
    exit 1
fi
echo "✓ Count matrix found"

# Check metadata exists
if [ ! -f "${METADATA}" ]; then
    echo "ERROR: Metadata not found!"
    echo "${METADATA}"
    exit 1
fi
echo "✓ Metadata found"
echo ""


# ══════════════════════════════════════════════════════════════════════════════
# RUN PYDESEQ2 ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════

# Make output directory
mkdir -p "${OUTPUT_DIR}"

echo "============================================================"
echo "RUNNING PYDESEQ2 ANALYSIS"
echo "============================================================"
echo ""
echo "This will:"
echo "  1. Normalize counts (size factors)"
echo "  2. Estimate dispersion (gene variability)"
echo "  3. Calculate log2FoldChange (root vs leaf)"
echo "  4. Calculate p-values (statistical test)"
echo "  5. Adjust p-values (padj for multiple testing)"
echo ""
echo "Running..."
echo ""


# TO_DO: TEST OUT LEAF VS ROOT
# Run the Python script
python3 "${CODE_DIR}/scripts/05_pydeseq2/pydeseq2_run_analysis.py" \
    "${COUNT_MATRIX}" \
    "${METADATA}" \
    -o "${OUTPUT_DIR}" \
    --contrast-factor "${CONTRAST_FACTOR}" \
    --contrast-A "${CONTRAST_A}" \
    --contrast-B "${CONTRAST_B}" \
    --min-counts "${MIN_COUNTS}"

# Check if it worked
if [ $? -ne 0 ]; then
    echo ""
    echo "ERROR: PyDESeq2 analysis failed"
    exit 1
fi


# ══════════════════════════════════════════════════════════════════════════════
# DONE!
# ══════════════════════════════════════════════════════════════════════════════

echo ""
echo "============================================================"
echo "STEP 1 COMPLETE!"
echo "============================================================"
echo ""
echo "Species: ${SPECIES_NAME}"
echo "Output: ${OUTPUT_DIR}"
echo ""
echo "Files created:"
echo "  • pydeseq2_results_UNFILTERED.tsv"
echo "    (ALL genes with stats - use for BLAST combine step)"
echo ""
echo "  • all_gene_ids.txt"
echo "    (Gene candidate list - use for extract_cds_from_gtf.py)"
echo ""
echo "  • qc_total_counts.pdf, log2foldchange_histogram.pdf, padj_histogram.pdf"
echo "    (Quality control and distribution plots)"
echo ""
echo "------------------------------------------------------------"
echo "NEXT STEPS (two paths depending on what you need)"
echo "------------------------------------------------------------"
echo ""
echo "PATH A: Quick plots (MA, volcano, PCA — no heatmaps)"
echo "  1. sbatch scripts/05_pydeseq2/run_step2_filter.sbatch ${SPECIES}"
echo "  2. sbatch scripts/05_pydeseq2/run_step3_plots.sbatch ${SPECIES}"
echo ""
echo "PATH B: Full plots INCLUDING CYP/OMT heatmaps (needs BLAST first)"
echo "  1. sbatch scripts/06_blast/run_blast_input.sbatch ${SPECIES}"
echo "  2. sbatch scripts/06_blast/run_blastp_discovery.sbatch ${SPECIES} swissprot"
echo "  3. sbatch scripts/06_blast/run_combine_blast_deseq.sbatch ${SPECIES} swissprot discovery"
echo "  4. sbatch scripts/05_pydeseq2/run_step3_plots.sbatch ${SPECIES} combined_annotated"
echo ""
echo "  Why? CYP/OMT heatmaps need BLAST descriptions to identify which"
echo "  genes are cytochrome P450s or O-methyltransferases. Without BLAST,"
echo "  the script doesn't know which genes to put in the heatmaps."
echo ""
echo "End time: $(date)"
echo "============================================================"

exit 0
