#!/bin/bash
################################################################################
# PYDESEQ2 ANALYSIS - STEP 3: Generate Plots
################################################################################
#
# What this does:
#   Takes filtered results and creates publication-quality plots
#   - MA plot (shows all genes, highlights significant ones)
#   - Volcano plot (shows fold change vs significance)
#   - Heatmap (shows CYP gene expression patterns)
#
# How to run:
#   sbatch run_pydeseq2_step3.sbatch DC    (for D. carota)
#   sbatch run_pydeseq2_step3.sbatch DG    (for D. glaber)
#
# What you get:
#   ma_plot.pdf       - MA plot (baseMean vs log2FC)
#   volcano_plot.pdf  - Volcano plot (log2FC vs -log10(padj))
#   cyp_heatmap.pdf   - CYP gene heatmap (if you have CYP genes)
#
################################################################################


# ══════════════════════════════════════════════════════════════════════════════
# SLURM SETTINGS - Change these if needed
# ══════════════════════════════════════════════════════════════════════════════

#SBATCH --job-name=pydeseq2_step3         # Job name
#SBATCH --account=tholl_lab_1             # Your lab account
#SBATCH --partition=normal_q              # Which queue to use
#SBATCH --cpus-per-task=4                 # Number of CPUs (4 for plotting)
#SBATCH --mem=16G                         # Memory (16GB for heatmaps)
#SBATCH --time=1:00:00                    # Max time (1 hour)
#SBATCH --mail-type=ALL                   # Email when done
#SBATCH --mail-user=daisycortesj@vt.edu   # Your email

# Where to save output and error messages
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis
#SBATCH -o 06_analysis/pydeseq2_step3_%j.out
#SBATCH -e 06_analysis/pydeseq2_step3_%j.err


# ══════════════════════════════════════════════════════════════════════════════
# PATHS - Change these if your folder structure is different
# ══════════════════════════════════════════════════════════════════════════════

BASE_DIR="/projects/tholl_lab_1/daisy_analysis"
CODE_DIR="${BASE_DIR}/05_rnaseq-code"


# ══════════════════════════════════════════════════════════════════════════════
# PLOT SETTINGS - CHANGE THESE TO CUSTOMIZE YOUR PLOTS
# ══════════════════════════════════════════════════════════════════════════════

# Heatmap scaling method
# "center" = subtract mean (shows high/low relative to gene average)
# "zscore" = show standard deviations from mean (more extreme colors)
SCALE_METHOD="center"

# Cluster rows in heatmap?
# "true"  = group similar genes together (default)
# "false" = keep genes in original order
CLUSTER_ROWS=true

# Comparison settings (should match Step 1)
CONTRAST_FACTOR="condition"
CONTRAST_A="R"
CONTRAST_B="L"


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
    echo "  sbatch run_pydeseq2_step3.sbatch DC"
    echo "  sbatch run_pydeseq2_step3.sbatch DG"
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

# Input files
FILTERED_FILE="${BASE_DIR}/06_analysis/pydeseq2_${SPECIES}_step2_filtered/pydeseq2_results_FILTERED.tsv"
COUNT_MATRIX="${BASE_DIR}/06_analysis/count_matrices_${SPECIES}/gene_count_matrix.tsv"
METADATA="${BASE_DIR}/06_analysis/count_matrices_${SPECIES}/sample_metadata.tsv"

# Output directory
OUTPUT_DIR="${BASE_DIR}/06_analysis/pydeseq2_${SPECIES}_step3_plots"


# ══════════════════════════════════════════════════════════════════════════════
# SHOW WHAT WE'RE DOING
# ══════════════════════════════════════════════════════════════════════════════

echo "============================================================"
echo "PyDESeq2 STEP 3: Generate Plots"
echo "============================================================"
echo ""
echo "Species: ${SPECIES_NAME} (${SPECIES})"
echo "Start time: $(date)"
echo ""
echo "------------------------------------------------------------"
echo "PLOT SETTINGS"
echo "------------------------------------------------------------"
echo "  Scaling method:  ${SCALE_METHOD}"
echo "  Cluster rows:    ${CLUSTER_ROWS}"
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

# Check filtered results exist
if [ ! -f "${FILTERED_FILE}" ]; then
    echo "ERROR: Filtered results not found!"
    echo "${FILTERED_FILE}"
    echo ""
    echo "Did you run Step 2 first?"
    echo "  sbatch scripts/run_pydeseq2_step2.sbatch ${SPECIES}"
    exit 1
fi
echo "✓ Filtered results found"

# Check count matrix exists
if [ ! -f "${COUNT_MATRIX}" ]; then
    echo "ERROR: Count matrix not found!"
    echo "${COUNT_MATRIX}"
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
# GENERATE PLOTS
# ══════════════════════════════════════════════════════════════════════════════

# Make output directory
mkdir -p "${OUTPUT_DIR}"

echo "============================================================"
echo "GENERATING PLOTS"
echo "============================================================"
echo ""
echo "This will create:"
echo "  • MA plot (baseMean vs log2FoldChange)"
echo "  • Volcano plot (log2FC vs significance)"
echo "  • CYP heatmap (if CYP genes found)"
echo ""
echo "Running..."
echo ""

# Build the command
PLOT_CMD="python3 \"${CODE_DIR}/pydeseq2_generate_plots.py\""
PLOT_CMD="${PLOT_CMD} \"${FILTERED_FILE}\""
PLOT_CMD="${PLOT_CMD} -o \"${OUTPUT_DIR}\""
PLOT_CMD="${PLOT_CMD} --count-matrix \"${COUNT_MATRIX}\""
PLOT_CMD="${PLOT_CMD} --metadata \"${METADATA}\""
PLOT_CMD="${PLOT_CMD} --contrast-factor \"${CONTRAST_FACTOR}\""
PLOT_CMD="${PLOT_CMD} --contrast-A \"${CONTRAST_A}\""
PLOT_CMD="${PLOT_CMD} --contrast-B \"${CONTRAST_B}\""
PLOT_CMD="${PLOT_CMD} --scale ${SCALE_METHOD}"

# Add no-cluster flag if needed
if [ "${CLUSTER_ROWS}" = "false" ]; then
    PLOT_CMD="${PLOT_CMD} --no-row-cluster"
fi

# Run the plotting
eval ${PLOT_CMD}

# Check if it worked
if [ $? -ne 0 ]; then
    echo "ERROR: Plot generation failed"
    exit 1
fi


# ══════════════════════════════════════════════════════════════════════════════
# DONE!
# ══════════════════════════════════════════════════════════════════════════════

echo ""
echo "============================================================"
echo "STEP 3 COMPLETE! ALL ANALYSIS DONE!"
echo "============================================================"
echo ""
echo "Species: ${SPECIES_NAME}"
echo "Output: ${OUTPUT_DIR}"
echo ""
echo "Files created:"
echo "  • ma_plot.pdf"
echo "  • volcano_plot.pdf"

# Check if CYP heatmap was created
if [ -f "${OUTPUT_DIR}/cyp_heatmap.pdf" ]; then
    echo "  • cyp_heatmap.pdf"
    echo "  • cyp_heatmap_matrix.tsv"
fi

echo ""
echo "------------------------------------------------------------"
echo "ALL PYDESEQ2 STEPS COMPLETE FOR ${SPECIES}!"
echo "------------------------------------------------------------"
echo ""
echo "You now have:"
echo "  1. Unfiltered results (for BLAST)"
echo "  2. Filtered results (significant genes)"
echo "  3. Publication plots"
echo ""
echo "Next: Review your plots and run BLAST annotation"
echo ""
echo "End time: $(date)"
echo "============================================================"

exit 0
