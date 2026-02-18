#!/bin/bash
################################################################################
# PYDESEQ2 ANALYSIS - STEP 3: Generate Plots
################################################################################
#
# What this does:
#   Creates publication-quality plots from DESeq2 results:
#   - MA plot (baseMean vs log2FoldChange)
#   - Volcano plot (log2FC vs -log10 padj)
#   - CYP heatmap (auto-detected from BLAST annotations or gene_family column)
#   - OMT heatmap (auto-detected from BLAST annotations or gene_family column)
#
#   CYP/OMT genes are identified automatically if the input has a
#   blast_description or gene_family column. HMMER domain labels are added
#   to heatmap rows if the HMMER file exists.
#
# How to run:
#   # From Step 2 filtered results (default):
#   sbatch scripts/run_pydeseq2_step3_plots.sbatch DC
#
#   # From combined BLAST+DESeq2 annotated file (skip Step 2):
#   sbatch scripts/run_pydeseq2_step3_plots.sbatch DC combined_annotated
#
#   # From gene-family verified file:
#   sbatch scripts/run_pydeseq2_step3_plots.sbatch DC gene_families
#
#   # From a specific file:
#   sbatch scripts/run_pydeseq2_step3_plots.sbatch DC /path/to/results.tsv
#
# What you get:
#   ma_plot.pdf          - MA plot
#   volcano_plot.pdf     - Volcano plot
#   cyp_heatmap.pdf      - CYP gene heatmap (if CYP genes detected)
#   omt_heatmap.pdf      - OMT gene heatmap (if OMT genes detected)
#   cyp_gene_list.tsv    - List of detected CYP genes
#   omt_gene_list.tsv    - List of detected OMT genes
#
################################################################################


# ══════════════════════════════════════════════════════════════════════════════
# SLURM SETTINGS
# ══════════════════════════════════════════════════════════════════════════════

#SBATCH --job-name=pydeseq2_step3
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu

#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis
#SBATCH -o 06_analysis/pydeseq2_step3_%j.out
#SBATCH -e 06_analysis/pydeseq2_step3_%j.err


# ══════════════════════════════════════════════════════════════════════════════
# PATHS
# ══════════════════════════════════════════════════════════════════════════════

BASE_DIR="/projects/tholl_lab_1/daisy_analysis"
CODE_DIR="${BASE_DIR}/05_rnaseq-code"


# ══════════════════════════════════════════════════════════════════════════════
# PLOT SETTINGS
# ══════════════════════════════════════════════════════════════════════════════

SCALE_METHOD="center"
CLUSTER_ROWS=true
CONTRAST_FACTOR="condition"
CONTRAST_A="R"
CONTRAST_B="L"


# ══════════════════════════════════════════════════════════════════════════════
# ARGUMENTS
# ══════════════════════════════════════════════════════════════════════════════

SPECIES="${1:-}"
INPUT_SOURCE="${2:-filtered}"

if [ -z "${SPECIES}" ]; then
    echo "ERROR: Missing species argument"
    echo ""
    echo "Usage: sbatch scripts/run_pydeseq2_step3_plots.sbatch SPECIES [INPUT_SOURCE]"
    echo ""
    echo "  SPECIES:      DC, DG, or MF"
    echo "  INPUT_SOURCE: filtered (default), combined_annotated, gene_families,"
    echo "                unfiltered, or a full file path"
    echo ""
    echo "Examples:"
    echo "  sbatch scripts/run_pydeseq2_step3_plots.sbatch DC"
    echo "  sbatch scripts/run_pydeseq2_step3_plots.sbatch DC combined_annotated"
    echo "  sbatch scripts/run_pydeseq2_step3_plots.sbatch DC gene_families"
    echo "  sbatch scripts/run_pydeseq2_step3_plots.sbatch DC unfiltered"
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
# RESOLVE INPUT FILE
# ══════════════════════════════════════════════════════════════════════════════

case "${INPUT_SOURCE}" in
    filtered)
        INPUT_FILE="${BASE_DIR}/06_analysis/pydeseq2_${SPECIES}_step2_filtered/pydeseq2_results_FILTERED.tsv"
        SOURCE_LABEL="Step 2 filtered results"
        OUTPUT_DIR="${BASE_DIR}/06_analysis/pydeseq2_${SPECIES}_step3_plots"
        ;;
    combined_annotated)
        INPUT_FILE="${BASE_DIR}/06_analysis/combined_${SPECIES}/${SPECIES}_swissprot_discovery_annotated.tsv"
        SOURCE_LABEL="Combined BLAST+DESeq2 annotated"
        OUTPUT_DIR="${BASE_DIR}/06_analysis/pydeseq2_${SPECIES}_step3_plots_annotated"
        ;;
    gene_families)
        INPUT_FILE="${BASE_DIR}/06_analysis/gene_families_${SPECIES}/${SPECIES}_CYP_OMT_combined.tsv"
        SOURCE_LABEL="Gene family verified (CYP+OMT)"
        OUTPUT_DIR="${BASE_DIR}/06_analysis/pydeseq2_${SPECIES}_step3_plots_families"
        ;;
    unfiltered)
        INPUT_FILE="${BASE_DIR}/06_analysis/pydeseq2_${SPECIES}_step1_unfiltered/pydeseq2_results_UNFILTERED.tsv"
        SOURCE_LABEL="Step 1 unfiltered results"
        OUTPUT_DIR="${BASE_DIR}/06_analysis/pydeseq2_${SPECIES}_step3_plots_unfiltered"
        ;;
    *)
        INPUT_FILE="${INPUT_SOURCE}"
        SOURCE_LABEL="Custom file"
        BASENAME=$(basename "${INPUT_SOURCE}" .tsv)
        OUTPUT_DIR="${BASE_DIR}/06_analysis/pydeseq2_${SPECIES}_step3_plots_${BASENAME}"
        ;;
esac

# Count matrix (from build_count_matrix.sbatch)
case "${SPECIES}" in
    DC) COUNT_FOLDER="00_1_DC" ;;
    DG) COUNT_FOLDER="00_2_DG" ;;
    MF) COUNT_FOLDER="00_3_MF" ;;
esac

COUNT_MATRIX="${BASE_DIR}/03_count_tables/${COUNT_FOLDER}/gene_count_matrix.tsv"
METADATA="${BASE_DIR}/03_count_tables/${COUNT_FOLDER}/sample_metadata.tsv"

# HMMER (auto-detected, optional)
HMMER_FILE="${BASE_DIR}/06_analysis/hmmer_${SPECIES}/all_genes_protein_pfam_domains.txt"


# ══════════════════════════════════════════════════════════════════════════════
# ACTIVATE ENVIRONMENT
# ══════════════════════════════════════════════════════════════════════════════

source ~/.bashrc
conda activate rnaseq

if [ -z "${CONDA_DEFAULT_ENV}" ]; then
    echo "ERROR: Conda environment didn't activate"
    exit 1
fi


# ══════════════════════════════════════════════════════════════════════════════
# VALIDATE
# ══════════════════════════════════════════════════════════════════════════════

echo "============================================================"
echo "PyDESeq2 STEP 3: Generate Plots"
echo "============================================================"
echo ""
echo "Species:      ${SPECIES_NAME} (${SPECIES})"
echo "Input source: ${SOURCE_LABEL}"
echo "Input file:   ${INPUT_FILE}"
echo "Count matrix: ${COUNT_MATRIX}"
echo "Output:       ${OUTPUT_DIR}"
echo ""
echo "Plot settings:"
echo "  Scaling:       ${SCALE_METHOD}"
echo "  Cluster rows:  ${CLUSTER_ROWS}"
echo ""
echo "Start time: $(date)"
echo ""

echo "------------------------------------------------------------"
echo "CHECKING FILES"
echo "------------------------------------------------------------"

if [ ! -f "${INPUT_FILE}" ]; then
    echo "ERROR: Input file not found!"
    echo "  ${INPUT_FILE}"
    echo ""
    case "${INPUT_SOURCE}" in
        filtered)
            echo "Run Step 2 first:"
            echo "  sbatch scripts/run_pydeseq2_step2_filter.sbatch ${SPECIES}"
            ;;
        combined_annotated)
            echo "Run combine+filter first:"
            echo "  sbatch scripts/run_combine_filter.sbatch ${SPECIES} swissprot discovery"
            ;;
        gene_families)
            echo "Run gene family extraction first:"
            echo "  sbatch scripts/run_gene_family_heatmap.sbatch ${SPECIES} swissprot discovery full"
            ;;
        unfiltered)
            echo "Run Step 1 first:"
            echo "  sbatch scripts/run_pydeseq2_step1_analysis.sbatch ${SPECIES}"
            ;;
        *)
            echo "Check the file path and try again."
            ;;
    esac
    exit 1
fi
echo "  Input file: found"

COUNT_ARG=""
META_ARG=""
if [ -f "${COUNT_MATRIX}" ] && [ -f "${METADATA}" ]; then
    echo "  Count matrix: found"
    echo "  Metadata: found"
    COUNT_ARG="--count-matrix \"${COUNT_MATRIX}\""
    META_ARG="--metadata \"${METADATA}\""
else
    echo "  Count matrix / metadata: NOT FOUND (heatmaps will be skipped)"
fi

HMMER_ARG=""
if [ -f "${HMMER_FILE}" ]; then
    echo "  HMMER domains: found (will add domain labels to heatmaps)"
    HMMER_ARG="--hmmer \"${HMMER_FILE}\""
else
    echo "  HMMER domains: not found (skipping domain labels)"
fi

echo ""


# ══════════════════════════════════════════════════════════════════════════════
# GENERATE PLOTS
# ══════════════════════════════════════════════════════════════════════════════

mkdir -p "${OUTPUT_DIR}"

echo "============================================================"
echo "GENERATING PLOTS"
echo "============================================================"
echo ""
echo "Creating:"
echo "  - MA plot (baseMean vs log2FoldChange)"
echo "  - Volcano plot (log2FC vs significance)"
echo "  - CYP heatmap (if CYP genes detected in input)"
echo "  - OMT heatmap (if OMT genes detected in input)"
echo ""

PLOT_CMD="python3 \"${CODE_DIR}/pydeseq2_generate_plots.py\""
PLOT_CMD="${PLOT_CMD} \"${INPUT_FILE}\""
PLOT_CMD="${PLOT_CMD} -o \"${OUTPUT_DIR}\""
PLOT_CMD="${PLOT_CMD} ${COUNT_ARG} ${META_ARG} ${HMMER_ARG}"
PLOT_CMD="${PLOT_CMD} --contrast-factor \"${CONTRAST_FACTOR}\""
PLOT_CMD="${PLOT_CMD} --contrast-A \"${CONTRAST_A}\""
PLOT_CMD="${PLOT_CMD} --contrast-B \"${CONTRAST_B}\""
PLOT_CMD="${PLOT_CMD} --scale ${SCALE_METHOD}"

if [ "${CLUSTER_ROWS}" = "false" ]; then
    PLOT_CMD="${PLOT_CMD} --no-row-cluster"
fi

eval ${PLOT_CMD}

if [ $? -ne 0 ]; then
    echo "ERROR: Plot generation failed"
    exit 1
fi


# ══════════════════════════════════════════════════════════════════════════════
# DONE
# ══════════════════════════════════════════════════════════════════════════════

echo ""
echo "============================================================"
echo "STEP 3 COMPLETE!"
echo "============================================================"
echo ""
echo "Species: ${SPECIES_NAME} (${SPECIES})"
echo "Output:  ${OUTPUT_DIR}"
echo ""
echo "Files:"
ls -la "${OUTPUT_DIR}"/ 2>/dev/null | grep -v "^total" | grep -v "^\." || true
echo ""
echo "------------------------------------------------------------"
echo "DOWNLOAD PLOTS"
echo "------------------------------------------------------------"
echo ""
echo "scp daisycortesj@tinkercliffs.arc.vt.edu:${OUTPUT_DIR}/*.pdf ~/Desktop/"
echo ""
echo "End time: $(date)"
echo "============================================================"

exit 0
