#!/bin/bash
################################################################################
# PYDESEQ2 ANALYSIS - STEP 3: Generate Plots
################################################################################
#
# What this does:
#   Creates publication-quality plots from DESeq2 results:
#   - MA plot (red/blue up/down, gene counts, threshold caption)
#   - Enhanced Volcano plot (4-color, boxed labels with connectors)
#   - PCA plot (top 500 variable genes, 95% confidence ellipses)
#   - Sample correlation heatmap (short names: DC1L, DC1R, etc.)
#   - CYP heatmap (gene locus labels, subfamily brackets: CYP71D, CYP81, etc.)
#   - OMT heatmap (gene locus labels, subfamily brackets: COMT, CCoAOMT, etc.)
#
#   CYP/OMT genes are identified using BLAST descriptions AND HMMER
#   Pfam domains (combined evidence). A gene is included if either
#   source identifies it. Only DE genes (padj + |log2FC| cutoffs) are plotted.
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
#   # Combined two-species heatmap (DC + DG on one plot):
#   sbatch scripts/run_pydeseq2_step3_plots.sbatch DC combined_annotated DG
#
#   # From a specific file:
#   sbatch scripts/run_pydeseq2_step3_plots.sbatch DC /path/to/results.tsv
#
#   # Custom filter cutoffs (override defaults):
#   PADJ_CUTOFF=0.03 LFC_CUTOFF=1.0 sbatch scripts/run_pydeseq2_step3_plots.sbatch DC
#
# What you get:
#   ma_plot.pdf                          - MA plot (red/blue up/down, gene counts)
#   volcano_plot.pdf                     - Enhanced Volcano (4-color, boxed labels)
#   pca_plot.pdf                         - PCA (top 500 genes, 95% ellipses)
#   sample_correlation_heatmap.pdf       - Sample distance heatmap
#   cyp_heatmap.pdf                      - CYP heatmap (locus labels, subfamily brackets)
#   omt_heatmap.pdf                      - OMT heatmap (locus labels, subfamily brackets)
#   cyp_heatmap_combined.pdf             - Combined CYP heatmap (if 2nd species given)
#   omt_heatmap_combined.pdf             - Combined OMT heatmap (if 2nd species given)
#   cyp_gene_list.tsv                    - List of detected CYP genes
#   omt_gene_list.tsv                    - List of detected OMT genes
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

SCALE_METHOD="zscore"
CLUSTER_ROWS=true
CONTRAST_FACTOR="condition"
CONTRAST_A="R"
CONTRAST_B="L"
TOP_N=""       # DE filter is the cutoff; set to a number (e.g. 50) to cap heatmap rows
PADJ_CUTOFF="${PADJ_CUTOFF:-0.05}"   # Only show genes with padj < this in heatmaps
LFC_CUTOFF="${LFC_CUTOFF:-2.0}"     # Only show genes with |log2FC| > this in heatmaps


# ══════════════════════════════════════════════════════════════════════════════
# ARGUMENTS
# ══════════════════════════════════════════════════════════════════════════════

SPECIES="${1:-}"
INPUT_SOURCE="${2:-filtered}"
SPECIES2="${3:-}"

if [ -z "${SPECIES}" ]; then
    echo "ERROR: Missing species argument"
    echo ""
    echo "Usage: sbatch scripts/run_pydeseq2_step3_plots.sbatch SPECIES [INPUT_SOURCE] [SPECIES2]"
    echo ""
    echo "  SPECIES:      DC, DG, or MF"
    echo "  INPUT_SOURCE: filtered (default), combined_annotated, gene_families,"
    echo "                unfiltered, or a full file path"
    echo "  SPECIES2:     Optional second species for combined heatmap (DC, DG, or MF)"
    echo ""
    echo "Examples:"
    echo "  sbatch scripts/run_pydeseq2_step3_plots.sbatch DC"
    echo "  sbatch scripts/run_pydeseq2_step3_plots.sbatch DC combined_annotated"
    echo "  sbatch scripts/run_pydeseq2_step3_plots.sbatch DC combined_annotated DG"
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

# Build a filter tag for the output directory so different cutoffs don't overwrite each other
# e.g. "padj05_lfc2" for defaults, "padj03_lfc1" for 0.03/1.0
FILTER_TAG="padj$(echo ${PADJ_CUTOFF} | tr -d '.')_lfc$(echo ${LFC_CUTOFF} | tr -d '.')"

case "${INPUT_SOURCE}" in
    filtered)
        INPUT_FILE="${BASE_DIR}/06_analysis/pydeseq2_${SPECIES}_step2_filtered/pydeseq2_results_FILTERED.tsv"
        SOURCE_LABEL="Step 2 filtered results"
        OUTPUT_DIR="${BASE_DIR}/06_analysis/pydeseq2_${SPECIES}_step3_plots_${FILTER_TAG}"
        ;;
    combined_annotated)
        INPUT_FILE="${BASE_DIR}/06_analysis/combined_${SPECIES}/${SPECIES}_swissprot_discovery_annotated.tsv"
        SOURCE_LABEL="Combined BLAST+DESeq2 annotated"
        OUTPUT_DIR="${BASE_DIR}/06_analysis/pydeseq2_${SPECIES}_step3_plots_annotated_${FILTER_TAG}"
        ;;
    gene_families)
        INPUT_FILE="${BASE_DIR}/06_analysis/gene_families_${SPECIES}/${SPECIES}_CYP_OMT_combined.tsv"
        SOURCE_LABEL="Gene family verified (CYP+OMT)"
        OUTPUT_DIR="${BASE_DIR}/06_analysis/pydeseq2_${SPECIES}_step3_plots_families_${FILTER_TAG}"
        ;;
    unfiltered)
        INPUT_FILE="${BASE_DIR}/06_analysis/pydeseq2_${SPECIES}_step1_unfiltered/pydeseq2_results_UNFILTERED.tsv"
        SOURCE_LABEL="Step 1 unfiltered results"
        OUTPUT_DIR="${BASE_DIR}/06_analysis/pydeseq2_${SPECIES}_step3_plots_unfiltered_${FILTER_TAG}"
        ;;
    *)
        INPUT_FILE="${INPUT_SOURCE}"
        SOURCE_LABEL="Custom file"
        BASENAME=$(basename "${INPUT_SOURCE}" .tsv)
        OUTPUT_DIR="${BASE_DIR}/06_analysis/pydeseq2_${SPECIES}_step3_plots_${BASENAME}_${FILTER_TAG}"
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

# Second species (optional, for combined heatmap)
COUNT_MATRIX2=""
METADATA2=""
SP2_LABEL=""

if [ -n "${SPECIES2}" ]; then
    SPECIES2=$(echo "${SPECIES2}" | tr '[:lower:]' '[:upper:]')
    if [[ ! "${SPECIES2}" =~ ^(DC|DG|MF)$ ]]; then
        echo "ERROR: Unknown second species: ${SPECIES2}"
        exit 1
    fi
    case "${SPECIES2}" in
        DC) COUNT_FOLDER2="00_1_DC" ;;
        DG) COUNT_FOLDER2="00_2_DG" ;;
        MF) COUNT_FOLDER2="00_3_MF" ;;
    esac
    COUNT_MATRIX2="${BASE_DIR}/03_count_tables/${COUNT_FOLDER2}/gene_count_matrix.tsv"
    METADATA2="${BASE_DIR}/03_count_tables/${COUNT_FOLDER2}/sample_metadata.tsv"
    SP2_LABEL="${SPECIES2}"
fi


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
echo "  Heatmap DE filter: padj < ${PADJ_CUTOFF}, |log2FC| > ${LFC_CUTOFF}"
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
    echo "  HMMER domains: found (domain info saved to gene list TSVs)"
    HMMER_ARG="--hmmer \"${HMMER_FILE}\""
else
    echo "  HMMER domains: not found (skipping)"
fi

# GTF annotation (auto-detected, optional)
GTF_FILE=""
for candidate in \
    "${BASE_DIR}/04_reference/${SPECIES,,}_genomic.gtf" \
    "${BASE_DIR}/04_reference/${SPECIES}_genomic.gtf" \
    "${BASE_DIR}/04_reference/genomic.gtf"; do
    if [ -f "${candidate}" ]; then
        GTF_FILE="${candidate}"
        break
    fi
done

GTF_ARG=""
if [ -n "${GTF_FILE}" ]; then
    echo "  GTF annotation: found"
    GTF_ARG="--gtf \"${GTF_FILE}\""
else
    echo "  GTF annotation: not found (skipping)"
fi

# Check second species files (combined heatmap mode)
COMBINED_ARG=""
if [ -n "${SPECIES2}" ]; then
    echo ""
    echo "  Combined heatmap mode: ${SPECIES} + ${SPECIES2}"
    if [ -f "${COUNT_MATRIX2}" ] && [ -f "${METADATA2}" ]; then
        echo "  ${SPECIES2} count matrix: found"
        echo "  ${SPECIES2} metadata: found"
        COMBINED_ARG="--count-matrix2 \"${COUNT_MATRIX2}\" --metadata2 \"${METADATA2}\" --species1 \"${SPECIES}\" --species2 \"${SPECIES2}\""
    else
        echo "  WARNING: ${SPECIES2} count matrix or metadata not found"
        echo "    Skipping combined heatmaps"
    fi
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
echo "  - MA plot (red/blue up/down, gene counts)"
echo "  - Enhanced Volcano plot (4-color, boxed labels)"
echo "  - PCA plot (top 500 variable genes, 95% ellipses)"
echo "  - Sample correlation heatmap (short names)"
echo "  - CYP heatmap (if CYP genes detected in input)"
echo "  - OMT heatmap (if OMT genes detected in input)"
if [ -n "${COMBINED_ARG}" ]; then
    echo "  - CYP combined heatmap (${SPECIES} + ${SPECIES2})"
    echo "  - OMT combined heatmap (${SPECIES} + ${SPECIES2})"
fi
echo ""

PLOT_CMD="python3 \"${CODE_DIR}/pydeseq2_generate_plots.py\""
PLOT_CMD="${PLOT_CMD} \"${INPUT_FILE}\""
PLOT_CMD="${PLOT_CMD} -o \"${OUTPUT_DIR}\""
PLOT_CMD="${PLOT_CMD} ${COUNT_ARG} ${META_ARG} ${HMMER_ARG} ${GTF_ARG}"
PLOT_CMD="${PLOT_CMD} ${COMBINED_ARG}"
PLOT_CMD="${PLOT_CMD} --contrast-factor \"${CONTRAST_FACTOR}\""
PLOT_CMD="${PLOT_CMD} --contrast-A \"${CONTRAST_A}\""
PLOT_CMD="${PLOT_CMD} --contrast-B \"${CONTRAST_B}\""
PLOT_CMD="${PLOT_CMD} --scale ${SCALE_METHOD}"

if [ "${CLUSTER_ROWS}" = "false" ]; then
    PLOT_CMD="${PLOT_CMD} --no-row-cluster"
fi

if [ -n "${TOP_N}" ]; then
    PLOT_CMD="${PLOT_CMD} --top-n ${TOP_N}"
fi

PLOT_CMD="${PLOT_CMD} --padj-cutoff ${PADJ_CUTOFF}"
PLOT_CMD="${PLOT_CMD} --lfc-cutoff ${LFC_CUTOFF}"

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
