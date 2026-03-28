#!/bin/bash
################################################################################
# Verify step 3 results against P450_list_RefSeq.txt (TWO RUNS)
################################################################################
#
# Runs verify_genelist.py TWICE to produce two comparison tables:
#   RUN 1: Filtered results (cyp_gene_list.tsv from Step 3 or Path A)
#   RUN 2: Unfiltered results (pydeseq2_results_UNFILTERED.tsv from Step 1)
#
# Compare the two tables to see which P450 genes were filtered out and why
# (padj too high? log2FC too small? not expressed at all?).
#
# How to run (defaults check DC data):
#   cd /projects/tholl_lab_1/daisy_analysis
#   sbatch 05_rnaseq-code/scripts/11_verify/run_verify_genelist.sbatch
#
# Run for SK:
#   RESULTS=06_analysis/pydeseq2_SK_step3_.../cyp_gene_list.tsv \
#   RESULTS_UNFILTERED=06_analysis/pydeseq2_SK_step1/pydeseq2_results_UNFILTERED.tsv \
#   OUTPUT=06_analysis/verify_SK_filtered.tsv \
#   OUTPUT_UNFILTERED=06_analysis/verify_SK_unfiltered.tsv \
#     sbatch 05_rnaseq-code/scripts/11_verify/run_verify_genelist.sbatch
#
# Skip the previous student comparisons (Check 1 only):
#   PREV_DE=none PREV_EXPR=none PREV_LFC=none PREV_UP=none PREV_DOWN=none \
#     sbatch 05_rnaseq-code/scripts/11_verify/run_verify_genelist.sbatch
#
# Six checks are performed (for each run):
#   CHECK 1: Are all genes in your results actually in P450_list_RefSeq.txt?
#   CHECK 2: DE overlap — which of previous student's DE genes appear in yours?
#            (P450_list_RefSeq_log2fold.txt)
#   CHECK 3: Expression comparison — correlate expression levels with
#            previous student's data (P450_expression_refseq.txt)
#   CHECK 4: Log2fold direction — do your genes go up/down in the same
#            direction as the previous student? Correlation of magnitudes?
#            (P450_expression_refseq_logfold2.txt)
#   CHECK 5: Upregulated agreement — are genes the previous student found
#            as upregulated also upregulated in your data?
#            (P450_upregulated.txt — note: IDs lack LOC prefix, handled automatically)
#   CHECK 6: Downregulated agreement — are genes the previous student found
#            as downregulated also downregulated in your data?
#            (P450_downregulated_list.txt — note: IDs may lack LOC prefix, handled automatically)
#
# Default contrast: R vs L (positive log2FC = up in Root)
# Override for DG (leaf focus):
#   CONTRAST_A=L CONTRAST_B=R RESULTS=... RESULTS_UNFILTERED=... \
#     sbatch 05_rnaseq-code/scripts/11_verify/run_verify_genelist.sbatch
#
# Output (two tables):
#   - verification_comparison_filtered.tsv   — your filtered results vs gene list
#   - verification_comparison_unfiltered.tsv — Step 1 unfiltered results vs gene list
#   Open both in Excel and compare to see what's missing!
#
################################################################################

#SBATCH --job-name=verify_genes
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu

#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis
#SBATCH -o 06_analysis/verify_genelist_%j.out
#SBATCH -e 06_analysis/verify_genelist_%j.err

set -euo pipefail

source ~/.bashrc
conda activate rnaseq

# ============================================================================
# PATHS — all relative paths are resolved against BASE_DIR
# ============================================================================
# ===== Load shared config =====
source "${SLURM_SUBMIT_DIR:-.}/scripts/config.sh"

SUKMAN="${BASE_DIR}/07_NRdatabase/sukman_database"

# Contrast direction (default: R vs L, same as DC)
# Override for DG: CONTRAST_A=L CONTRAST_B=R
CONTRAST_A="${CONTRAST_A:-R}"
CONTRAST_B="${CONTRAST_B:-L}"

GENE_LIST="${GENE_LIST:-${SUKMAN}/P450_list_RefSeq.txt}"

# Two results paths: filtered (from step 3 or Path A) and unfiltered (from step 1)
RESULTS="${RESULTS:-06_analysis/pydeseq2_DC_step3_plots_geneious_candidates_padj005_lfc20/cyp_gene_list.tsv}"
RESULTS_UNFILTERED="${RESULTS_UNFILTERED:-06_analysis/pydeseq2_DC_step1_unfiltered/pydeseq2_results_UNFILTERED.tsv}"

PREV_DE="${PREV_DE:-${SUKMAN}/P450_list_RefSeq_log2fold.txt}"
PREV_EXPR="${PREV_EXPR:-${SUKMAN}/P450_expression_refseq.txt}"
PREV_LFC="${PREV_LFC:-${SUKMAN}/P450_expression_refseq_logfold2.txt}"
PREV_UP="${PREV_UP:-${SUKMAN}/P450_upregulated.txt}"
PREV_DOWN="${PREV_DOWN:-${SUKMAN}/P450_downregulated_list.txt}"

# Two output tables — one for each results file
OUTPUT="${OUTPUT:-${SUKMAN}/verification_comparison_filtered.tsv}"
OUTPUT_UNFILTERED="${OUTPUT_UNFILTERED:-${SUKMAN}/verification_comparison_unfiltered.tsv}"

# Make paths absolute if they aren't already
[[ "$RESULTS"            != /* ]] && RESULTS="${BASE_DIR}/${RESULTS}"
[[ "$RESULTS_UNFILTERED" != /* ]] && RESULTS_UNFILTERED="${BASE_DIR}/${RESULTS_UNFILTERED}"
[[ "$GENE_LIST"          != /* ]] && GENE_LIST="${BASE_DIR}/${GENE_LIST}"

# ============================================================================
# VALIDATION
# ============================================================================
echo "========================================================================"
echo "  Verify Gene List — Double Check"
echo "  $(date)"
echo "========================================================================"
echo ""
LABEL_A=$( [ "${CONTRAST_A}" = "R" ] && echo "root" || echo "leaf" )
LABEL_B=$( [ "${CONTRAST_B}" = "L" ] && echo "leaf" || echo "root" )
echo "  Contrast:               ${CONTRAST_A} vs ${CONTRAST_B} (positive log2FC = up in ${LABEL_A})"
echo "  Results (filtered):     ${RESULTS}"
echo "  Results (unfiltered):   ${RESULTS_UNFILTERED}"
echo "  Gene list:              ${GENE_LIST}"
echo "  Prev DE list:           ${PREV_DE}"
echo "  Prev expression:        ${PREV_EXPR}"
echo "  Prev log2fold:          ${PREV_LFC}"
echo "  Prev upregulated:       ${PREV_UP}"
echo "  Prev downreg:           ${PREV_DOWN}"
echo "  Output (filtered):      ${OUTPUT}"
echo "  Output (unfiltered):    ${OUTPUT_UNFILTERED}"
echo ""

if [ ! -f "$GENE_LIST" ]; then
    echo "ERROR: Gene list not found: $GENE_LIST"
    echo ""
    echo "  Check the file exists:"
    echo "    ls -la ${GENE_LIST}"
    exit 1
fi

if [ ! -f "$RESULTS" ]; then
    echo "ERROR: Filtered results file not found: $RESULTS"
    echo ""
    echo "  Did step 3 / Path A finish? Check with:"
    echo "    ls -la ${BASE_DIR}/06_analysis/pydeseq2_*_step3_plots_*/"
    exit 1
fi

if [ ! -f "$RESULTS_UNFILTERED" ]; then
    echo "ERROR: Unfiltered results file not found: $RESULTS_UNFILTERED"
    echo ""
    echo "  Did PyDESeq2 step 1 finish? Check with:"
    echo "    ls -la ${BASE_DIR}/06_analysis/pydeseq2_*_step1/"
    exit 1
fi

# ============================================================================
# RUN
# ============================================================================

PREV_DE_ARG=""
PREV_EXPR_ARG=""
PREV_LFC_ARG=""
PREV_UP_ARG=""
PREV_DOWN_ARG=""

if [ -f "$PREV_DE" ]; then
    PREV_DE_ARG="--prev-de ${PREV_DE}"
else
    echo "  NOTE: ${PREV_DE} not found, skipping Check 2"
fi

if [ -f "$PREV_EXPR" ]; then
    PREV_EXPR_ARG="--prev-expression ${PREV_EXPR}"
else
    echo "  NOTE: ${PREV_EXPR} not found, skipping Check 3"
fi

if [ -f "$PREV_LFC" ]; then
    PREV_LFC_ARG="--prev-log2fold ${PREV_LFC}"
else
    echo "  NOTE: ${PREV_LFC} not found, skipping Check 4"
fi

if [ -f "$PREV_UP" ]; then
    PREV_UP_ARG="--prev-upregulated ${PREV_UP}"
else
    echo "  NOTE: ${PREV_UP} not found, skipping Check 5"
fi

if [ -f "$PREV_DOWN" ]; then
    PREV_DOWN_ARG="--prev-downregulated ${PREV_DOWN}"
else
    echo "  NOTE: ${PREV_DOWN} not found, skipping Check 6"
fi
echo ""

# ────────────────────────────────────────────────────────────────────────────
# RUN 1: Filtered results (cyp_gene_list.tsv from Step 3 or Path A)
# ────────────────────────────────────────────────────────────────────────────
echo "╔══════════════════════════════════════════════════════════════════════╗"
echo "║  RUN 1 OF 2: FILTERED RESULTS                                      ║"
echo "║  ${RESULTS}"
echo "╚══════════════════════════════════════════════════════════════════════╝"
echo ""

python "${CODE_DIR}/scripts/11_verify/verify_genelist.py" \
    --results "$RESULTS" \
    --gene-list "$GENE_LIST" \
    --output "$OUTPUT" \
    $PREV_DE_ARG \
    $PREV_EXPR_ARG \
    $PREV_LFC_ARG \
    $PREV_UP_ARG \
    $PREV_DOWN_ARG

echo ""

# ────────────────────────────────────────────────────────────────────────────
# RUN 2: Unfiltered results (pydeseq2_results_UNFILTERED.tsv from Step 1)
# ────────────────────────────────────────────────────────────────────────────
echo "╔══════════════════════════════════════════════════════════════════════╗"
echo "║  RUN 2 OF 2: UNFILTERED RESULTS (Step 1)                           ║"
echo "║  ${RESULTS_UNFILTERED}"
echo "╚══════════════════════════════════════════════════════════════════════╝"
echo ""

python "${CODE_DIR}/scripts/11_verify/verify_genelist.py" \
    --results "$RESULTS_UNFILTERED" \
    --gene-list "$GENE_LIST" \
    --output "$OUTPUT_UNFILTERED" \
    $PREV_DE_ARG \
    $PREV_EXPR_ARG \
    $PREV_LFC_ARG \
    $PREV_UP_ARG \
    $PREV_DOWN_ARG

echo ""
echo "========================================================================"
echo "  BOTH RUNS COMPLETE"
echo "========================================================================"
echo ""
echo "  Table 1 (filtered):    ${OUTPUT}"
echo "  Table 2 (unfiltered):  ${OUTPUT_UNFILTERED}"
echo ""
echo "  Download both to your laptop:"
echo "    scp daisycortesj@tinkercliffs.arc.vt.edu:${OUTPUT} ~/Desktop/"
echo "    scp daisycortesj@tinkercliffs.arc.vt.edu:${OUTPUT_UNFILTERED} ~/Desktop/"
echo ""
echo "  Compare the two tables in Excel to see what's missing!"
echo ""
echo "Done: $(date)"
