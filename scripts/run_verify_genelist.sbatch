#!/bin/bash
################################################################################
# Verify step 3 results against P450_list_RefSeq.txt
################################################################################
#
# Double-checks that every gene in your step 3 output actually appears in
# the original P450 gene list. Run AFTER step 3 plots finish.
#
# How to run (defaults check step 3 geneious output):
#   cd /projects/tholl_lab_1/daisy_analysis
#   sbatch 05_rnaseq-code/scripts/run_verify_genelist.sbatch
#
# Check geneious_candidates.tsv (Path A output, before step 3):
#   RESULTS=07_NRdatabase/sukman_database/geneious_candidates.tsv \
#     sbatch 05_rnaseq-code/scripts/run_verify_genelist.sbatch
#
# Check cyp_strict step 3 output:
#   RESULTS=06_analysis/pydeseq2_DC_step3_plots_cyp_strict_padj005_lfc20/cyp_gene_list.tsv \
#     sbatch 05_rnaseq-code/scripts/run_verify_genelist.sbatch
#
# Skip the previous student comparisons (Check 1 only):
#   PREV_DE=none PREV_EXPR=none PREV_LFC=none PREV_UP=none \
#     sbatch 05_rnaseq-code/scripts/run_verify_genelist.sbatch
#
# Five checks are performed:
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
#
# Output:
#   - PASS/FAIL printed to job log
#   - verification_comparison.tsv — side-by-side table with all genes,
#     your stats + prev student stats + flags (open in Excel)
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
BASE_DIR="/projects/tholl_lab_1/daisy_analysis"
CODE_DIR="${BASE_DIR}/05_rnaseq-code"

SUKMAN="${BASE_DIR}/07_NRdatabase/sukman_database"

GENE_LIST="${GENE_LIST:-${SUKMAN}/P450_list_RefSeq.txt}"
RESULTS="${RESULTS:-06_analysis/pydeseq2_DC_step3_plots_geneious_candidates_padj005_lfc20/cyp_gene_list.tsv}"
PREV_DE="${PREV_DE:-${SUKMAN}/P450_list_RefSeq_log2fold.txt}"
PREV_EXPR="${PREV_EXPR:-${SUKMAN}/P450_expression_refseq.txt}"
PREV_LFC="${PREV_LFC:-${SUKMAN}/P450_expression_refseq_logfold2.txt}"
PREV_UP="${PREV_UP:-${SUKMAN}/P450_upregulated.txt}"
OUTPUT="${OUTPUT:-${SUKMAN}/verification_comparison.tsv}"

# Make paths absolute if they aren't already
[[ "$RESULTS"   != /* ]] && RESULTS="${BASE_DIR}/${RESULTS}"
[[ "$GENE_LIST" != /* ]] && GENE_LIST="${BASE_DIR}/${GENE_LIST}"

# ============================================================================
# VALIDATION
# ============================================================================
echo "========================================================================"
echo "  Verify Gene List — Double Check"
echo "  $(date)"
echo "========================================================================"
echo ""
echo "  Results file:     ${RESULTS}"
echo "  Gene list:        ${GENE_LIST}"
echo "  Prev DE list:     ${PREV_DE}"
echo "  Prev expression:  ${PREV_EXPR}"
echo "  Prev log2fold:    ${PREV_LFC}"
echo "  Prev upregulated: ${PREV_UP}"
echo "  Output table:     ${OUTPUT}"
echo ""

if [ ! -f "$GENE_LIST" ]; then
    echo "ERROR: Gene list not found: $GENE_LIST"
    echo ""
    echo "  Check the file exists:"
    echo "    ls -la ${GENE_LIST}"
    exit 1
fi

if [ ! -f "$RESULTS" ]; then
    echo "ERROR: Results file not found: $RESULTS"
    echo ""
    echo "  Did step 3 finish? Check with:"
    echo "    ls -la ${BASE_DIR}/06_analysis/pydeseq2_DC_step3_plots_*/"
    exit 1
fi

# ============================================================================
# RUN
# ============================================================================

PREV_DE_ARG=""
PREV_EXPR_ARG=""
PREV_LFC_ARG=""
PREV_UP_ARG=""

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
echo ""

python "${CODE_DIR}/scripts/verify_genelist.py" \
    --results "$RESULTS" \
    --gene-list "$GENE_LIST" \
    --output "$OUTPUT" \
    $PREV_DE_ARG \
    $PREV_EXPR_ARG \
    $PREV_LFC_ARG \
    $PREV_UP_ARG

echo ""
echo "Output table: ${OUTPUT}"
echo "  Download to your laptop:"
echo "    scp daisycortesj@tinkercliffs.arc.vt.edu:${OUTPUT} ~/Desktop/"
echo ""
echo "Done: $(date)"
