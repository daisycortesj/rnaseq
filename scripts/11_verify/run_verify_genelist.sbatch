#!/bin/bash
################################################################################
# Verify step 3 results against P450_list_RefSeq.txt (TWO RUNS)
################################################################################
#
# This script compares YOUR PyDESeq2 results against a previous student's
# P450 gene data.  It supports two databases:
#
#   DATABASE=sukman  (default)  — Sukman's 5 reference files (.txt)
#   DATABASE=ahmed              — Ahmed's upregulated P450 CSV
#
# Runs verify_genelist.py TWICE to produce two comparison tables:
#   RUN 1: Filtered results (cyp_gene_list.tsv from Step 3 or Path A)
#   RUN 2: Unfiltered results (pydeseq2_results_UNFILTERED.tsv from Step 1)
#
# Compare the two tables to see which P450 genes were filtered out and why
# (padj too high? log2FC too small? not expressed at all?).
#
# ──────────────────────────────────────────────────────────────────────
# HOW TO RUN
# ──────────────────────────────────────────────────────────────────────
#
# Compare DC data against Sukman (the default):
#   cd /projects/tholl_lab_1/daisy_analysis
#   sbatch 05_rnaseq-code/scripts/11_verify/run_verify_genelist.sbatch
#
# CHANGED: Compare DC data against Ahmed instead:
#   DATABASE=ahmed \
#     sbatch 05_rnaseq-code/scripts/11_verify/run_verify_genelist.sbatch
#
# Compare SK data against Sukman:
#   RESULTS=06_analysis/pydeseq2_SK_step3_.../cyp_gene_list.tsv \
#   RESULTS_UNFILTERED=06_analysis/pydeseq2_SK_step1/pydeseq2_results_UNFILTERED.tsv \
#   OUTPUT=07_NRdatabase/sukman_database/verify_SK_filtered_sukman.tsv \
#   OUTPUT_UNFILTERED=07_NRdatabase/sukman_database/verify_SK_unfiltered_sukman.tsv \
#     sbatch 05_rnaseq-code/scripts/11_verify/run_verify_genelist.sbatch
#
# Compare SK data against Ahmed:
#   DATABASE=ahmed \
#   RESULTS=06_analysis/pydeseq2_SK_step3_.../cyp_gene_list.tsv \
#   RESULTS_UNFILTERED=06_analysis/pydeseq2_SK_step1/pydeseq2_results_UNFILTERED.tsv \
#   OUTPUT=07_NRdatabase/ahmed_database/verify_SK_filtered_ahmed.tsv \
#   OUTPUT_UNFILTERED=07_NRdatabase/ahmed_database/verify_SK_unfiltered_ahmed.tsv \
#     sbatch 05_rnaseq-code/scripts/11_verify/run_verify_genelist.sbatch
#
# Skip all previous student comparisons (Check 1 only):
#   PREV_DE=none PREV_EXPR=none PREV_LFC=none PREV_UP=none PREV_DOWN=none \
#     sbatch 05_rnaseq-code/scripts/11_verify/run_verify_genelist.sbatch
#
# ──────────────────────────────────────────────────────────────────────
# SIX CHECKS (for each run)
# ──────────────────────────────────────────────────────────────────────
#
#   CHECK 1: Are all genes in your results in the gene list?
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
#   CONTRAST_A=L CONTRAST_B=R \
#     sbatch 05_rnaseq-code/scripts/11_verify/run_verify_genelist.sbatch
#
# DG vs Ahmed (leaf focus):
#   DATABASE=ahmed CONTRAST_A=L CONTRAST_B=R \
#   RESULTS=06_analysis/pydeseq2_DG_step3_plots_geneious_candidates_DG_padj005_lfc20/cyp_gene_list.tsv \
#   RESULTS_UNFILTERED=06_analysis/pydeseq2_DG_step1_unfiltered/pydeseq2_results_UNFILTERED.tsv \
#   OUTPUT=07_NRdatabase/ahmed_database/verify_DG_filtered_ahmed.tsv \
#   OUTPUT_UNFILTERED=07_NRdatabase/ahmed_database/verify_DG_unfiltered_ahmed.tsv \
#     sbatch 05_rnaseq-code/scripts/11_verify/run_verify_genelist.sbatch
#
# ──────────────────────────────────────────────────────────────────────
# OUTPUT FILES  (names end with which database was compared)
# ──────────────────────────────────────────────────────────────────────
#
#   verify_DC_filtered_sukman.tsv           — filtered results vs Sukman
#   verify_DC_filtered_SUMMARY_sukman.txt   — human-readable summary
#   verify_DC_filtered_ahmed.tsv            — filtered results vs Ahmed
#   verify_DC_filtered_SUMMARY_ahmed.txt    — human-readable summary
#   (same pattern for unfiltered)
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

# ── Contrast direction (default: R vs L, same as DC) ──
# Override for DG: CONTRAST_A=L CONTRAST_B=R
CONTRAST_A="${CONTRAST_A:-R}"
CONTRAST_B="${CONTRAST_B:-L}"

# ── Your results (same regardless of which database you compare to) ──
# Filtered = step 3 / Path A output.  Unfiltered = step 1 raw output.
RESULTS="${RESULTS:-06_analysis/pydeseq2_DC_step3_plots_geneious_candidates_padj005_lfc20/cyp_gene_list.tsv}"
RESULTS_UNFILTERED="${RESULTS_UNFILTERED:-06_analysis/pydeseq2_DC_step1_unfiltered/pydeseq2_results_UNFILTERED.tsv}"

# ============================================================================
# CHANGED: DATABASE SELECTION — pick which previous student to compare against
# ============================================================================
#
# Set DATABASE before running sbatch to choose the comparison database:
#   DATABASE=sukman  (default)  — uses Sukman's separate .txt reference files
#   DATABASE=ahmed              — uses Ahmed's single upregulated P450 CSV
#
# To add a NEW database later, just add another "elif" block below
# following the same pattern.  Copy the "ahmed" block and change:
#   1. The DB_PATH to your new folder
#   2. The file paths (GENE_LIST, PREV_DE, etc.)
#   3. The OUTPUT names (end with your new database name)
# ============================================================================

DATABASE="${DATABASE:-sukman}"

# --------------------------------------------------------------------------
# This if/elif sets all the file paths depending on which DATABASE you chose.
# --------------------------------------------------------------------------

if [ "$DATABASE" = "sukman" ]; then
    # ── Sukman database: 5 separate .txt files ──
    DB_PATH="${BASE_DIR}/07_NRdatabase/sukman_database"

    # Check 1: gene list to verify membership
    GENE_LIST="${GENE_LIST:-${DB_PATH}/P450_list_RefSeq.txt}"

    # Checks 2-6: previous student reference files
    PREV_DE="${PREV_DE:-${DB_PATH}/P450_list_RefSeq_log2fold.txt}"
    PREV_EXPR="${PREV_EXPR:-${DB_PATH}/P450_expression_refseq.txt}"
    PREV_LFC="${PREV_LFC:-${DB_PATH}/P450_expression_refseq_logfold2.txt}"
    PREV_UP="${PREV_UP:-${DB_PATH}/P450_upregulated.txt}"
    PREV_DOWN="${PREV_DOWN:-${DB_PATH}/P450_downregulated_list.txt}"

    # Output tables — filenames end with _sukman so you know the source
    # CHANGED: was "verification_comparison_filtered.tsv"; now ends with
    # _sukman to match the ahmed naming convention.
    OUTPUT="${OUTPUT:-${DB_PATH}/verify_DC_filtered_sukman.tsv}"
    OUTPUT_UNFILTERED="${OUTPUT_UNFILTERED:-${DB_PATH}/verify_DC_unfiltered_sukman.tsv}"

elif [ "$DATABASE" = "ahmed" ]; then
    # ── CHANGED: Ahmed database — single CSV with upregulated P450 genes ──
    # Ahmed's file has full DESeq2 columns (gene_id, baseMean,
    # log2FoldChange, lfcSE, stat, pvalue, padj), so the same CSV
    # works for multiple checks.
    DB_PATH="${BASE_DIR}/07_NRdatabase/ahmed_database"
    AHMED_CSV="${DB_PATH}/Upregulated_P450_values_ahmed.csv"

    # Check 1: use Ahmed's gene IDs as the gene list
    GENE_LIST="${GENE_LIST:-${AHMED_CSV}}"

    # Checks 2-5: the same CSV provides DE list, expression (baseMean),
    # log2fold (log2FoldChange), and upregulated gene IDs.
    PREV_DE="${PREV_DE:-${AHMED_CSV}}"
    PREV_EXPR="${PREV_EXPR:-${AHMED_CSV}}"
    PREV_LFC="${PREV_LFC:-${AHMED_CSV}}"
    PREV_UP="${PREV_UP:-${AHMED_CSV}}"

    # Check 6: Ahmed only has upregulated genes — no downregulated file.
    # Set to empty so Check 6 is skipped automatically.
    PREV_DOWN="${PREV_DOWN:-}"

    # Output tables — filenames end with _ahmed
    OUTPUT="${OUTPUT:-${DB_PATH}/verify_DC_filtered_ahmed.tsv}"
    OUTPUT_UNFILTERED="${OUTPUT_UNFILTERED:-${DB_PATH}/verify_DC_unfiltered_ahmed.tsv}"

else
    # If someone types a database name that doesn't exist, stop with a
    # helpful error message listing the valid options.
    echo "ERROR: Unknown DATABASE='${DATABASE}'"
    echo ""
    echo "  Valid options:"
    echo "    DATABASE=sukman   (Sukman's P450 reference files)"
    echo "    DATABASE=ahmed    (Ahmed's upregulated P450 CSV)"
    echo ""
    echo "  Example:"
    echo "    DATABASE=ahmed sbatch 05_rnaseq-code/scripts/11_verify/run_verify_genelist.sbatch"
    exit 1
fi

# Make paths absolute if they aren't already
[[ "$RESULTS"            != /* ]] && RESULTS="${BASE_DIR}/${RESULTS}"
[[ "$RESULTS_UNFILTERED" != /* ]] && RESULTS_UNFILTERED="${BASE_DIR}/${RESULTS_UNFILTERED}"
[[ "$GENE_LIST"          != /* ]] && GENE_LIST="${BASE_DIR}/${GENE_LIST}"

# ============================================================================
# VALIDATION — make sure all required files exist before running
# ============================================================================
echo "========================================================================"
echo "  Verify Gene List — ${DATABASE} database"
echo "  $(date)"
echo "========================================================================"
echo ""
LABEL_A=$( [ "${CONTRAST_A}" = "R" ] && echo "root" || echo "leaf" )
LABEL_B=$( [ "${CONTRAST_B}" = "L" ] && echo "leaf" || echo "root" )
echo "  Database:               ${DATABASE}"
echo "  Contrast:               ${CONTRAST_A} vs ${CONTRAST_B} (positive log2FC = up in ${LABEL_A})"
echo "  Results (filtered):     ${RESULTS}"
echo "  Results (unfiltered):   ${RESULTS_UNFILTERED}"
echo "  Gene list:              ${GENE_LIST}"
echo "  Prev DE list:           ${PREV_DE:-<none>}"
echo "  Prev expression:        ${PREV_EXPR:-<none>}"
echo "  Prev log2fold:          ${PREV_LFC:-<none>}"
echo "  Prev upregulated:       ${PREV_UP:-<none>}"
echo "  Prev downreg:           ${PREV_DOWN:-<none — skipped>}"
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
# BUILD OPTIONAL ARGUMENTS — skip checks whose files are missing or empty
# ============================================================================

PREV_DE_ARG=""
PREV_EXPR_ARG=""
PREV_LFC_ARG=""
PREV_UP_ARG=""
PREV_DOWN_ARG=""

# For each check: if the variable is non-empty AND the file exists, pass it.
# Otherwise print a note and skip that check.

if [ -n "${PREV_DE}" ] && [ -f "$PREV_DE" ]; then
    PREV_DE_ARG="--prev-de ${PREV_DE}"
else
    echo "  NOTE: Skipping Check 2 (DE overlap) — no file"
fi

if [ -n "${PREV_EXPR}" ] && [ -f "$PREV_EXPR" ]; then
    PREV_EXPR_ARG="--prev-expression ${PREV_EXPR}"
else
    echo "  NOTE: Skipping Check 3 (expression) — no file"
fi

if [ -n "${PREV_LFC}" ] && [ -f "$PREV_LFC" ]; then
    PREV_LFC_ARG="--prev-log2fold ${PREV_LFC}"
else
    echo "  NOTE: Skipping Check 4 (log2fold) — no file"
fi

if [ -n "${PREV_UP}" ] && [ -f "$PREV_UP" ]; then
    PREV_UP_ARG="--prev-upregulated ${PREV_UP}"
else
    echo "  NOTE: Skipping Check 5 (upregulated) — no file"
fi

# CHANGED: added "-n" check so an empty PREV_DOWN (Ahmed) is skipped cleanly
# instead of testing -f on an empty string.
if [ -n "${PREV_DOWN}" ] && [ -f "$PREV_DOWN" ]; then
    PREV_DOWN_ARG="--prev-downregulated ${PREV_DOWN}"
else
    echo "  NOTE: Skipping Check 6 (downregulated) — no file"
fi
echo ""

# ────────────────────────────────────────────────────────────────────────────
# RUN 1: Filtered results (cyp_gene_list.tsv from Step 3 or Path A)
# ────────────────────────────────────────────────────────────────────────────
echo "╔══════════════════════════════════════════════════════════════════════╗"
echo "║  RUN 1 OF 2: FILTERED RESULTS  (vs ${DATABASE} database)"
echo "║  ${RESULTS}"
echo "╚══════════════════════════════════════════════════════════════════════╝"
echo ""

# CHANGED: added --database-label so the summary file ends with
# _SUMMARY_sukman.txt or _SUMMARY_ahmed.txt (tells you which database).
python "${CODE_DIR}/scripts/11_verify/verify_genelist.py" \
    --results "$RESULTS" \
    --gene-list "$GENE_LIST" \
    --output "$OUTPUT" \
    --contrast-a "$CONTRAST_A" \
    --contrast-b "$CONTRAST_B" \
    --database-label "$DATABASE" \
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
echo "║  RUN 2 OF 2: UNFILTERED RESULTS (Step 1)  (vs ${DATABASE} database)"
echo "║  ${RESULTS_UNFILTERED}"
echo "╚══════════════════════════════════════════════════════════════════════╝"
echo ""

python "${CODE_DIR}/scripts/11_verify/verify_genelist.py" \
    --results "$RESULTS_UNFILTERED" \
    --gene-list "$GENE_LIST" \
    --output "$OUTPUT_UNFILTERED" \
    --contrast-a "$CONTRAST_A" \
    --contrast-b "$CONTRAST_B" \
    --database-label "$DATABASE" \
    $PREV_DE_ARG \
    $PREV_EXPR_ARG \
    $PREV_LFC_ARG \
    $PREV_UP_ARG \
    $PREV_DOWN_ARG

echo ""
echo "========================================================================"
echo "  BOTH RUNS COMPLETE — compared against ${DATABASE} database"
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
