#!/bin/bash
################################################################################
# Filter count matrix by gene list + run PyDESeq2 for expression stats
################################################################################
#
# What this does:
#   1. Loads the count matrix (all ~33K genes)
#   2. Runs PyDESeq2 on ALL genes (uses full dataset for normalization)
#   3. Filters to only genes in your gene list (e.g. P450_list_RefSeq.txt)
#   4. Applies DE cutoffs (padj < 0.05, |log2FC| > 2.0)
#   5. Outputs curated candidate list → ready for step 3 plots
#
# How to run — Geneious P450 list:
#   sbatch scripts/run_filter_count_genelist.sbatch
#
# How to run — your CYP master list:
#   GENE_LIST=07_NRdatabase/cyp450_database/cyp_master_list.csv \
#   OUT_NAME=cyp_candidates.tsv \
#     sbatch scripts/run_filter_count_genelist.sbatch
#
# How to run — any gene list:
#   GENE_LIST=/path/to/gene_list.txt OUT_NAME=my_candidates.tsv \
#     sbatch scripts/run_filter_count_genelist.sbatch
#
# Custom cutoffs:
#   PADJ=0.01 LFC=3.0 sbatch scripts/run_filter_count_genelist.sbatch
#
# Match previous student's min-count threshold (20 instead of default 10):
#   MIN_COUNTS=20 OUT_NAME=geneious_candidates_mincount20.tsv \
#     sbatch scripts/run_filter_count_genelist.sbatch
#
# Output:
#   07_NRdatabase/sukman_database/geneious_candidates.tsv  (default)
#   → ready for: sbatch scripts/run_pydeseq2_step3_plots.sbatch DC <output>
#
################################################################################

#SBATCH --job-name=filter_genelist
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu

#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis
#SBATCH -o 07_NRdatabase/filter_genelist_%j.out
#SBATCH -e 07_NRdatabase/filter_genelist_%j.err

set -euo pipefail

# ===== Load shared config =====
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

source ~/.bashrc
conda activate rnaseq

# ============================================================================
# PATHS — override with env vars if needed
# ============================================================================
COUNT_MATRIX="${BASE_DIR}/03_count_tables/00_1_DC/gene_count_matrix.tsv"
METADATA="${BASE_DIR}/03_count_tables/00_1_DC/sample_metadata.tsv"

GENE_LIST="${GENE_LIST:-${BASE_DIR}/07_NRdatabase/sukman_database/P450_list_RefSeq.txt}"
OUT_DIR="${BASE_DIR}/07_NRdatabase/sukman_database"
OUT_NAME="${OUT_NAME:-geneious_candidates.tsv}"
OUTPUT="${OUT_DIR}/${OUT_NAME}"

PADJ="${PADJ:-0.05}"
LFC="${LFC:-2.0}"
CONTRAST_A="${CONTRAST_A:-R}"
CONTRAST_B="${CONTRAST_B:-L}"
MIN_COUNTS="${MIN_COUNTS:-10}"

# ============================================================================
# VALIDATION
# ============================================================================
echo "========================================================================"
echo "  Filter count matrix by gene list + PyDESeq2"
echo "  Started: $(date)"
echo "========================================================================"
echo ""
echo "  Count matrix:  ${COUNT_MATRIX}"
echo "  Metadata:      ${METADATA}"
echo "  Gene list:     ${GENE_LIST}"
echo "  DE cutoffs:    padj < ${PADJ}, |log2FC| > ${LFC}"
echo "  Min counts:    ${MIN_COUNTS}  (prev student used 20)"
echo "  Contrast:      ${CONTRAST_A} vs ${CONTRAST_B}"
echo "  Output:        ${OUTPUT}"
echo ""

for f in "$COUNT_MATRIX" "$METADATA" "$GENE_LIST"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: File not found: $f"
        exit 1
    fi
done

mkdir -p "$OUT_DIR"

# ============================================================================
# RUN
# ============================================================================

python "${CODE_DIR}/scripts/08_gene_families/filter_count_by_genelist.py" \
    --counts "$COUNT_MATRIX" \
    --gene-list "$GENE_LIST" \
    --metadata "$METADATA" \
    --contrast-a "$CONTRAST_A" \
    --contrast-b "$CONTRAST_B" \
    --min-counts "$MIN_COUNTS" \
    --padj-cutoff "$PADJ" \
    --lfc-cutoff "$LFC" \
    --output "$OUTPUT"

# ============================================================================
# SUMMARY
# ============================================================================
echo ""
echo "========================================================================"
echo "  DONE! $(date)"
echo ""
echo "  Output: ${OUTPUT}"
if [ -f "$OUTPUT" ]; then
    N_GENES=$(tail -n +2 "$OUTPUT" | wc -l)
    echo "  DE candidates: ${N_GENES} genes"
fi
echo ""
echo "  NEXT — generate plots:"
echo "    sbatch scripts/run_pydeseq2_step3_plots.sbatch DC ${OUTPUT}"
echo "========================================================================"
