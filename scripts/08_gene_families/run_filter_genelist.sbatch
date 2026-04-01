#!/bin/bash
################################################################################
# Filter count matrix by gene list + run PyDESeq2 for expression stats
################################################################################
#
# What this does:
#   1. Loads the count matrix (all ~33K genes)
#   2. Runs PyDESeq2 on ALL genes (uses full dataset for normalization)
#   3. Filters to only genes in your gene list (CYP or OMT)
#   4. Applies DE cutoffs (padj ≤ 0.05, |log2FC| ≥ 2.0)
#   5. Tags genes with a gene_family label (so heatmaps work without BLAST)
#   6. Outputs curated candidate list → ready for step 3 plots
#
# How to run (CYP — the default):
#   sbatch scripts/08_gene_families/run_filter_genelist.sbatch DC
#   sbatch scripts/08_gene_families/run_filter_genelist.sbatch DC CYP
#
# How to run (OMT — O-methyltransferases):
#   sbatch scripts/08_gene_families/run_filter_genelist.sbatch DC OMT
#   sbatch scripts/08_gene_families/run_filter_genelist.sbatch SK OMT
#   sbatch scripts/08_gene_families/run_filter_genelist.sbatch DG OMT
#
# Custom gene list (overrides family default):
#   GENE_LIST=/path/to/gene_list.txt \
#     sbatch scripts/08_gene_families/run_filter_genelist.sbatch SK
#
# Custom cutoffs:
#   PADJ=0.01 LFC=3.0 sbatch scripts/08_gene_families/run_filter_genelist.sbatch DC OMT
#
# Custom output name:
#   OUT_NAME=omt_candidates.tsv \
#     sbatch scripts/08_gene_families/run_filter_genelist.sbatch DC OMT
#
# Match previous student's min-count threshold (20 instead of default 10):
#   MIN_COUNTS=20 sbatch scripts/08_gene_families/run_filter_genelist.sbatch DC
#
# Output:
#   CYP: 07_NRdatabase/sukman_database/geneious_candidates_DC.tsv  (default)
#   OMT: 07_NRdatabase/sukman_database/omt_candidates_DC.tsv
#   → ready for: sbatch scripts/05_pydeseq2/run_step3_plots.sbatch DC <output>
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
source "${SLURM_SUBMIT_DIR:-.}/scripts/config.sh"

source ~/.bashrc
conda activate rnaseq


# ============================================================================
# GET SPECIES + GENE FAMILY FROM COMMAND LINE
# ============================================================================

# CHANGED: 2026-04-01 — Added 2nd positional arg (CYP or OMT) so you don't
# need to set GENE_FAMILY + GENE_LIST env vars manually for OMT.
SPECIES="${1:-}"
FAMILY_ARG="${2:-}"

if [ -z "${SPECIES}" ]; then
    echo "ERROR: You need to tell me which species!"
    echo ""
    echo "Usage:"
    echo "  sbatch scripts/08_gene_families/run_filter_genelist.sbatch SPECIES [FAMILY]"
    echo ""
    echo "  SPECIES: DC, SK, DG, MF"
    echo "  FAMILY:  CYP (default) or OMT"
    echo ""
    echo "Examples (CYP — default):"
    echo "  sbatch scripts/08_gene_families/run_filter_genelist.sbatch DC CYP"
    echo ""
    echo "Examples (OMT):"
    echo "  sbatch scripts/08_gene_families/run_filter_genelist.sbatch DC OMT"
    echo "  sbatch scripts/08_gene_families/run_filter_genelist.sbatch SK OMT"
    echo ""
    echo "Custom gene list (overrides family default):"
    echo "  GENE_LIST=/path/to/genes.txt sbatch scripts/08_gene_families/run_filter_genelist.sbatch SK"
    exit 1
fi

SPECIES=$(echo "${SPECIES}" | tr '[:lower:]' '[:upper:]')
get_sample_info "${SPECIES}" || exit 1

# Normalize family argument to uppercase (CYP or OMT)
if [ -n "${FAMILY_ARG}" ]; then
    FAMILY_ARG=$(echo "${FAMILY_ARG}" | tr '[:lower:]' '[:upper:]')
fi


# ============================================================================
# PATHS — uses species + gene family from command line
# ============================================================================

COUNT_MATRIX="${BASE_DIR}/03_count_tables/${SPECIES_DIR}/gene_count_matrix.tsv"
METADATA="${BASE_DIR}/03_count_tables/${SPECIES_DIR}/sample_metadata.tsv"

# CHANGED: 2026-04-01 — Gene family selection via 2nd positional arg.
# The 2nd positional arg (FAMILY_ARG) sets defaults for GENE_LIST, GENE_FAMILY,
# and OUT_NAME. Environment variables still override if set explicitly.
#
#   CYP (default): P450_list_RefSeq.txt → geneious_candidates_DC.tsv
#   OMT:           Methyltransferase_list.txt → omt_candidates_DC.tsv
if [ "${FAMILY_ARG}" = "OMT" ]; then
    GENE_LIST="${GENE_LIST:-${BASE_DIR}/07_NRdatabase/sukman_database/Methyltransferase_list.txt}"
    GENE_FAMILY="${GENE_FAMILY:-OMT}"
    OUT_NAME="${OUT_NAME:-omt_candidates_${SPECIES}.tsv}"
else
    GENE_LIST="${GENE_LIST:-${BASE_DIR}/07_NRdatabase/sukman_database/P450_list_RefSeq.txt}"
    GENE_FAMILY="${GENE_FAMILY:-CYP}"
    OUT_NAME="${OUT_NAME:-geneious_candidates_${SPECIES}.tsv}"
fi

OUT_DIR="${BASE_DIR}/07_NRdatabase/sukman_database"
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
echo "  Species:       ${SPECIES_NAME} (${SPECIES})"
echo "  Count matrix:  ${COUNT_MATRIX}"
echo "  Metadata:      ${METADATA}"
echo "  Gene list:     ${GENE_LIST}"
echo "  Gene family:   ${GENE_FAMILY}  (used for heatmap gene detection)"
echo "  DE cutoffs:    padj ≤ ${PADJ}, |log2FC| ≥ ${LFC}"
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
    --gene-family "$GENE_FAMILY" \
    --output "$OUTPUT"


# ============================================================================
# SUMMARY
# ============================================================================

echo ""
echo "========================================================================"
echo "  DONE! $(date)"
echo ""
echo "  Species: ${SPECIES_NAME} (${SPECIES})"
echo "  Output:  ${OUTPUT}"
if [ -f "$OUTPUT" ]; then
    N_GENES=$(tail -n +2 "$OUTPUT" | wc -l)
    echo "  DE candidates: ${N_GENES} genes"
fi
echo ""
echo "  NEXT — generate plots (including ${GENE_FAMILY} heatmap):"
echo "    sbatch scripts/05_pydeseq2/run_step3_plots.sbatch ${SPECIES} ${OUTPUT}"
echo "========================================================================"
