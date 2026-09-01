#!/bin/bash
################################################################################
# Build the REPLACEMENT candidate TSV for heatmaps (longest isoform / locus)
################################################################################
#
# What this does:
#   Reads the full DE candidate table (e.g. cyp_candidates_MF.tsv, ~304 genes).
#   For each Trinity component (TRINITY_DN###_c#), keeps ONE gene: the gene
#   whose longest isoform is the longest.
#   Writes a NEW replacement TSV. The original 304-gene file is NOT changed.
#
# Why:
#   Trinity often splits one locus into several gene IDs (_g1, _g3, _c1_g1).
#   Plotting all of them makes the heatmap look overcrowded. This replacement
#   file is what you pass to run_step3_plots.sbatch instead.
#
# How to run (from 05_rnaseq-code):
#   sbatch scripts/08_gene_families/run_collapse_longest_isoform.sbatch CYP
#   sbatch scripts/08_gene_families/run_collapse_longest_isoform.sbatch OMT
#   sbatch scripts/08_gene_families/run_collapse_longest_isoform.sbatch both
#
# Output (same folder as the original candidates):
#   07_NRdatabase/00_2_MF/cyp_candidates_MF_longest.tsv          ← plot THIS
#   07_NRdatabase/00_2_MF/cyp_candidates_MF_longest_dropped.tsv  (removed genes)
#   07_NRdatabase/00_2_MF/cyp_candidates_MF.tsv                  (unchanged)
#
# Next:
#   sbatch scripts/05_pydeseq2/run_step3_plots.sbatch MF \
#     /projects/tholl_lab_1/daisy_analysis/07_NRdatabase/00_2_MF/cyp_candidates_MF_longest.tsv
#
################################################################################

#SBATCH --job-name=collapse_iso
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu

#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis
#SBATCH -o 07_NRdatabase/00_2_MF/collapse_iso_%j.out
#SBATCH -e 07_NRdatabase/00_2_MF/collapse_iso_%j.err

set -euo pipefail

# ===== Load shared config =====
source "${SLURM_SUBMIT_DIR:-.}/scripts/config.sh"

source ~/.bashrc
conda activate rnaseq


# ============================================================================
# ARGUMENTS
# ============================================================================

# FAMILY: CYP (default), OMT, or both
FAMILY_ARG="${1:-CYP}"
FAMILY_ARG=$(echo "${FAMILY_ARG}" | tr '[:lower:]' '[:upper:]')

# Grouping: component = TRINITY_DN###_c# (default). Override with GROUP=dn
GROUP="${GROUP:-component}"

ASSEMBLY="${PROCESSED_DIR}/00_6_cdhit/MF_trinity_cdhit95.fasta"
CAND_DIR="${MF_CANDIDATE_DIR}"
SCRIPT="${CODE_DIR}/scripts/08_gene_families/collapse_longest_isoform.py"


# ============================================================================
# CHECKS
# ============================================================================

if [ "${FAMILY_ARG}" != "CYP" ] && [ "${FAMILY_ARG}" != "OMT" ] && [ "${FAMILY_ARG}" != "BOTH" ]; then
    echo "ERROR: FAMILY must be CYP, OMT, or both (got: ${FAMILY_ARG})"
    echo "  sbatch scripts/08_gene_families/run_collapse_longest_isoform.sbatch CYP"
    exit 1
fi

if [ ! -f "${SCRIPT}" ]; then
    echo "ERROR: Python script not found: ${SCRIPT}"
    echo "  git pull in 05_rnaseq-code so collapse_longest_isoform.py is on the cluster."
    exit 1
fi

if [ ! -f "${ASSEMBLY}" ]; then
    echo "ERROR: Assembly FASTA not found: ${ASSEMBLY}"
    exit 1
fi

mkdir -p "${CAND_DIR}"

echo "========================================================================"
echo "  Build replacement candidate TSV (longest isoform per component)"
echo "  Started: $(date)"
echo "========================================================================"
echo ""
echo "  Family:     ${FAMILY_ARG}"
echo "  Group:      ${GROUP}"
echo "  FASTA:      ${ASSEMBLY}"
echo "  Candidate:  ${CAND_DIR}"
echo ""


# ============================================================================
# RUN
# ============================================================================

# --family both is lowercase "both" in the Python script
if [ "${FAMILY_ARG}" = "BOTH" ]; then
    PY_FAMILY="both"
else
    PY_FAMILY="${FAMILY_ARG}"
fi

python "${SCRIPT}" \
    --family "${PY_FAMILY}" \
    --fasta "${ASSEMBLY}" \
    --group "${GROUP}" \
    --out-dir "${CAND_DIR}"


# ============================================================================
# SUMMARY
# ============================================================================

echo ""
echo "========================================================================"
echo "  DONE! $(date)"
echo ""
if [ "${FAMILY_ARG}" = "CYP" ] || [ "${FAMILY_ARG}" = "BOTH" ]; then
    echo "  Original:    ${CAND_DIR}/cyp_candidates_MF.tsv"
    echo "  Replacement: ${CAND_DIR}/cyp_candidates_MF_longest.tsv"
    if [ -f "${CAND_DIR}/cyp_candidates_MF_longest.tsv" ]; then
        echo "  Replacement rows: $(tail -n +2 "${CAND_DIR}/cyp_candidates_MF_longest.tsv" | wc -l)"
    fi
    echo ""
    echo "  NEXT — CYP heatmap from the replacement file:"
    echo "    sbatch scripts/05_pydeseq2/run_step3_plots.sbatch MF \\"
    echo "      ${CAND_DIR}/cyp_candidates_MF_longest.tsv"
    echo ""
fi
if [ "${FAMILY_ARG}" = "OMT" ] || [ "${FAMILY_ARG}" = "BOTH" ]; then
    echo "  Original:    ${CAND_DIR}/omt_candidates_MF.tsv"
    echo "  Replacement: ${CAND_DIR}/omt_candidates_MF_longest.tsv"
    if [ -f "${CAND_DIR}/omt_candidates_MF_longest.tsv" ]; then
        echo "  Replacement rows: $(tail -n +2 "${CAND_DIR}/omt_candidates_MF_longest.tsv" | wc -l)"
    fi
    echo ""
    echo "  NEXT — OMT heatmap from the replacement file:"
    echo "    sbatch scripts/05_pydeseq2/run_step3_plots.sbatch MF \\"
    echo "      ${CAND_DIR}/omt_candidates_MF_longest.tsv"
    echo ""
fi
echo "========================================================================"
