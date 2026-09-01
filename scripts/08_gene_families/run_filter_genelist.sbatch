#!/bin/bash
################################################################################
# Filter count matrix by gene list + run PyDESeq2 for expression stats
################################################################################
#
# What this does:
#   1. Loads the count matrix (carrot ~33K genes; MF Trinity ~330K genes)
#   2. Runs PyDESeq2 on ALL genes (uses full dataset for normalization)
#      MF Trinity LFC fitting can take >1 hour — do not use a 1h walltime.
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
# How to run (MF Trinity — nutmeg fruit vs leaf):
#   Needs 03_count_tables/00_5_MF_trinity/gene_count_matrix.tsv first
#   (sbatch scripts/04_counting/run_count_matrix.sbatch MF --type rsem)
#   Default contrast for MF: fruit (F) vs leaf (L)
#     F = test group (numerator)  → positive log2FC = higher in fruit
#     L = baseline  (denominator)
#   Default gene list: 06_analysis/00_3_MF/combined_cyp450_genes.txt
#   sbatch scripts/08_gene_families/run_filter_genelist.sbatch MF CYP
#   sbatch scripts/08_gene_families/run_filter_genelist.sbatch MF OMT
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
#   MF CYP: 07_NRdatabase/00_2_MF/cyp_candidates_MF.tsv           (full DE list)
#           07_NRdatabase/00_2_MF/cyp_candidates_MF_longest.tsv  (heatmap replacement)
#   MF OMT: 07_NRdatabase/00_2_MF/omt_candidates_MF.tsv
#           07_NRdatabase/00_2_MF/omt_candidates_MF_longest.tsv
#   DC CYP: 07_NRdatabase/sukman_database/geneious_candidates_DC.tsv  (default)
#   DC OMT: 07_NRdatabase/sukman_database/omt_candidates_DC.tsv
#   → plot MF from the *_longest.tsv replacement:
#     sbatch scripts/05_pydeseq2/run_step3_plots.sbatch MF \
#       /projects/tholl_lab_1/daisy_analysis/07_NRdatabase/00_2_MF/cyp_candidates_MF_longest.tsv
#
################################################################################

#SBATCH --job-name=filter_genelist
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
# Job 7298602 died at 1h during LFC fit on ~126k Trinity genes.
# This should finish in ~2–4h. 8h is a buffer.
# Do NOT copy BLASTX's 6-day limit: long requests wait longer in the queue.
#SBATCH --time=8:00:00
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

# Pick the count-table folder (MF Trinity → 00_5_MF_trinity, not 00_3_MF).
# SPECIES_DIR is the raw-data folder. COUNT_FOLDER is where PyDESeq2 counts live.
# Override: COUNT_FOLDER=00_3_MF sbatch ... MF
resolve_count_folder

# Pick contrast: A = test group (numerator), B = baseline (denominator).
# MF default: F vs L so positive log2FC = higher in fruit.
# Override: CONTRAST_A=L CONTRAST_B=F sbatch ... MF
resolve_contrast

# Normalize family argument to uppercase (CYP or OMT)
if [ -n "${FAMILY_ARG}" ]; then
    FAMILY_ARG=$(echo "${FAMILY_ARG}" | tr '[:lower:]' '[:upper:]')
fi


# ============================================================================
# PATHS — uses species + gene family from command line
# ============================================================================

COUNT_MATRIX="${BASE_DIR}/03_count_tables/${COUNT_FOLDER}/gene_count_matrix.tsv"
METADATA="${BASE_DIR}/03_count_tables/${COUNT_FOLDER}/sample_metadata.tsv"

# CHANGED: 2026-04-01 — Gene family selection via 2nd positional arg.
# The 2nd positional arg (FAMILY_ARG) sets defaults for GENE_LIST, GENE_FAMILY,
# and OUT_NAME. Environment variables still override if set explicitly.
#
#   CYP (default): P450_list_RefSeq.txt → geneious_candidates_DC.tsv
#   OMT:           Methyltransferase_list.txt → omt_candidates_DC.tsv
#   MF CYP:        06_analysis/<SPECIES_DIR>/combined_cyp450_genes.txt
#   MF OMT:        06_analysis/<SPECIES_DIR>/combined_omt_genes.txt
if [ "${FAMILY_ARG}" = "OMT" ]; then
    GENE_FAMILY="${GENE_FAMILY:-OMT}"
    OUT_NAME="${OUT_NAME:-omt_candidates_${SPECIES}.tsv}"
    if [ "${SPECIES}" = "MF" ]; then
        GENE_LIST="${GENE_LIST:-${ANALYSIS_DIR}/${SPECIES_DIR}/combined_omt_genes.txt}"
    else
        GENE_LIST="${GENE_LIST:-${BASE_DIR}/07_NRdatabase/sukman_database/Methyltransferase_list.txt}"
    fi
else
    GENE_FAMILY="${GENE_FAMILY:-CYP}"
    if [ "${SPECIES}" = "MF" ]; then
        GENE_LIST="${GENE_LIST:-${ANALYSIS_DIR}/${SPECIES_DIR}/combined_cyp450_genes.txt}"
        OUT_NAME="${OUT_NAME:-cyp_candidates_${SPECIES}.tsv}"
    else
        GENE_LIST="${GENE_LIST:-${BASE_DIR}/07_NRdatabase/sukman_database/P450_list_RefSeq.txt}"
        OUT_NAME="${OUT_NAME:-geneious_candidates_${SPECIES}.tsv}"
    fi
fi

# MF candidate tables live with the nutmeg database folder, not Sukman's carrot lists.
# Override: OUT_DIR=/path sbatch ... MF CYP
if [ "${SPECIES}" = "MF" ]; then
    OUT_DIR="${OUT_DIR:-${MF_CANDIDATE_DIR}}"
else
    OUT_DIR="${OUT_DIR:-${DB_DIR}/sukman_database}"
fi
OUTPUT="${OUT_DIR}/${OUT_NAME}"

PADJ="${PADJ:-0.05}"
LFC="${LFC:-2.0}"
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
echo "                 ${CONTRAST_A} = test group (numerator)"
echo "                 ${CONTRAST_B} = baseline (denominator)"
echo "                 positive log2FC = higher in ${CONTRAST_A}"
echo "                 negative log2FC = higher in ${CONTRAST_B}"
echo "  Output:        ${OUTPUT}"
echo ""

# The gene list (CYP/OMT IDs) is NOT enough. PyDESeq2 also needs the
# full count matrix + metadata so it can normalize using ALL genes.
if [ ! -f "$COUNT_MATRIX" ]; then
    echo "ERROR: File not found: $COUNT_MATRIX"
    echo ""
    echo "  That is the COUNT MATRIX (reads per gene per sample), not the gene list."
    echo "  Your gene list can exist and this step still fails without counts."
    echo ""
    if [ -n "${TRINITY_COUNT_DIR:-}" ]; then
        echo "  For ${SPECIES} Trinity, the matrix should be here:"
        echo "    ${BASE_DIR}/03_count_tables/${TRINITY_COUNT_DIR}/gene_count_matrix.tsv"
        echo ""
        echo "  If that file is missing, build it first (after RSEM):"
        echo "    sbatch scripts/04_counting/run_count_matrix.sbatch ${SPECIES} --type rsem"
    else
        echo "  Build it with featureCounts:"
        echo "    sbatch scripts/04_counting/run_featurecounts.sbatch ${SPECIES}"
    fi
    echo ""
    exit 1
fi

for f in "$METADATA" "$GENE_LIST"; do
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
# MF: replacement heatmap TSV (one gene per Trinity component)
# ============================================================================
# The full DE table above can have several Trinity gene IDs per locus
# (_g1, _g3, _c1_g1). For heatmaps we keep the gene whose longest isoform
# is longest. Original ${OUTPUT} is not overwritten.
REPLACEMENT="${OUTPUT%.tsv}_longest.tsv"
if [ "${SPECIES}" = "MF" ] && [ -f "$OUTPUT" ]; then
    ASSEMBLY="${PROCESSED_DIR}/00_6_cdhit/${SPECIES_CODE_UPPER}_trinity_cdhit95.fasta"
    echo ""
    echo "  Building heatmap replacement (longest isoform per component)..."
    python "${CODE_DIR}/scripts/08_gene_families/collapse_longest_isoform.py" \
        --candidates "$OUTPUT" \
        --fasta "$ASSEMBLY" \
        --group component \
        --out-dir "$OUT_DIR" \
        --output "$REPLACEMENT"
fi


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
if [ -f "$REPLACEMENT" ]; then
    N_KEEP=$(tail -n +2 "$REPLACEMENT" | wc -l)
    echo "  Heatmap replacement: ${REPLACEMENT}"
    echo "  Replacement genes:   ${N_KEEP}"
    echo ""
    echo "  NEXT — generate plots from the replacement file:"
    echo "    sbatch scripts/05_pydeseq2/run_step3_plots.sbatch ${SPECIES} ${REPLACEMENT}"
else
    echo ""
    echo "  NEXT — generate plots (including ${GENE_FAMILY} heatmap):"
    echo "    sbatch scripts/05_pydeseq2/run_step3_plots.sbatch ${SPECIES} ${OUTPUT}"
fi
echo "========================================================================"
