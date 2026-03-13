#!/bin/bash
################################################################################
# GENOMIC CLUSTER ANALYSIS
################################################################################
#
# PURPOSE:
#   Find gene clusters on chromosomes for a gene family (CYP, OMT, etc.).
#   Parses gene coordinates from the GTF, calculates inter-gene distances,
#   identifies clusters within a threshold, and generates plots.
#
#   This is the Python equivalent of the previous student's R workflow
#   that used rentrez + manual distance calculation.
#
# PREREQUISITES (run these first):
#   1. featurecounts.sbatch DC
#      → gene_count_matrix.tsv + sample_metadata.tsv
#
#   2. run_pydeseq2_step1_analysis.sbatch DC
#      → pydeseq2_results_UNFILTERED.tsv
#
#   3. run_cyp_express_extract.sbatch
#      → cyp_expressed_list.tsv + cyp_proteins.fasta
#      (add PREV_LIST=...P450_list_RefSeq.txt for geneious_expressed_list.tsv)
#
# HOW TO RUN:
#   # CYP genes (default):
#   sbatch scripts/run_genomic_clustering.sbatch DC cyp
#
#   # Geneious P450 list:
#   sbatch scripts/run_genomic_clustering.sbatch DC geneious
#
#   # OMT / methyltransferases:
#   sbatch scripts/run_genomic_clustering.sbatch DC omt
#
#   # Custom threshold (default 50 kb):
#   THRESHOLD=100000 sbatch scripts/run_genomic_clustering.sbatch DC cyp
#
#   # Fetch NCBI gene descriptions (requires internet):
#   FETCH_NCBI=1 sbatch scripts/run_genomic_clustering.sbatch DC cyp
#
#   # Chain after extract step:
#   JOB1=$(sbatch --parsable scripts/run_cyp_express_extract.sbatch)
#   sbatch --dependency=afterok:$JOB1 scripts/run_genomic_clustering.sbatch DC cyp
#
# INPUT FILES (where they come from):
#
#   File                          Created by                                HPC path
#   ───────────────────────────── ───────────────────────────────────────── ─────────────────────────────────────
#   dc_genomic.gtf                NCBI reference download                  04_reference/
#   pydeseq2_results_UNFILTERED   run_pydeseq2_step1_analysis.sbatch DC    06_analysis/pydeseq2_DC_step1_unfiltered/
#   cyp_expressed_list.tsv        run_cyp_express_extract.sbatch           07_NRdatabase/cyp450_database/
#   geneious_expressed_list.tsv   PREV_LIST=... run_cyp_express_extract    07_NRdatabase/cyp450_database/
#   omt_expressed_list.tsv        (same process for OMT family)            07_NRdatabase/omt_database/
#
# OUTPUT:
#   06_analysis/genomic_clusters_DC_CYP/
#     CYP_genomic_coordinates.tsv     — gene positions on chromosomes (+ expression stats)
#     CYP_clusters.tsv                — cluster assignments (cluster_id, size, position)
#     CYP_distance_matrix.tsv         — pairwise distances between genes
#     CYP_chromosome_map.pdf/png      — overview: all genes on chromosomes, clusters in red
#     CYP_distance_*.pdf              — distance-from-anchor plots per chromosome
#     CYP_cluster_*_detail.pdf        — zoomed view of each cluster
#
################################################################################

#SBATCH --job-name=genomic_cluster
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=30:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu

#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis
#SBATCH -o 06_analysis/genomic_cluster_%j.out
#SBATCH -e 06_analysis/genomic_cluster_%j.err

set -euo pipefail

BASE_DIR="/projects/tholl_lab_1/daisy_analysis"
CODE_DIR="${BASE_DIR}/05_rnaseq-code"

# ── Arguments ────────────────────────────────────────────────────────────────
SPECIES="${1:-DC}"
FAMILY_KEY="${2:-cyp}"
SPECIES=$(echo "${SPECIES}" | tr '[:lower:]' '[:upper:]')
FAMILY_KEY=$(echo "${FAMILY_KEY}" | tr '[:upper:]' '[:lower:]')

# ── Settings (all overridable via env vars) ──────────────────────────────────
THRESHOLD="${THRESHOLD:-50000}"
FETCH_NCBI="${FETCH_NCBI:-0}"
EMAIL="${EMAIL:-daisycortesj@vt.edu}"

# ── Resolve family-specific paths ────────────────────────────────────────────
case "${FAMILY_KEY}" in
    cyp)
        GENE_LIST="${BASE_DIR}/07_NRdatabase/cyp450_database/cyp_expressed_list.tsv"
        FAMILY_LABEL="CYP"
        ;;
    geneious)
        GENE_LIST="${BASE_DIR}/07_NRdatabase/cyp450_database/geneious_expressed_list.tsv"
        FAMILY_LABEL="CYP_geneious"
        ;;
    omt|mt)
        GENE_LIST="${BASE_DIR}/07_NRdatabase/omt_database/omt_expressed_list.tsv"
        FAMILY_LABEL="OMT"
        ;;
    *)
        GENE_LIST="${FAMILY_KEY}"
        FAMILY_LABEL="custom"
        ;;
esac

GENE_LIST="${GENE_LIST_OVERRIDE:-${GENE_LIST}}"

# ── Resolve species-specific paths ───────────────────────────────────────────
case "${SPECIES}" in
    DC|DG)
        GTF="${BASE_DIR}/04_reference/dc_genomic.gtf"
        ;;
    *)
        echo "ERROR: No GTF for species ${SPECIES}"
        exit 1
        ;;
esac

DESEQ_RESULTS="${BASE_DIR}/06_analysis/pydeseq2_${SPECIES}_step1_unfiltered/pydeseq2_results_UNFILTERED.tsv"
OUTPUT_DIR="${BASE_DIR}/06_analysis/genomic_clusters_${SPECIES}_${FAMILY_LABEL}"

# ── Activate environment ─────────────────────────────────────────────────────
source ~/.bashrc
conda activate rnaseq

echo "============================================================"
echo "Genomic Cluster Analysis — ${FAMILY_LABEL} (${SPECIES})"
echo "============================================================"
echo "  Gene list:   ${GENE_LIST}"
echo "  GTF:         ${GTF}"
echo "  Threshold:   ${THRESHOLD} bp ($(( THRESHOLD / 1000 )) kb)"
echo "  Output:      ${OUTPUT_DIR}"
[ "${FETCH_NCBI}" = "1" ] && echo "  NCBI fetch:  enabled"
echo ""

# ── Validate inputs ──────────────────────────────────────────────────────────
MISSING=0
if [ ! -f "${GENE_LIST}" ]; then
    echo "ERROR: Gene list not found: ${GENE_LIST}"
    echo "  → Run: sbatch scripts/run_cyp_express_extract.sbatch"
    MISSING=1
fi

if [ ! -f "${GTF}" ]; then
    echo "ERROR: GTF not found: ${GTF}"
    MISSING=1
fi

[ "${MISSING}" = "1" ] && exit 1

mkdir -p "${OUTPUT_DIR}"

# ── Build command ────────────────────────────────────────────────────────────
CMD=(python3 "${CODE_DIR}/scripts/10_post_analysis/genomic_cluster_analysis.py"
     --gene-list "${GENE_LIST}"
     --gtf "${GTF}"
     --output-dir "${OUTPUT_DIR}"
     --threshold "${THRESHOLD}"
     --family "${FAMILY_LABEL}")

if [ -f "${DESEQ_RESULTS}" ]; then
    CMD+=(--deseq "${DESEQ_RESULTS}")
fi

if [ "${FETCH_NCBI}" = "1" ]; then
    CMD+=(--fetch-ncbi --email "${EMAIL}")
fi

echo "  Running: ${CMD[*]}"
echo ""
"${CMD[@]}"

echo ""
echo "============================================================"
echo "COMPLETE — Genomic Cluster Analysis"
echo "============================================================"
echo "Output: ${OUTPUT_DIR}/"
ls -la "${OUTPUT_DIR}/"
echo ""
echo "Download:"
echo "  scp -r daisycortesj@tinkercliffs.arc.vt.edu:${OUTPUT_DIR}/ ~/Desktop/genomic_clusters/"
echo ""
echo "End time: $(date)"
