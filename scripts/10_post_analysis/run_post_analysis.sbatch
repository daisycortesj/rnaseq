#!/bin/bash
################################################################################
# POST-ANALYSIS — Phylo Trees + Genomic Clustering
################################################################################
#
# PURPOSE:
#   Run phylogenetic tree + heatmap and genomic clustering AFTER the CYP
#   heatmap pipeline paths (A, B, or C) are complete. Uses the expressed
#   gene lists and protein FASTAs produced by earlier pipeline steps.
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
#   sbatch scripts/post_pydeq2.sbatch DC cyp
#
#   # Geneious P450 list:
#   sbatch scripts/post_pydeq2.sbatch DC geneious
#
#   # OMT / methyltransferases:
#   sbatch scripts/post_pydeq2.sbatch DC omt
#
#   # Skip phylo tree (only genomic clustering):
#   SKIP_PHYLO=1 sbatch scripts/post_pydeq2.sbatch DC cyp
#
#   # Skip genomic clustering (only phylo tree):
#   SKIP_GENOMIC=1 sbatch scripts/post_pydeq2.sbatch DC cyp
#
#   # Fetch NCBI gene descriptions (requires internet):
#   FETCH_NCBI=1 sbatch scripts/post_pydeq2.sbatch DC cyp
#
#   # Custom genomic cluster threshold (default 50 kb):
#   THRESHOLD=100000 sbatch scripts/post_pydeq2.sbatch DC cyp
#
#   # Only significant / upregulated / downregulated:
#   SUBSET=significant sbatch scripts/post_pydeq2.sbatch DC cyp
#   SUBSET=upregulated sbatch scripts/post_pydeq2.sbatch DC cyp
#
#   # Chain after extract step:
#   JOB1=$(sbatch --parsable scripts/run_cyp_express_extract.sbatch)
#   sbatch --dependency=afterok:$JOB1 scripts/post_pydeq2.sbatch DC cyp
#
# INPUT FILES (where they come from):
#
#   File                          Created by                                HPC path
#   ───────────────────────────── ───────────────────────────────────────── ─────────────────────────────────────
#   gene_count_matrix.tsv         featurecounts.sbatch DC                  03_count_tables/00_1_DC/
#   sample_metadata.tsv           featurecounts.sbatch DC                  03_count_tables/00_1_DC/
#   dc_genomic.gtf                NCBI reference download                  04_reference/
#   cyp_expressed_list.tsv        run_cyp_express_extract.sbatch           07_NRdatabase/cyp450_database/
#   cyp_proteins.fasta            run_cyp_express_extract.sbatch           07_NRdatabase/cyp450_database/
#   geneious_expressed_list.tsv   PREV_LIST=... run_cyp_express_extract    07_NRdatabase/cyp450_database/
#   geneious_proteins.fasta       PREV_LIST=... run_cyp_express_extract    07_NRdatabase/cyp450_database/
#   omt_expressed_list.tsv        (same process for OMT family)            07_NRdatabase/omt_database/
#   omt_proteins.fasta            (same process for OMT family)            07_NRdatabase/omt_database/
#
# OUTPUT:
#   06_analysis/post_analysis_DC_CYP/
#     phylo/
#       CYP_phylo_heatmap.pdf/png       — phylo tree + expression heatmap
#       CYP_clustered_heatmap.pdf/png   — expression-clustered heatmap
#       msa_alignment.fasta             — MAFFT alignment
#       phylo_tree.nwk                  — Newick tree (reusable)
#     genomic/
#       CYP_genomic_coordinates.tsv     — gene positions on chromosomes
#       CYP_clusters.tsv                — cluster assignments
#       CYP_chromosome_map.pdf/png      — full chromosome overview
#       CYP_distance_*.pdf              — distance-from-anchor plots
#
################################################################################

#SBATCH --job-name=post_pydeq2
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu

#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis
#SBATCH -o 06_analysis/post_pydeq2_%j.out
#SBATCH -e 06_analysis/post_pydeq2_%j.err

set -euo pipefail

BASE_DIR="/projects/tholl_lab_1/daisy_analysis"
CODE_DIR="${BASE_DIR}/05_rnaseq-code"

SPECIES="${1:-DC}"
FAMILY_KEY="${2:-cyp}"
SPECIES=$(echo "${SPECIES}" | tr '[:lower:]' '[:upper:]')
FAMILY_KEY=$(echo "${FAMILY_KEY}" | tr '[:upper:]' '[:lower:]')

# ── Settings ─────────────────────────────────────────────────────────────────
THRESHOLD="${THRESHOLD:-50000}"
SKIP_PHYLO="${SKIP_PHYLO:-0}"
SKIP_GENOMIC="${SKIP_GENOMIC:-0}"
FETCH_NCBI="${FETCH_NCBI:-0}"
SUBSET="${SUBSET:-}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"
EMAIL="${EMAIL:-daisycortesj@vt.edu}"

# ── Resolve family-specific paths ────────────────────────────────────────────
case "${FAMILY_KEY}" in
    cyp)
        EXPRESSED="${BASE_DIR}/07_NRdatabase/cyp450_database/cyp_expressed_list.tsv"
        FASTA="${BASE_DIR}/07_NRdatabase/cyp450_database/cyp_proteins.fasta"
        FAMILY_LABEL="CYP"
        ;;
    geneious)
        EXPRESSED="${BASE_DIR}/07_NRdatabase/cyp450_database/geneious_expressed_list.tsv"
        FASTA="${BASE_DIR}/07_NRdatabase/cyp450_database/geneious_proteins.fasta"
        FAMILY_LABEL="CYP_geneious"
        ;;
    omt|mt)
        EXPRESSED="${BASE_DIR}/07_NRdatabase/omt_database/omt_expressed_list.tsv"
        FASTA="${BASE_DIR}/07_NRdatabase/omt_database/omt_proteins.fasta"
        FAMILY_LABEL="OMT"
        ;;
    *)
        echo "ERROR: Unknown family key: ${FAMILY_KEY}"
        echo "Valid: cyp, geneious, omt"
        exit 1
        ;;
esac

# ── Resolve species-specific paths ───────────────────────────────────────────
case "${SPECIES}" in
    DC) COUNT_FOLDER="00_1_DC" ;;
    DG) COUNT_FOLDER="00_2_DG" ;;
    MF) COUNT_FOLDER="00_3_MF" ;;
    *)  echo "ERROR: Unknown species: ${SPECIES}"; exit 1 ;;
esac

COUNT_MATRIX="${BASE_DIR}/03_count_tables/${COUNT_FOLDER}/gene_count_matrix.tsv"
METADATA="${BASE_DIR}/03_count_tables/${COUNT_FOLDER}/sample_metadata.tsv"
GTF="${BASE_DIR}/04_reference/dc_genomic.gtf"

SUFFIX="${FAMILY_LABEL}"
[ -n "${SUBSET}" ] && SUFFIX="${SUFFIX}_${SUBSET}"
OUTPUT_BASE="${BASE_DIR}/06_analysis/post_analysis_${SPECIES}_${SUFFIX}"
PHYLO_DIR="${OUTPUT_BASE}/phylo"
GENOMIC_DIR="${OUTPUT_BASE}/genomic"

# ── Activate environment ─────────────────────────────────────────────────────
source ~/.bashrc
conda activate rnaseq
module load MAFFT/7.505-GCC-11.3.0-with-extensions 2>/dev/null || true
module load FastTree/2.1.11-GCCcore-11.3.0 2>/dev/null || true

echo "============================================================"
echo "Post-Analysis — ${FAMILY_LABEL} (${SPECIES})"
echo "============================================================"
echo "  Expressed list: ${EXPRESSED}"
echo "  Protein FASTA:  ${FASTA}"
echo "  Count matrix:   ${COUNT_MATRIX}"
echo "  GTF:            ${GTF}"
echo "  Output:         ${OUTPUT_BASE}"
[ -n "${SUBSET}" ] && echo "  Subset:         ${SUBSET}"
echo ""

# ── Validate inputs ──────────────────────────────────────────────────────────
MISSING=0
for f in "${COUNT_MATRIX}" "${METADATA}"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Required file not found: $f"
        echo "  → Run: sbatch scripts/featurecounts.sbatch ${SPECIES}"
        MISSING=1
    fi
done

if [ ! -f "${EXPRESSED}" ]; then
    echo "ERROR: Expressed gene list not found: ${EXPRESSED}"
    echo "  → Run: sbatch scripts/run_cyp_express_extract.sbatch"
    MISSING=1
fi

if [ ! -f "${GTF}" ]; then
    echo "ERROR: GTF not found: ${GTF}"
    MISSING=1
fi

[ "${MISSING}" = "1" ] && exit 1

# ══════════════════════════════════════════════════════════════════════════════
# STEP 1: PHYLOGENETIC TREE + HEATMAP
# ══════════════════════════════════════════════════════════════════════════════
if [ "${SKIP_PHYLO}" != "1" ]; then
    echo "────────────────────────────────────────────────────────────"
    echo "PHYLO TREE + HEATMAP"
    echo "────────────────────────────────────────────────────────────"
    mkdir -p "${PHYLO_DIR}"

    PHYLO_CMD=(python3 "${CODE_DIR}/scripts/10_post_analysis/phylo_heatmap.py"
               --counts "${COUNT_MATRIX}"
               --metadata "${METADATA}"
               --expression "${EXPRESSED}"
               --output-dir "${PHYLO_DIR}"
               --family "${FAMILY_LABEL}"
               --threads "${THREADS}")

    if [ -f "${FASTA}" ]; then
        PHYLO_CMD+=(--fasta "${FASTA}")
    else
        echo "  WARNING: Protein FASTA not found (${FASTA})"
        echo "  Running in cluster-only mode (no phylo tree)"
        PHYLO_CMD+=(--cluster-only)
    fi

    [ -n "${SUBSET}" ] && PHYLO_CMD+=(--subset "${SUBSET}")

    echo "  Running: ${PHYLO_CMD[*]}"
    echo ""
    "${PHYLO_CMD[@]}" || echo "  WARNING: Phylo step had errors (continuing)"
    echo ""
else
    echo "  Skipping phylo tree (SKIP_PHYLO=1)"
fi

# ══════════════════════════════════════════════════════════════════════════════
# STEP 2: GENOMIC CLUSTERING
# ══════════════════════════════════════════════════════════════════════════════
if [ "${SKIP_GENOMIC}" != "1" ]; then
    echo "────────────────────────────────────────────────────────────"
    echo "GENOMIC CLUSTERING"
    echo "────────────────────────────────────────────────────────────"
    mkdir -p "${GENOMIC_DIR}"

    GENOMIC_CMD=(python3 "${CODE_DIR}/scripts/10_post_analysis/genomic_cluster_analysis.py"
                 --gene-list "${EXPRESSED}"
                 --gtf "${GTF}"
                 --output-dir "${GENOMIC_DIR}"
                 --family "${FAMILY_LABEL}"
                 --threshold "${THRESHOLD}")

    if [ "${FETCH_NCBI}" = "1" ]; then
        GENOMIC_CMD+=(--fetch-ncbi --email "${EMAIL}")
    fi

    echo "  Running: ${GENOMIC_CMD[*]}"
    echo ""
    "${GENOMIC_CMD[@]}" || echo "  WARNING: Genomic clustering had errors (continuing)"
    echo ""
else
    echo "  Skipping genomic clustering (SKIP_GENOMIC=1)"
fi

# ══════════════════════════════════════════════════════════════════════════════
# DONE
# ══════════════════════════════════════════════════════════════════════════════
echo ""
echo "============================================================"
echo "POST-ANALYSIS COMPLETE — ${FAMILY_LABEL} (${SPECIES})"
echo "============================================================"
echo ""
echo "Output:"
[ -d "${PHYLO_DIR}" ] && echo "  Phylo:    ${PHYLO_DIR}/" && ls "${PHYLO_DIR}/" 2>/dev/null
echo ""
[ -d "${GENOMIC_DIR}" ] && echo "  Genomic:  ${GENOMIC_DIR}/" && ls "${GENOMIC_DIR}/" 2>/dev/null
echo ""
echo "Download:"
echo "  scp -r daisycortesj@tinkercliffs.arc.vt.edu:${OUTPUT_BASE}/ ~/Desktop/post_analysis/"
echo ""
echo "End time: $(date)"
