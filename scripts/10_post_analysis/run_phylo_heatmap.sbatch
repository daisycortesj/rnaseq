#!/bin/bash
################################################################################
# PHYLOGENETIC TREE + EXPRESSION HEATMAP
################################################################################
#
# PURPOSE:
#   Build a phylogenetic tree from protein sequences and combine it with an
#   expression heatmap — the Python equivalent of ggtree + gheatmap in R.
#
#   Three modes:
#     1. Build tree: FASTA → MAFFT alignment → FastTree → tree + heatmap
#     2. Load tree:  Newick file → tree + heatmap
#     3. Cluster:    No tree, hierarchical clustering of expression data
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
#   # Build tree from CYP proteins + show expression:
#   sbatch scripts/run_phylo_heatmap.sbatch DC cyp
#
#   # Build tree from Geneious proteins:
#   sbatch scripts/run_phylo_heatmap.sbatch DC geneious
#
#   # OMT / methyltransferases:
#   sbatch scripts/run_phylo_heatmap.sbatch DC omt
#
#   # Only significant genes:
#   SUBSET=significant sbatch scripts/run_phylo_heatmap.sbatch DC cyp
#
#   # Only upregulated:
#   SUBSET=upregulated sbatch scripts/run_phylo_heatmap.sbatch DC cyp
#
#   # Use a pre-built Newick tree:
#   TREE=/path/to/tree.nwk sbatch scripts/run_phylo_heatmap.sbatch DC cyp
#
#   # Expression-only clustering (no sequences/tree):
#   CLUSTER_ONLY=1 sbatch scripts/run_phylo_heatmap.sbatch DC cyp
#
#   # Chain after extract step:
#   JOB1=$(sbatch --parsable scripts/run_cyp_express_extract.sbatch)
#   sbatch --dependency=afterok:$JOB1 scripts/run_phylo_heatmap.sbatch DC cyp
#
# INPUT FILES (where they come from):
#
#   File                          Created by                                HPC path
#   ───────────────────────────── ───────────────────────────────────────── ─────────────────────────────────────
#   gene_count_matrix.tsv         featurecounts.sbatch DC                  03_count_tables/00_1_DC/
#   sample_metadata.tsv           featurecounts.sbatch DC                  03_count_tables/00_1_DC/
#   cyp_expressed_list.tsv        run_cyp_express_extract.sbatch           07_NRdatabase/cyp450_database/
#   cyp_proteins.fasta            run_cyp_express_extract.sbatch           07_NRdatabase/cyp450_database/
#   geneious_expressed_list.tsv   PREV_LIST=... run_cyp_express_extract    07_NRdatabase/cyp450_database/
#   geneious_proteins.fasta       PREV_LIST=... run_cyp_express_extract    07_NRdatabase/cyp450_database/
#   omt_expressed_list.tsv        (same process for OMT family)            07_NRdatabase/omt_database/
#   omt_proteins.fasta            (same process for OMT family)            07_NRdatabase/omt_database/
#
# OUTPUT:
#   06_analysis/phylo_heatmap_DC_CYP/
#     CYP_phylo_heatmap.pdf/png       — tree + heatmap combined
#     CYP_clustered_heatmap.pdf/png   — expression-clustered heatmap
#     msa_alignment.fasta             — MAFFT alignment (if tree was built)
#     phylo_tree.nwk                  — Newick tree (if tree was built, reusable)
#
################################################################################

#SBATCH --job-name=phylo_heatmap
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu

#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis
#SBATCH -o 06_analysis/phylo_heatmap_%j.out
#SBATCH -e 06_analysis/phylo_heatmap_%j.err

set -euo pipefail

BASE_DIR="/projects/tholl_lab_1/daisy_analysis"
CODE_DIR="${BASE_DIR}/05_rnaseq-code"

# ── Arguments ────────────────────────────────────────────────────────────────
SPECIES="${1:-DC}"
FAMILY_KEY="${2:-cyp}"
SPECIES=$(echo "${SPECIES}" | tr '[:lower:]' '[:upper:]')
FAMILY_KEY=$(echo "${FAMILY_KEY}" | tr '[:upper:]' '[:lower:]')

# ── Settings (all overridable via env vars) ──────────────────────────────────
SUBSET="${SUBSET:-}"
TREE="${TREE:-}"
CLUSTER_ONLY="${CLUSTER_ONLY:-0}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"

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

SUFFIX="${FAMILY_LABEL}"
[ -n "${SUBSET}" ] && SUFFIX="${SUFFIX}_${SUBSET}"
OUTPUT_DIR="${BASE_DIR}/06_analysis/phylo_heatmap_${SPECIES}_${SUFFIX}"

# ── Activate environment ─────────────────────────────────────────────────────
source ~/.bashrc
conda activate rnaseq
module load MAFFT/7.505-GCC-11.3.0-with-extensions 2>/dev/null || true
module load FastTree/2.1.11-GCCcore-11.3.0 2>/dev/null || true

echo "============================================================"
echo "Phylo + Heatmap — ${FAMILY_LABEL} (${SPECIES})"
echo "============================================================"
echo "  Expressed list: ${EXPRESSED}"
echo "  Protein FASTA:  ${FASTA}"
echo "  Count matrix:   ${COUNT_MATRIX}"
echo "  Output:         ${OUTPUT_DIR}"
[ -n "${SUBSET}" ] && echo "  Subset:         ${SUBSET}"
[ -n "${TREE}" ] && echo "  Tree:           ${TREE}"
[ "${CLUSTER_ONLY}" = "1" ] && echo "  Mode:           cluster-only (no tree)"
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

[ "${MISSING}" = "1" ] && exit 1

mkdir -p "${OUTPUT_DIR}"

# ── Build command ────────────────────────────────────────────────────────────
CMD=(python3 "${CODE_DIR}/scripts/10_post_analysis/phylo_heatmap.py"
     --counts "${COUNT_MATRIX}"
     --metadata "${METADATA}"
     --expression "${EXPRESSED}"
     --output-dir "${OUTPUT_DIR}"
     --family "${FAMILY_LABEL}"
     --threads "${THREADS}")

if [ "${CLUSTER_ONLY}" = "1" ]; then
    CMD+=(--cluster-only)
elif [ -n "${TREE}" ] && [ -f "${TREE}" ]; then
    CMD+=(--tree "${TREE}")
elif [ -f "${FASTA}" ]; then
    CMD+=(--fasta "${FASTA}")
else
    echo "  WARNING: Protein FASTA not found (${FASTA})"
    echo "  Running in cluster-only mode (no phylo tree)"
    CMD+=(--cluster-only)
fi

[ -n "${SUBSET}" ] && CMD+=(--subset "${SUBSET}")

echo "  Running: ${CMD[*]}"
echo ""
"${CMD[@]}"

echo ""
echo "============================================================"
echo "COMPLETE — Phylo + Heatmap"
echo "============================================================"
echo "Output: ${OUTPUT_DIR}/"
ls -la "${OUTPUT_DIR}/"
echo ""
echo "Download:"
echo "  scp -r daisycortesj@tinkercliffs.arc.vt.edu:${OUTPUT_DIR}/ ~/Desktop/phylo_heatmap/"
echo ""
echo "End time: $(date)"
