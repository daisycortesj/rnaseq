#!/bin/bash
################################################################################
# RUN R_pydeq2.py 
################################################################################
#
# PURPOSE:
#   Runs the complete student R DESeq2 workflow in Python — one script:
#   Load counts → PyDESeq2 → QC plots → filter significant → heatmaps →
#   volcano → MA → P450/MT family analysis → phylo trees → NCBI fetch →
#   genomic clustering
#
# HOW TO RUN:
#   # Default (DC, Root vs Leaf — positive log2FC = up in Root):
#   sbatch scripts/R_pydeseq2/run_R_pydeq2.sbatch DC
#
#   # Run DG separately:
#   sbatch scripts/R_pydeseq2/run_R_pydeq2.sbatch DG
#
#   # Combined DC + DG (auto-merges count matrices, more replicates):
#   sbatch scripts/R_pydeseq2/run_R_pydeq2.sbatch DCDG
#
#   # With custom contrast (e.g. Ethylene vs Normal — Normal is baseline):
#   CONTRAST_A=E CONTRAST_B=N sbatch scripts/R_pydeseq2/run_R_pydeq2.sbatch DC
#
#   # Skip slow steps (phylo tree / NCBI):
#   SKIP_PHYLO=1 SKIP_NCBI=1 sbatch scripts/R_pydeseq2/run_R_pydeq2.sbatch DC
#
#   # With methyltransferase list:
#   MT_LIST=07_NRdatabase/omt_database/omt_list.txt \
#     sbatch scripts/R_pydeseq2/run_R_pydeq2.sbatch DC
#
#   # Custom cutoffs:
#   PADJ=0.01 LFC=3.0 MIN_COUNTS=10 sbatch scripts/R_pydeseq2/run_R_pydeq2.sbatch DC
#
# PACKAGES LOADED:
#   conda env:  rnaseq  (numpy, pandas, matplotlib, seaborn, scipy, pydeseq2, biopython)
#   modules:    MAFFT, FastTree  (for phylogenetic tree building)
#
# OUTPUT:
#   06_analysis/R_pydeseq2_DC/
#     normalized_counts_refseq.csv    — normalized counts
#     res_all.csv                     — all DE results
#     res_all_significant.csv         — padj ≤ 0.05 & |log2FC| ≥ 2
#     res_Upregulated.csv             — upregulated genes only
#     res_Downregulated.csv           — downregulated genes only
#     pca_plot.pdf/png                — PCA (like vst + plotPCA)
#     dispersion_plot.pdf/png         — dispersion (like plotDispEsts)
#     volcano_plot.pdf/png            — volcano (like EnhancedVolcano)
#     ma_plot.pdf/png                 — MA plot (like plotMA)
#     heatmap_all_significant_genes.* — all sig DE heatmap
#     heatmap_P450.*                  — P450 family heatmap (all ~50 genes)
#     heatmap_P450_upregulated.*      — P450 upregulated (from this analysis)
#     heatmap_P450_DE_prev_student.*  — P450 DE genes (previous student's list)
#     heatmap_P450_upregulated_prev_student.*  — P450 upregulated (previous student)
#     heatmap_P450_downregulated_prev_student.* — P450 downregulated (previous student)
#     heatmap_MT.*                    — methyltransferase heatmap
#     heatmap_phylogeny_P450.*        — phylo tree + expression
#     P450_genomic_locations_description.csv — NCBI gene info
#     Genes_with_distance_P450.csv    — genomic proximity
#     P450_distance_*.png             — distance plots
#
################################################################################

#SBATCH --job-name=R_pydeq2
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu

#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis
#SBATCH -o 06_analysis/R_pydeq2_%j.out
#SBATCH -e 06_analysis/R_pydeq2_%j.err

set -euo pipefail

# ===== Load shared config =====
source "${SLURM_SUBMIT_DIR:-.}/scripts/config.sh"

SPECIES="${1:-DC}"
SPECIES=$(echo "${SPECIES}" | tr '[:lower:]' '[:upper:]')

# ── Settings (all overridable via env vars) ──────────────────────────────────
CONTRAST_A="${CONTRAST_A:-R}"
CONTRAST_B="${CONTRAST_B:-L}"
CONTRAST_FACTOR="${CONTRAST_FACTOR:-condition}"
MIN_COUNTS="${MIN_COUNTS:-20}"
PADJ="${PADJ:-0.05}"
LFC="${LFC:-2.0}"
THRESHOLD="${THRESHOLD:-50000}"
EMAIL="${EMAIL:-daisycortesj@vt.edu}"
SKIP_PHYLO="${SKIP_PHYLO:-0}"
SKIP_NCBI="${SKIP_NCBI:-0}"

# ── Resolve species-specific paths ───────────────────────────────────────────
get_sample_info "${SPECIES}" || exit 1
COUNT_FOLDER="${SPECIES_DIR}"

# For combined species (DCDG), merge count matrices first
if [[ "${SPECIES}" == "DCDG" ]]; then
    MERGE_DIR="${BASE_DIR}/03_count_tables/${COUNT_FOLDER}"
    DC_DIR="${BASE_DIR}/03_count_tables/00_1_DC"
    DG_DIR="${BASE_DIR}/03_count_tables/00_2_DG"

    echo "── Merging DC + DG count matrices ──"
    python3 "${CODE_DIR}/scripts/05_pydeseq2/R_pydeseq2/merge_count_matrices.py" \
        --dirs "${DC_DIR}" "${DG_DIR}" \
        --output-dir "${MERGE_DIR}"
    echo ""
fi

COUNT_MATRIX="${BASE_DIR}/03_count_tables/${COUNT_FOLDER}/gene_count_matrix.tsv"
METADATA="${BASE_DIR}/03_count_tables/${COUNT_FOLDER}/sample_metadata.tsv"
GTF="${BASE_DIR}/04_reference/dc_genomic.gtf"
PROTEIN_FASTA="${BASE_DIR}/04_reference/GCF_001625215.2_DH1_v3.0_protein.faa"
P450_LIST="${P450_LIST:-${BASE_DIR}/07_NRdatabase/sukman_database/P450_list_RefSeq.txt}"
MT_LIST="${MT_LIST:-${BASE_DIR}/07_NRdatabase/sukman_database/MT_expression_refseq.txt}"
OUTPUT_DIR="${BASE_DIR}/06_analysis/R_pydeseq2_${SPECIES}"

# ── Load packages / activate environment ─────────────────────────────────────
source ~/.bashrc
conda activate rnaseq

# Load external tools for phylogenetic analysis (if available on HPC)
module load MAFFT/7.505-GCC-11.3.0-with-extensions 2>/dev/null || true
module load FastTree/2.1.11-GCCcore-11.3.0 2>/dev/null || true

echo "============================================================"
echo "R_pydeq2 — Student R Analysis → Python"
echo "============================================================"
echo "Species:  ${SPECIES}"
echo "Contrast: ${CONTRAST_A} vs ${CONTRAST_B}"
echo "Cutoffs:  padj ≤ ${PADJ}, |log2FC| ≥ ${LFC}, min_counts ≥ ${MIN_COUNTS}"
echo "Output:   ${OUTPUT_DIR}"
echo ""

# ── Show loaded package versions ─────────────────────────────────────────────
echo "── Package versions ──"
python3 -c "
import numpy, pandas, matplotlib, seaborn, scipy, pydeseq2
print(f'  numpy:      {numpy.__version__}')
print(f'  pandas:     {pandas.__version__}')
print(f'  matplotlib: {matplotlib.__version__}')
print(f'  seaborn:    {seaborn.__version__}')
print(f'  scipy:      {scipy.__version__}')
print(f'  pydeseq2:   {pydeseq2.__version__}')
try:
    import Bio; print(f'  biopython:  {Bio.__version__}')
except ImportError:
    print('  biopython:  NOT INSTALLED (NCBI steps will be skipped)')
"
echo ""

# ── Validate inputs ──────────────────────────────────────────────────────────
for f in "${COUNT_MATRIX}" "${METADATA}"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Required file not found: $f"
        exit 1
    fi
done

mkdir -p "${OUTPUT_DIR}"

# ── Build and run command ────────────────────────────────────────────────────
CMD=(python3 "${CODE_DIR}/scripts/05_pydeseq2/R_pydeseq2/R_pydeq2.py"
     --counts "${COUNT_MATRIX}"
     --metadata "${METADATA}"
     --gtf "${GTF}"
     --protein-fasta "${PROTEIN_FASTA}"
     --output-dir "${OUTPUT_DIR}"
     --contrast-factor "${CONTRAST_FACTOR}"
     --contrast-A "${CONTRAST_A}"
     --contrast-B "${CONTRAST_B}"
     --min-counts "${MIN_COUNTS}"
     --padj "${PADJ}"
     --lfc "${LFC}"
     --threshold "${THRESHOLD}"
     --email "${EMAIL}")

[ -f "${P450_LIST}" ] && CMD+=(--p450-list "${P450_LIST}")
[ -n "${MT_LIST}" ] && [ -f "${MT_LIST}" ] && CMD+=(--mt-list "${MT_LIST}")
[ "${SKIP_PHYLO}" = "1" ] && CMD+=(--skip-phylo)
[ "${SKIP_NCBI}" = "1" ] && CMD+=(--skip-ncbi)

echo "Running: ${CMD[*]}"
echo ""

"${CMD[@]}"

echo ""
echo "============================================================"
echo "COMPLETE"
echo "============================================================"
echo "Output: ${OUTPUT_DIR}/"
echo ""
echo "Download to laptop:"
echo "  scp -r daisycortesj@tinkercliffs.arc.vt.edu:${OUTPUT_DIR}/ ~/Desktop/R_pydeseq2_results/"
echo ""
echo "End time: $(date)"
