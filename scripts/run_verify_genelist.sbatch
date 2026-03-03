#!/bin/bash
################################################################################
# Verify step 3 results against P450_list_RefSeq.txt
################################################################################
#
# Double-checks that every gene in your step 3 output actually appears in
# the original P450 gene list. Run AFTER step 3 plots finish.
#
# How to run (defaults check step 3 geneious output):
#   sbatch scripts/run_verify_genelist.sbatch
#
# Check geneious_candidates.tsv (Path A output, before step 3):
#   RESULTS=07_NRdatabase/sukman_database/geneious_candidates.tsv \
#     sbatch scripts/run_verify_genelist.sbatch
#
# Check cyp_strict step 3 output:
#   RESULTS=06_analysis/pydeseq2_DC_step3_plots_cyp_strict_padj005_lfc20/cyp_gene_list.tsv \
#     sbatch scripts/run_verify_genelist.sbatch
#
# Use a different gene list:
#   GENE_LIST=07_NRdatabase/cyp450_database/cyp_master_list.csv \
#     sbatch scripts/run_verify_genelist.sbatch
#
# Output: prints PASS/FAIL to the job log — no files created
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
#SBATCH -o 07_NRdatabase/verify_genelist_%j.out
#SBATCH -e 07_NRdatabase/verify_genelist_%j.err

set -euo pipefail

source ~/.bashrc
conda activate rnaseq

# ============================================================================
# PATHS
# ============================================================================
BASE_DIR="/projects/tholl_lab_1/daisy_analysis"
CODE_DIR="${BASE_DIR}/05_rnaseq-code"

GENE_LIST="${GENE_LIST:-${BASE_DIR}/07_NRdatabase/sukman_database/P450_list_RefSeq.txt}"
RESULTS="${RESULTS:-${BASE_DIR}/06_analysis/pydeseq2_DC_step3_plots_geneious_candidates_padj005_lfc20/cyp_gene_list.tsv}"

# ============================================================================
# VALIDATION
# ============================================================================
echo "========================================================================"
echo "  Verify Gene List — Double Check"
echo "  $(date)"
echo "========================================================================"
echo ""
echo "  Results file: ${RESULTS}"
echo "  Gene list:    ${GENE_LIST}"
echo ""

if [ ! -f "$GENE_LIST" ]; then
    echo "ERROR: Gene list not found: $GENE_LIST"
    exit 1
fi

if [ ! -f "$RESULTS" ]; then
    echo "ERROR: Results file not found: $RESULTS"
    echo ""
    echo "  Did step 3 finish? Check with:"
    echo "    ls -la 06_analysis/pydeseq2_DC_step3_plots_*/"
    exit 1
fi

# ============================================================================
# RUN
# ============================================================================
python "${CODE_DIR}/scripts/verify_genelist.py" \
    --results "$RESULTS" \
    --gene-list "$GENE_LIST"

echo ""
echo "Done: $(date)"
