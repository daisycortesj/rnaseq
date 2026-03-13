#!/bin/bash
################################################################################
# BUILD CYP450 DATABASE — GFF keyword search + HMMER domain scan + Combine
################################################################################
#
# What this does:
#   Step 1: Searches the GTF annotation for CYP450 keywords
#   Step 2: Runs HMMER domain scan (PF00067) against all carrot proteins
#   Step 3: Combines both lists into a master CYP450 gene database
#
# How to run:
#   sbatch scripts/run_cyp450_database.sbatch
#
# What you get (in 08_cyp450_database/):
#   cyp_gtf_candidates.csv   — genes found by keyword search
#   cyp_hmmer_confirmed.csv  — genes confirmed by P450 domain
#   cyp_master_list.csv      — combined master list with evidence flags
#
################################################################################


# ══════════════════════════════════════════════════════════════════════════════
# SLURM SETTINGS
# ══════════════════════════════════════════════════════════════════════════════

#SBATCH --job-name=cyp450_db
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daisycortesj@vt.edu

#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis
#SBATCH -o 07_NRdatabase/cyp450_db_%j.out
#SBATCH -e 07_NRdatabase/cyp450_db_%j.err


# ══════════════════════════════════════════════════════════════════════════════
# PATHS
# ══════════════════════════════════════════════════════════════════════════════

BASE_DIR="/projects/tholl_lab_1/daisy_analysis"
CODE_DIR="${BASE_DIR}/05_rnaseq-code"

set -e

# ══════════════════════════════════════════════════════════════════════════════
# SETUP
# ══════════════════════════════════════════════════════════════════════════════

source ~/.bashrc
conda activate rnaseq

cd "${BASE_DIR}"
mkdir -p 07_NRdatabase/cyp450_database

echo "========================================================================"
echo "  CYP450 DATABASE BUILD"
echo "  Started: $(date)"
echo "  Base dir: ${BASE_DIR}"
echo "========================================================================"
echo ""

# ══════════════════════════════════════════════════════════════════════════════
# STEP 1: GFF KEYWORD SEARCH
# ══════════════════════════════════════════════════════════════════════════════

echo "STEP 1: GFF keyword search..."
echo "------------------------------------------------------------------------"

python "${CODE_DIR}/scripts/08_gene_families/cyp_gff_keyword_search.py" \
    --gff "${BASE_DIR}/04_reference/dc_genomic.gtf" \
    --output "${BASE_DIR}/07_NRdatabase/cyp450_database/cyp_gtf_candidates.csv"

echo ""

# ══════════════════════════════════════════════════════════════════════════════
# STEP 2: HMMER DOMAIN SCAN
# ══════════════════════════════════════════════════════════════════════════════

echo "STEP 2: HMMER domain scan..."
echo "------------------------------------------------------------------------"

python "${CODE_DIR}/scripts/08_gene_families/cyp_hmmer_scan.py" \
    --proteins "${BASE_DIR}/04_reference/GCF_001625215.2_DH1_v3.0_protein.faa" \
    --output "${BASE_DIR}/07_NRdatabase/cyp450_database/cyp_hmmer_confirmed.csv" \
    --evalue 1e-5 \
    --cpus 4

echo ""

# ══════════════════════════════════════════════════════════════════════════════
# STEP 3: COMBINE RESULTS
# ══════════════════════════════════════════════════════════════════════════════

echo "STEP 3: Combining results..."
echo "------------------------------------------------------------------------"

python "${CODE_DIR}/scripts/08_gene_families/cyp_combine_results.py" \
    --gtf-csv "${BASE_DIR}/07_NRdatabase/cyp450_database/cyp_gtf_candidates.csv" \
    --hmmer-csv "${BASE_DIR}/07_NRdatabase/cyp450_database/cyp_hmmer_confirmed.csv" \
    --gff "${BASE_DIR}/04_reference/dc_genomic.gtf" \
    --output "${BASE_DIR}/07_NRdatabase/cyp450_database/cyp_master_list.csv"

echo ""
echo "========================================================================"
echo "  DONE! $(date)"
echo "  Results in: ${BASE_DIR}/07_NRdatabase/cyp450_database/"
echo "========================================================================"
