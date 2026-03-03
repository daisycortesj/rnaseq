#!/bin/bash
################################################################################
# CYP450 PIPELINE — Step 1: Intersect + Step 2: Extract proteins
################################################################################
#
# Prerequisite:
#   Run run_cyp450_database.sbatch FIRST to produce cyp_master_list.csv
#
# What this does (2 steps):
#   1. Intersect CYP master list with PyDESeq2 results → cyp_expressed_list.tsv
#   2. Extract CYP protein sequences from the full FASTA → cyp_proteins.fasta
#
# After this job finishes, continue with:
#   Step 3 — BLAST:
#     sbatch scripts/blastp_discoveryfilter.sbatch DC swissprot \
#         /projects/tholl_lab_1/daisy_analysis/07_NRdatabase/cyp450_database/cyp_proteins.fasta
#
# How to run:
#   sbatch scripts/run_cyp_express_extract.sbatch
#
# Output (in 07_NRdatabase/cyp450_database/):
#   cyp_expressed_list.tsv    — CYPs in count matrix with expression stats
#   cyp_proteins.fasta        — subset protein sequences (gene_id headers)
#
################################################################################

#SBATCH --job-name=cyp_extract
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=30:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu

#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis
#SBATCH -o 07_NRdatabase/cyp_extract_%j.out
#SBATCH -e 07_NRdatabase/cyp_extract_%j.err

set -euo pipefail

source ~/.bashrc
conda activate rnaseq

# ============================================================================
# PATHS
# ============================================================================
BASE_DIR="/projects/tholl_lab_1/daisy_analysis"
CODE_DIR="${BASE_DIR}/05_rnaseq-code"
OUTDIR="${BASE_DIR}/07_NRdatabase/cyp450_database"

CYP_LIST="${OUTDIR}/cyp_master_list.csv"
GFF="${BASE_DIR}/04_reference/dc_genomic.gtf"
PROTEIN_FASTA="${BASE_DIR}/04_reference/GCF_001625215.2_DH1_v3.0_protein.faa"
DESEQ="${BASE_DIR}/06_analysis/pydeseq2_DC_step1_unfiltered/pydeseq2_results_UNFILTERED.tsv"

EXPRESSED="${OUTDIR}/cyp_expressed_list.tsv"
CYP_FASTA="${OUTDIR}/cyp_proteins.fasta"

# ============================================================================
# VALIDATION
# ============================================================================
echo "========================================================================"
echo "  CYP450 PIPELINE — Intersect + Extract"
echo "  Started: $(date)"
echo "========================================================================"

for f in "$CYP_LIST" "$PROTEIN_FASTA" "$DESEQ" "$GFF"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Required file not found: $f"
        exit 1
    fi
done

mkdir -p "$OUTDIR"
cd "$BASE_DIR"

# ============================================================================
# STEP 1: Intersect CYP master list with PyDESeq2
# ============================================================================
echo ""
echo "STEP 1/2: Intersecting CYP master list with PyDESeq2..."
echo "------------------------------------------------------------------------"

python "${CODE_DIR}/scripts/cyp_intersect_pydeseq2.py" \
    --cyp-list "$CYP_LIST" \
    --deseq "$DESEQ" \
    --gff "$GFF" \
    --output "$EXPRESSED"

echo ""

# ============================================================================
# STEP 2: Extract CYP protein sequences
# ============================================================================
echo "STEP 2/2: Extracting CYP protein sequences..."
echo "------------------------------------------------------------------------"

python "${CODE_DIR}/scripts/cyp_extract_proteins.py" \
    --expressed "$EXPRESSED" \
    --fasta "$PROTEIN_FASTA" \
    --output "$CYP_FASTA"

N_PROTS=$(grep -c "^>" "$CYP_FASTA" || true)
echo "  Proteins in subset FASTA: $N_PROTS"

echo ""
echo "========================================================================"
echo "  DONE! $(date)"
echo "  Results in: ${OUTDIR}/"
echo ""
echo "  Output files:"
echo "    ${EXPRESSED}     — CYP genes with expression stats"
echo "    ${CYP_FASTA}     — subset protein FASTA ($N_PROTS sequences)"
echo ""
echo "  NEXT STEP — Run BLAST on the CYP proteins:"
echo "    sbatch scripts/blastp_discoveryfilter.sbatch DC swissprot ${CYP_FASTA}"
echo "========================================================================"
