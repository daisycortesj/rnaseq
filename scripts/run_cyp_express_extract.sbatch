#!/bin/bash
################################################################################
# CYP450 PIPELINE — Step 1: Intersect + Step 2: Extract proteins
################################################################################
#
# Prerequisite:
#   Run run_cyp450_database.sbatch FIRST to produce cyp_master_list.csv
#
# What this does (2 steps, run separately for each list):
#   1. Intersect gene list with PyDESeq2 results → expressed list TSV
#   2. Extract protein sequences from the full FASTA → proteins FASTA
#
# How to run (your list only — default):
#   sbatch scripts/run_cyp_express_extract.sbatch
#
# How to run (also run Geneious list separately):
#   PREV_LIST=07_NRdatabase/sukman_database/P450_list_RefSeq.txt \
#     sbatch scripts/run_cyp_express_extract.sbatch
#
# Output (in 07_NRdatabase/cyp450_database/):
#   cyp_expressed_list.tsv       — YOUR CYPs in count matrix with expression stats
#   cyp_proteins.fasta           — YOUR subset protein sequences
#   geneious_expressed_list.tsv  — Geneious P450s in count matrix (only if PREV_LIST set)
#   geneious_proteins.fasta      — Geneious subset protein sequences (only if PREV_LIST set)
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

# Optional: set PREV_LIST to also run the Geneious list separately
# e.g.  PREV_LIST=07_NRdatabase/sukman_database/P450_list_RefSeq.txt sbatch scripts/run_cyp_express_extract.sbatch
PREV_LIST="${PREV_LIST:-}"

EXPRESSED="${OUTDIR}/cyp_expressed_list.tsv"
CYP_FASTA="${OUTDIR}/cyp_proteins.fasta"

GENEIOUS_EXPRESSED="${OUTDIR}/geneious_expressed_list.tsv"
GENEIOUS_FASTA="${OUTDIR}/geneious_proteins.fasta"

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

if [ -n "$PREV_LIST" ] && [ ! -f "$PREV_LIST" ]; then
    echo "ERROR: Previous student list not found: $PREV_LIST"
    exit 1
fi

mkdir -p "$OUTDIR"
cd "$BASE_DIR"

# ============================================================================
# YOUR LIST — Intersect + Extract
# ============================================================================
echo ""
echo "──── YOUR HMMER/GTF LIST ────────────────────────────────────────────────"
echo "  Intersecting your CYP master list with PyDESeq2..."

python "${CODE_DIR}/scripts/cyp_intersect_pydeseq2.py" \
    --gene-list "$CYP_LIST" \
    --deseq "$DESEQ" \
    --gff "$GFF" \
    --output "$EXPRESSED" \
    --database CYP

echo ""
echo "  Extracting protein sequences..."

python "${CODE_DIR}/scripts/cyp_extract_proteins.py" \
    --expressed "$EXPRESSED" \
    --fasta "$PROTEIN_FASTA" \
    --output "$CYP_FASTA"

N_YOURS=$(grep -c "^>" "$CYP_FASTA" || true)
echo "  Your proteins: $N_YOURS"

# ============================================================================
# GENEIOUS LIST — Intersect + Extract (only if PREV_LIST is set)
# ============================================================================
if [ -n "$PREV_LIST" ]; then
    echo ""
    echo "──── GENEIOUS LIST ──────────────────────────────────────────────────────"
    echo "  Intersecting Geneious list with the SAME PyDESeq2 results..."

    python "${CODE_DIR}/scripts/cyp_intersect_pydeseq2.py" \
        --gene-list "$PREV_LIST" \
        --deseq "$DESEQ" \
        --gff "$GFF" \
        --output "$GENEIOUS_EXPRESSED" \
        --database CYP \
        --evidence all

    echo ""
    echo "  Extracting protein sequences..."

    python "${CODE_DIR}/scripts/cyp_extract_proteins.py" \
        --expressed "$GENEIOUS_EXPRESSED" \
        --fasta "$PROTEIN_FASTA" \
        --output "$GENEIOUS_FASTA"

    N_GENEIOUS=$(grep -c "^>" "$GENEIOUS_FASTA" || true)
    echo "  Geneious proteins: $N_GENEIOUS"
fi

# ============================================================================
# SUMMARY
# ============================================================================
echo ""
echo "========================================================================"
echo "  DONE! $(date)"
echo "  Results in: ${OUTDIR}/"
echo ""
echo "  Your list:"
echo "    ${EXPRESSED}"
echo "    ${CYP_FASTA}  ($N_YOURS proteins)"
if [ -n "$PREV_LIST" ]; then
    echo ""
    echo "  Geneious list:"
    echo "    ${GENEIOUS_EXPRESSED}"
    echo "    ${GENEIOUS_FASTA}  ($N_GENEIOUS proteins)"
fi
echo ""
echo "  NEXT — go straight to plots (short path):"
echo "    sbatch scripts/run_pydeseq2_step3_plots.sbatch DC cyp_expressed"
if [ -n "$PREV_LIST" ]; then
    echo "    sbatch scripts/run_pydeseq2_step3_plots.sbatch DC ${GENEIOUS_EXPRESSED}"
fi
echo ""
echo "  OR — run BLAST first (full path):"
echo "    sbatch scripts/blastp_discoveryfilter.sbatch DC swissprot ${CYP_FASTA}"
echo "========================================================================"
