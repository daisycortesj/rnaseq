#!/bin/bash
################################################################################
# Gene-Family PIPELINE — Step 1: Intersect + Step 2: Extract proteins
################################################################################
#
# Supports: CYP (cytochrome P450) and OMT (O-methyltransferase)
#
# Prerequisite:
#   Run the appropriate database builder FIRST:
#     CYP: sbatch scripts/run_cyp450_database.sbatch  → cyp_master_list.csv
#     OMT: sbatch scripts/run_omt_database.sbatch     → omt_master_list.csv
#
# What this does (2 steps, run separately for each list):
#   1. Intersect gene list with PyDESeq2 results → expressed list TSV
#   2. Extract protein sequences from the full FASTA → proteins FASTA
#
# How to run (FAMILY is required: cyp or omt):
#   sbatch scripts/run_gene_extract.sbatch cyp
#   sbatch scripts/run_gene_extract.sbatch omt
#
# How to run (also run Geneious list separately):
#   PREV_LIST=07_NRdatabase/sukman_database/P450_list_RefSeq.txt \
#     sbatch scripts/run_gene_extract.sbatch cyp
#
# Output (in 07_NRdatabase/<family>_database/):
#   <family>_expressed_list.tsv       — YOUR genes in count matrix with expression stats
#   <family>_proteins.fasta           — YOUR subset protein sequences
#   geneious_expressed_list.tsv       — Geneious list in count matrix (only if PREV_LIST set)
#   geneious_proteins.fasta           — Geneious subset protein sequences (only if PREV_LIST set)
#
# CHANGED: 2026-03-31 — Renamed from run_cyp_extract.sbatch → run_gene_extract.sbatch.
#   Now accepts FAMILY arg (cyp or omt) instead of being CYP-only.
#   All CYP-hardcoded paths replaced with FAMILY-driven paths.
#   Python script calls updated: cyp_intersect_pydeseq2.py → gene_intersect_pydeseq2.py,
#                                 cyp_extract_proteins.py → gene_extract_proteins.py.
################################################################################

#SBATCH --job-name=family_extract
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=30:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu

#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis
#SBATCH -o 07_NRdatabase/family_extract_%j.out
#SBATCH -e 07_NRdatabase/family_extract_%j.err

set -euo pipefail

# ===== Load shared config =====
source "${SLURM_SUBMIT_DIR:-.}/scripts/config.sh"

source ~/.bashrc
conda activate rnaseq

# ============================================================================
# FAMILY ARGUMENT (required: cyp or omt)
# ============================================================================
FAMILY="${1:-}"

if [ -z "${FAMILY}" ]; then
    echo "ERROR: You need to specify a gene family!"
    echo ""
    echo "Usage:"
    echo "  sbatch scripts/run_gene_extract.sbatch cyp"
    echo "  sbatch scripts/run_gene_extract.sbatch omt"
    echo ""
    echo "Optional: also run a Geneious list:"
    echo "  PREV_LIST=07_NRdatabase/sukman_database/P450_list_RefSeq.txt \\"
    echo "    sbatch scripts/run_gene_extract.sbatch cyp"
    exit 1
fi

# Normalize to lowercase for paths, uppercase for labels
FAMILY=$(echo "${FAMILY}" | tr '[:upper:]' '[:lower:]')
FAMILY_UPPER=$(echo "${FAMILY}" | tr '[:lower:]' '[:upper:]')

if [ "${FAMILY}" != "cyp" ] && [ "${FAMILY}" != "omt" ]; then
    echo "ERROR: Unknown gene family: ${FAMILY}"
    echo "  Supported families: cyp, omt"
    exit 1
fi

# ============================================================================
# PATHS — set based on FAMILY
# ============================================================================

# Family-specific directory and files
if [ "${FAMILY}" = "cyp" ]; then
    OUTDIR="${BASE_DIR}/07_NRdatabase/cyp450_database"
    MASTER_LIST="${OUTDIR}/cyp_master_list.csv"
    FAMILY_LABEL="CYP450 (cytochrome P450)"
elif [ "${FAMILY}" = "omt" ]; then
    OUTDIR="${BASE_DIR}/07_NRdatabase/omt_database"
    MASTER_LIST="${OUTDIR}/omt_master_list.csv"
    FAMILY_LABEL="OMT (O-methyltransferase)"
fi

GFF="${BASE_DIR}/04_reference/dc_genomic.gtf"
PROTEIN_FASTA="${BASE_DIR}/04_reference/GCF_001625215.2_DH1_v3.0_protein.faa"
DESEQ="${BASE_DIR}/06_analysis/pydeseq2_DC_step1_unfiltered/pydeseq2_results_UNFILTERED.tsv"

# Optional: set PREV_LIST to also run the Geneious list separately
# e.g.  PREV_LIST=07_NRdatabase/sukman_database/P450_list_RefSeq.txt sbatch scripts/run_gene_extract.sbatch cyp
PREV_LIST="${PREV_LIST:-}"

# Output file names use the family prefix (cyp_expressed_list.tsv, omt_expressed_list.tsv, etc.)
EXPRESSED="${OUTDIR}/${FAMILY}_expressed_list.tsv"
FAMILY_FASTA="${OUTDIR}/${FAMILY}_proteins.fasta"

GENEIOUS_EXPRESSED="${OUTDIR}/geneious_expressed_list.tsv"
GENEIOUS_FASTA="${OUTDIR}/geneious_proteins.fasta"

# ============================================================================
# VALIDATION
# ============================================================================
echo "========================================================================"
echo "  ${FAMILY_UPPER} PIPELINE — Intersect + Extract"
echo "  Family: ${FAMILY_LABEL}"
echo "  Started: $(date)"
echo "========================================================================"

for f in "$MASTER_LIST" "$PROTEIN_FASTA" "$DESEQ" "$GFF"; do
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
echo "  Intersecting your ${FAMILY_UPPER} master list with PyDESeq2..."

# CHANGED: 2026-03-31 — renamed cyp_intersect_pydeseq2.py → gene_intersect_pydeseq2.py
python "${CODE_DIR}/scripts/08_gene_families/gene_intersect_pydeseq2.py" \
    --gene-list "$MASTER_LIST" \
    --deseq "$DESEQ" \
    --gff "$GFF" \
    --output "$EXPRESSED" \
    --database "$FAMILY_UPPER"

echo ""
echo "  Extracting protein sequences..."

# CHANGED: 2026-03-31 — renamed cyp_extract_proteins.py → gene_extract_proteins.py
python "${CODE_DIR}/scripts/08_gene_families/gene_extract_proteins.py" \
    --expressed "$EXPRESSED" \
    --fasta "$PROTEIN_FASTA" \
    --output "$FAMILY_FASTA"

N_YOURS=$(grep -c "^>" "$FAMILY_FASTA" || true)
echo "  Your proteins: $N_YOURS"

# ============================================================================
# GENEIOUS LIST — Intersect + Extract (only if PREV_LIST is set)
# ============================================================================
if [ -n "$PREV_LIST" ]; then
    echo ""
    echo "──── GENEIOUS LIST ──────────────────────────────────────────────────────"
    echo "  Intersecting Geneious list with the SAME PyDESeq2 results..."

    python "${CODE_DIR}/scripts/08_gene_families/gene_intersect_pydeseq2.py" \
        --gene-list "$PREV_LIST" \
        --deseq "$DESEQ" \
        --gff "$GFF" \
        --output "$GENEIOUS_EXPRESSED" \
        --database "$FAMILY_UPPER" \
        --evidence all

    echo ""
    echo "  Extracting protein sequences..."

    python "${CODE_DIR}/scripts/08_gene_families/gene_extract_proteins.py" \
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
echo "  Family: ${FAMILY_LABEL}"
echo "  Results in: ${OUTDIR}/"
echo ""
echo "  Your list:"
echo "    ${EXPRESSED}"
echo "    ${FAMILY_FASTA}  ($N_YOURS proteins)"
if [ -n "$PREV_LIST" ]; then
    echo ""
    echo "  Geneious list:"
    echo "    ${GENEIOUS_EXPRESSED}"
    echo "    ${GENEIOUS_FASTA}  ($N_GENEIOUS proteins)"
fi
echo ""
echo "  NEXT — go straight to plots (short path B):"
echo "    sbatch scripts/05_pydeseq2/run_step3_plots.sbatch DC ${FAMILY}_expressed"
if [ -n "$PREV_LIST" ]; then
    echo "    sbatch scripts/05_pydeseq2/run_step3_plots.sbatch DC ${GENEIOUS_EXPRESSED}"
fi
echo ""
echo "  OR — run BLAST first (full path C):"
echo "    sbatch scripts/blastp_discoveryfilter.sbatch DC swissprot ${FAMILY_FASTA}"
echo "    sbatch scripts/06_blast/run_combine_blast_deseq.sbatch DC swissprot discovery standard ${FAMILY}"
echo "    sbatch scripts/05_pydeseq2/run_step3_plots.sbatch DC ${FAMILY}_discovery"
echo "========================================================================"
