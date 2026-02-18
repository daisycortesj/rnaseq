#!/bin/bash
################################################################################
# Compare differential expression between two species
################################################################################
#
# Merges annotated results from two species into a single comparison table.
# Both species must use the same reference genome (same gene IDs).
#
# Requires:
#   - Annotated results for both species (from run_combine_filter.sbatch)
#
# How to run:
#   sbatch scripts/run_compare_species.sbatch DC DG swissprot discovery
#   sbatch scripts/run_compare_species.sbatch DC DG nr discovery
#   sbatch scripts/run_compare_species.sbatch DC DG swissprot discovery 0.01 2.5
#
# Arguments:
#   $1 = Species 1 (e.g., DC)
#   $2 = Species 2 (e.g., DG)
#   $3 = Database (e.g., swissprot, nr)
#   $4 = BLAST mode (discovery or strict)
#   $5 = padj cutoff (optional, default: 0.05)
#   $6 = log2FC cutoff (optional, default: 2.0)
#
# Output:
#   06_analysis/comparison_{SP1}_vs_{SP2}/
#     {SP1}_vs_{SP2}_{DB}_{MODE}_comparison.tsv        (full table)
#     {SP1}_vs_{SP2}_{DB}_{MODE}_comparison_shared_same_direction.tsv
#     {SP1}_vs_{SP2}_{DB}_{MODE}_comparison_{SP1}_only.tsv
#     {SP1}_vs_{SP2}_{DB}_{MODE}_comparison_{SP2}_only.tsv
#
################################################################################

#SBATCH --job-name=compare_species
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=30:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu

#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis
#SBATCH -o 06_analysis/compare_species_%j.out
#SBATCH -e 06_analysis/compare_species_%j.err

set -euo pipefail

source ~/.bashrc
conda activate rnaseq

# ============================================================================
# ARGUMENTS
# ============================================================================
SP1="${1:-}"
SP2="${2:-}"
DBNAME="${3:-}"
BLAST_MODE="${4:-discovery}"
PADJ="${5:-0.05}"
LFC="${6:-2.0}"

if [ -z "$SP1" ] || [ -z "$SP2" ] || [ -z "$DBNAME" ]; then
    echo "ERROR: Missing required arguments"
    echo ""
    echo "Usage: sbatch scripts/run_compare_species.sbatch SP1 SP2 DATABASE [BLAST_MODE] [PADJ] [LFC]"
    echo ""
    echo "  SP1:        First species  (DC, DG, or MF)"
    echo "  SP2:        Second species (DC, DG, or MF)"
    echo "  DATABASE:   swissprot, nr, refseq_protein, or DCdb"
    echo "  BLAST_MODE: discovery (default) or strict"
    echo "  PADJ:       padj cutoff (default: 0.05)"
    echo "  LFC:        |log2FC| cutoff (default: 2.0)"
    echo ""
    echo "Examples:"
    echo "  sbatch scripts/run_compare_species.sbatch DC DG swissprot discovery"
    echo "  sbatch scripts/run_compare_species.sbatch DC DG nr discovery 0.01 2.5"
    exit 1
fi

SP1=$(echo "$SP1" | tr '[:lower:]' '[:upper:]')
SP2=$(echo "$SP2" | tr '[:lower:]' '[:upper:]')

for SP in "$SP1" "$SP2"; do
    if [[ ! "$SP" =~ ^(DC|DG|MF)$ ]]; then
        echo "ERROR: Invalid species '$SP'. Must be DC, DG, or MF"
        exit 1
    fi
done

if [ "$SP1" = "$SP2" ]; then
    echo "ERROR: Species 1 and 2 must be different (got '$SP1' for both)"
    exit 1
fi

# ============================================================================
# PATHS
# ============================================================================
BASE_DIR="/projects/tholl_lab_1/daisy_analysis"
CODE_DIR="${BASE_DIR}/05_rnaseq-code"

SP1_FILE="${BASE_DIR}/06_analysis/combined_${SP1}/${SP1}_${DBNAME}_${BLAST_MODE}_annotated.tsv"
SP2_FILE="${BASE_DIR}/06_analysis/combined_${SP2}/${SP2}_${DBNAME}_${BLAST_MODE}_annotated.tsv"

OUTPUT_DIR="${BASE_DIR}/06_analysis/comparison_${SP1}_vs_${SP2}"
OUTPUT_FILE="${OUTPUT_DIR}/${SP1}_vs_${SP2}_${DBNAME}_${BLAST_MODE}_comparison.tsv"

mkdir -p "$OUTPUT_DIR"

# ============================================================================
# VALIDATE
# ============================================================================
if [ ! -f "$SP1_FILE" ]; then
    echo "ERROR: ${SP1} annotated results not found: $SP1_FILE"
    echo ""
    echo "Run combine+filter first:"
    echo "  sbatch scripts/run_combine_filter.sbatch $SP1 $DBNAME $BLAST_MODE"
    exit 1
fi

if [ ! -f "$SP2_FILE" ]; then
    echo "ERROR: ${SP2} annotated results not found: $SP2_FILE"
    echo ""
    echo "Run combine+filter first:"
    echo "  sbatch scripts/run_combine_filter.sbatch $SP2 $DBNAME $BLAST_MODE"
    exit 1
fi

echo "=========================================="
echo "Cross-Species Comparison: ${SP1} vs ${SP2}"
echo "=========================================="
echo ""
echo "Species 1:    $SP1  ($SP1_FILE)"
echo "Species 2:    $SP2  ($SP2_FILE)"
echo "Database:     $DBNAME ($BLAST_MODE)"
echo "Thresholds:   padj < $PADJ, |log2FC| > $LFC"
echo "Output:       $OUTPUT_DIR"
echo ""
date
echo "Hostname: $(hostname)"
echo "------------------------------------------"

# ============================================================================
# RUN COMPARISON
# ============================================================================
python3 "${CODE_DIR}/scripts/compare_species.py" \
    --sp1 "$SP1_FILE" \
    --sp2 "$SP2_FILE" \
    --sp1-name "$SP1" \
    --sp2-name "$SP2" \
    --output "$OUTPUT_FILE" \
    --padj "$PADJ" \
    --lfc "$LFC"

if [ $? -ne 0 ]; then
    echo "ERROR: Comparison failed"
    exit 1
fi

echo ""
echo "------------------------------------------"
echo "NEXT STEPS"
echo "------------------------------------------"
echo ""
echo "1. View the full comparison:"
echo "   less $OUTPUT_FILE"
echo ""
echo "2. See only shared DE genes:"
echo "   grep 'shared' $OUTPUT_FILE | head -20"
echo ""
echo "3. See ${SP1}-specific genes:"
echo "   grep '${SP1}_only' $OUTPUT_FILE | head -20"
echo ""
echo "4. Search for gene families:"
echo "   grep -i 'terpene' $OUTPUT_FILE"
echo "   grep -i 'cytochrome' $OUTPUT_FILE"
echo ""
date
echo "=========================================="
