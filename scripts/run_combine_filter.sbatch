#!/bin/bash
################################################################################
# Combine BLAST annotations with DESeq2 statistics, then filter
################################################################################
#
# Two-step pipeline:
#   Step 1: Merge BLAST results with PyDESeq2 expression data
#           → One table per gene: expression stats + functional annotation
#   Step 2: Filter to significant, well-annotated candidates
#
# Requires:
#   - PyDESeq2 step 1 results (run_pydeseq2_step1_analysis.sbatch)
#   - BLAST results (blastp_discoveryfilter.sbatch or blastp_strictfilter.sbatch)
#
# How to run:
#   sbatch scripts/run_combine_filter.sbatch DC swissprot discovery
#   sbatch scripts/run_combine_filter.sbatch DC nr discovery
#   sbatch scripts/run_combine_filter.sbatch DC swissprot strict
#   sbatch scripts/run_combine_filter.sbatch DG swissprot discovery
#   sbatch scripts/run_combine_filter.sbatch MF nr discovery
#
# Filter modes (optional 4th argument, default: standard):
#   standard  - padj<0.05, |log2FC|>2.0 (recommended)
#   strict    - padj<0.01, |log2FC|>2.5 (publication-ready)
#   lenient   - padj<0.1,  |log2FC|>1.0 (exploratory)
#   blast_only - any DE, but must have BLAST hit
#
# Examples:
#   sbatch scripts/run_combine_filter.sbatch DC swissprot discovery strict
#   sbatch scripts/run_combine_filter.sbatch DC nr discovery lenient
#
# Output:
#   06_analysis/combined_{SPECIES}/
#     {SPECIES}_{DB}_{MODE}_annotated.tsv       (combined, all genes)
#     {SPECIES}_{DB}_{MODE}_filtered_{FILTER}.tsv (filtered candidates)
#     {SPECIES}_{DB}_{MODE}_filtered_{FILTER}_CYP_only.tsv (CYP genes)
#
################################################################################

#SBATCH --job-name=combine_filter
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=30:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu

#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis
#SBATCH -o 06_analysis/combine_filter_%j.out
#SBATCH -e 06_analysis/combine_filter_%j.err

set -euo pipefail

source ~/.bashrc
conda activate rnaseq

# ============================================================================
# ARGUMENTS
# ============================================================================
SPECIES="${1:-}"
DBNAME="${2:-}"
BLAST_MODE="${3:-discovery}"
FILTER_MODE="${4:-standard}"

if [ -z "$SPECIES" ] || [ -z "$DBNAME" ]; then
    echo "ERROR: Missing required arguments"
    echo ""
    echo "Usage: sbatch scripts/run_combine_filter.sbatch SPECIES DATABASE [BLAST_MODE] [FILTER_MODE]"
    echo ""
    echo "  SPECIES:    DC, DG, or MF"
    echo "  DATABASE:   swissprot, nr, refseq_protein, or DCdb"
    echo "  BLAST_MODE: discovery (default) or strict"
    echo "  FILTER_MODE: standard (default), strict, lenient, or blast_only"
    echo ""
    echo "Examples:"
    echo "  sbatch scripts/run_combine_filter.sbatch DC swissprot discovery"
    echo "  sbatch scripts/run_combine_filter.sbatch DC nr discovery strict"
    echo "  sbatch scripts/run_combine_filter.sbatch DG swissprot discovery lenient"
    exit 1
fi

SPECIES=$(echo "$SPECIES" | tr '[:lower:]' '[:upper:]')

if [[ ! "$SPECIES" =~ ^(DC|DG|MF)$ ]]; then
    echo "ERROR: Invalid species '$SPECIES'. Must be DC, DG, or MF"
    exit 1
fi

if [[ ! "$DBNAME" =~ ^(nr|swissprot|refseq_protein|DCdb)$ ]]; then
    echo "ERROR: Invalid database '$DBNAME'. Must be nr, swissprot, refseq_protein, or DCdb"
    exit 1
fi

if [[ ! "$BLAST_MODE" =~ ^(discovery|strict)$ ]]; then
    echo "ERROR: Invalid BLAST mode '$BLAST_MODE'. Must be discovery or strict"
    exit 1
fi

if [[ ! "$FILTER_MODE" =~ ^(standard|strict|lenient|blast_only)$ ]]; then
    echo "ERROR: Invalid filter mode '$FILTER_MODE'. Must be standard, strict, lenient, or blast_only"
    exit 1
fi

# ============================================================================
# PATHS
# ============================================================================
BASE_DIR="/projects/tholl_lab_1/daisy_analysis"
CODE_DIR="${BASE_DIR}/05_rnaseq-code"

DESEQ_FILE="${BASE_DIR}/06_analysis/pydeseq2_${SPECIES}_step1_unfiltered/pydeseq2_results_UNFILTERED.tsv"

# Determine BLAST results path based on blast type (blastp vs blastx)
# Default to blastp; if not found, try blastx
BLASTP_FILE="${BASE_DIR}/06_analysis/blastp_${SPECIES}/blastp_${SPECIES}_${DBNAME}_${BLAST_MODE}.tsv"
BLASTX_FILE="${BASE_DIR}/06_analysis/blastx_${SPECIES}/blastx_${SPECIES}_${DBNAME}_${BLAST_MODE}.tsv"

if [ -f "$BLASTP_FILE" ]; then
    BLAST_FILE="$BLASTP_FILE"
    BLAST_TYPE="blastp"
elif [ -f "$BLASTX_FILE" ]; then
    BLAST_FILE="$BLASTX_FILE"
    BLAST_TYPE="blastx"
else
    echo "ERROR: No BLAST results found. Looked for:"
    echo "  $BLASTP_FILE"
    echo "  $BLASTX_FILE"
    echo ""
    echo "Run BLAST first:"
    echo "  sbatch scripts/blastp_discoveryfilter.sbatch $SPECIES $DBNAME"
    exit 1
fi

OUTPUT_DIR="${BASE_DIR}/06_analysis/combined_${SPECIES}"
COMBINED_FILE="${OUTPUT_DIR}/${SPECIES}_${DBNAME}_${BLAST_MODE}_annotated.tsv"
FILTERED_FILE="${OUTPUT_DIR}/${SPECIES}_${DBNAME}_${BLAST_MODE}_filtered_${FILTER_MODE}.tsv"

mkdir -p "$OUTPUT_DIR"

# ============================================================================
# VALIDATE INPUTS
# ============================================================================
if [ ! -f "$DESEQ_FILE" ]; then
    echo "ERROR: PyDESeq2 results not found: $DESEQ_FILE"
    echo ""
    echo "Run PyDESeq2 step 1 first:"
    echo "  sbatch scripts/run_pydeseq2_step1_analysis.sbatch $SPECIES"
    exit 1
fi

echo "=========================================="
echo "Combine BLAST + DESeq2, then Filter"
echo "=========================================="
echo ""
echo "Species:      $SPECIES"
echo "Database:     $DBNAME"
echo "BLAST mode:   $BLAST_MODE ($BLAST_TYPE)"
echo "Filter mode:  $FILTER_MODE"
echo ""
echo "DESeq2 file:  $DESEQ_FILE"
echo "BLAST file:   $BLAST_FILE"
echo "Output dir:   $OUTPUT_DIR"
echo ""
date
echo "Hostname: $(hostname)"
echo "------------------------------------------"

# ============================================================================
# STEP 1: COMBINE BLAST + DESeq2
# ============================================================================
echo ""
echo "=========================================="
echo "STEP 1/2: Combining BLAST + DESeq2"
echo "=========================================="
echo ""

python3 "${CODE_DIR}/scripts/combine_blast_deseq.py" \
    --deseq "$DESEQ_FILE" \
    --blast "$BLAST_FILE" \
    --output "$COMBINED_FILE"

if [ $? -ne 0 ]; then
    echo "ERROR: Combine step failed"
    exit 1
fi

echo "✓ Combined file: $COMBINED_FILE"

# ============================================================================
# STEP 2: FILTER COMBINED RESULTS
# ============================================================================
echo ""
echo "=========================================="
echo "STEP 2/2: Filtering ($FILTER_MODE mode)"
echo "=========================================="
echo ""

python3 "${CODE_DIR}/scripts/filter_combined_results.py" \
    --input "$COMBINED_FILE" \
    --output "$FILTERED_FILE" \
    --mode "$FILTER_MODE"

if [ $? -ne 0 ]; then
    echo "ERROR: Filter step failed"
    exit 1
fi

echo "✓ Filtered file: $FILTERED_FILE"

# ============================================================================
# SUMMARY
# ============================================================================
echo ""
echo "=========================================="
echo "COMBINE + FILTER COMPLETE!"
echo "=========================================="
echo ""
echo "Output files in: $OUTPUT_DIR"
echo ""

if [ -f "$COMBINED_FILE" ]; then
    TOTAL=$(tail -n +2 "$COMBINED_FILE" | grep -cv "^#" || true)
    echo "  All genes (annotated):  $COMBINED_FILE"
    echo "    → $TOTAL genes total"
fi

if [ -f "$FILTERED_FILE" ]; then
    FILTERED=$(tail -n +2 "$FILTERED_FILE" | grep -cv "^#" || true)
    echo "  Filtered candidates:    $FILTERED_FILE"
    echo "    → $FILTERED genes pass $FILTER_MODE filters"
fi

CYP_FILE="${OUTPUT_DIR}/${SPECIES}_${DBNAME}_${BLAST_MODE}_annotated_CYP_only.tsv"
if [ -f "$CYP_FILE" ]; then
    CYP_COUNT=$(tail -n +2 "$CYP_FILE" | wc -l)
    echo "  CYP genes:              $CYP_FILE"
    echo "    → $CYP_COUNT CYP genes found"
fi

echo ""
echo "------------------------------------------"
echo "NEXT STEPS"
echo "------------------------------------------"
echo ""
echo "1. View your top candidates:"
echo "   head -20 $FILTERED_FILE | column -t -s\$'\\t'"
echo ""
echo "2. Search for specific gene families:"
echo "   grep -i 'terpene' $COMBINED_FILE"
echo "   grep -i 'cytochrome' $COMBINED_FILE"
echo ""
echo "3. Try different filter stringency:"
echo "   sbatch scripts/run_combine_filter.sbatch $SPECIES $DBNAME $BLAST_MODE strict"
echo "   sbatch scripts/run_combine_filter.sbatch $SPECIES $DBNAME $BLAST_MODE lenient"
echo ""
echo "4. Combine with a different BLAST database:"
echo "   sbatch scripts/run_combine_filter.sbatch $SPECIES nr $BLAST_MODE"
echo ""
date
echo "=========================================="
