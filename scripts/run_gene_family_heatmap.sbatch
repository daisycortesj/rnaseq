#!/bin/bash
################################################################################
# Extract verified CYP + OMT gene families and generate heatmaps
################################################################################
#
# Two-step pipeline:
#   Step 1: Extract CYP and OMT genes from BLAST (+ HMMER if available)
#   Step 2: Generate heatmaps from verified gene lists + count matrix
#
# HMMER is optional. Run with BLAST only first, then re-run when HMMER
# finishes to upgrade confidence scores.
#
# Requires:
#   - Combined annotated file (from run_combine_filter.sbatch)
#   - Count matrix (from build_count_matrix.sbatch)
#
# How to run:
#   # BLAST only:
#   sbatch scripts/run_gene_family_heatmap.sbatch DC swissprot discovery
#
#   # BLAST + HMMER:
#   sbatch scripts/run_gene_family_heatmap.sbatch DC swissprot discovery full
#
#   # For DG:
#   sbatch scripts/run_gene_family_heatmap.sbatch DG swissprot discovery
#
# Output:
#   06_analysis/gene_families_{SPECIES}/
#     {SP}_CYP_verified.tsv
#     {SP}_OMT_verified.tsv
#     {SP}_CYP_OMT_combined.tsv
#     {SP}_gene_family_summary.txt
#     cyp_heatmap.pdf / .png / .tsv
#     omt_heatmap.pdf / .png / .tsv
#     cyp_omt_combined_heatmap.pdf / .png / .tsv
#
################################################################################

#SBATCH --job-name=gene_family
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=30:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu

#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis
#SBATCH -o 06_analysis/gene_family_%j.out
#SBATCH -e 06_analysis/gene_family_%j.err

set -euo pipefail

source ~/.bashrc
conda activate rnaseq

# ============================================================================
# ARGUMENTS
# ============================================================================
SPECIES="${1:-}"
DBNAME="${2:-swissprot}"
BLAST_MODE="${3:-discovery}"
EVIDENCE="${4:-blast}"

if [ -z "$SPECIES" ]; then
    echo "ERROR: Missing species argument"
    echo ""
    echo "Usage: sbatch scripts/run_gene_family_heatmap.sbatch SPECIES [DATABASE] [BLAST_MODE] [EVIDENCE]"
    echo ""
    echo "  SPECIES:   DC, DG, or MF"
    echo "  DATABASE:  swissprot (default), nr, refseq_protein"
    echo "  BLAST_MODE: discovery (default) or strict"
    echo "  EVIDENCE:  blast (default) or full (includes HMMER)"
    echo ""
    echo "Examples:"
    echo "  sbatch scripts/run_gene_family_heatmap.sbatch DC swissprot discovery"
    echo "  sbatch scripts/run_gene_family_heatmap.sbatch DC swissprot discovery full"
    exit 1
fi

SPECIES=$(echo "$SPECIES" | tr '[:lower:]' '[:upper:]')

if [[ ! "$SPECIES" =~ ^(DC|DG|MF)$ ]]; then
    echo "ERROR: Invalid species '$SPECIES'. Must be DC, DG, or MF"
    exit 1
fi

# ============================================================================
# PATHS
# ============================================================================
BASE_DIR="/projects/tholl_lab_1/daisy_analysis"
CODE_DIR="${BASE_DIR}/05_rnaseq-code"

BLAST_FILE="${BASE_DIR}/06_analysis/combined_${SPECIES}/${SPECIES}_${DBNAME}_${BLAST_MODE}_annotated.tsv"
OUTPUT_DIR="${BASE_DIR}/06_analysis/gene_families_${SPECIES}"

# Count matrix: try step1 location first, then step3 location
case "$SPECIES" in
    DC) COUNT_FOLDER="00_1_DC" ;;
    DG) COUNT_FOLDER="00_2_DG" ;;
    MF) COUNT_FOLDER="00_3_MF" ;;
esac

COUNT_MATRIX="${BASE_DIR}/03_count_tables/${COUNT_FOLDER}/gene_count_matrix.tsv"
METADATA="${BASE_DIR}/03_count_tables/${COUNT_FOLDER}/sample_metadata.tsv"

# HMMER (optional)
HMMER_FILE="${BASE_DIR}/06_analysis/hmmer_${SPECIES}/all_genes_protein_pfam_domains.txt"

# ============================================================================
# VALIDATE
# ============================================================================
if [ ! -f "$BLAST_FILE" ]; then
    echo "ERROR: Combined annotated file not found: $BLAST_FILE"
    echo ""
    echo "Run combine+filter first:"
    echo "  sbatch scripts/run_combine_filter.sbatch $SPECIES $DBNAME $BLAST_MODE"
    exit 1
fi

if [ ! -f "$COUNT_MATRIX" ]; then
    echo "ERROR: Count matrix not found: $COUNT_MATRIX"
    exit 1
fi

if [ ! -f "$METADATA" ]; then
    echo "ERROR: Metadata not found: $METADATA"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

echo "=========================================="
echo "Gene Family Extraction + Heatmap"
echo "=========================================="
echo ""
echo "Species:      $SPECIES"
echo "Database:     $DBNAME ($BLAST_MODE)"
echo "Evidence:     $EVIDENCE"
echo "BLAST file:   $BLAST_FILE"
echo "Count matrix: $COUNT_MATRIX"
echo "Output:       $OUTPUT_DIR"
echo ""

# Check HMMER availability
HMMER_ARG=""
if [ "$EVIDENCE" = "full" ]; then
    if [ -f "$HMMER_FILE" ]; then
        echo "HMMER:        $HMMER_FILE (found)"
        HMMER_ARG="--hmmer $HMMER_FILE"
    else
        echo "HMMER:        not found (skipping)"
    fi
else
    echo "HMMER:        skipped (use 'full' argument to include)"
fi

echo ""
date
echo "------------------------------------------"

# ============================================================================
# STEP 1: EXTRACT GENE FAMILIES
# ============================================================================
echo ""
echo "=========================================="
echo "STEP 1/2: Extracting CYP + OMT Genes"
echo "=========================================="
echo ""

python3 "${CODE_DIR}/scripts/extract_gene_families.py" \
    --blast "$BLAST_FILE" \
    $HMMER_ARG \
    --species "$SPECIES" \
    --output "$OUTPUT_DIR"

if [ $? -ne 0 ]; then
    echo "ERROR: Gene family extraction failed"
    exit 1
fi

# ============================================================================
# STEP 2: GENERATE HEATMAPS
# ============================================================================
echo ""
echo "=========================================="
echo "STEP 2/2: Generating Heatmaps"
echo "=========================================="
echo ""

# CYP heatmap
CYP_FILE="${OUTPUT_DIR}/${SPECIES}_CYP_verified.tsv"
if [ -f "$CYP_FILE" ]; then
    CYP_COUNT=$(tail -n +2 "$CYP_FILE" | wc -l)
    if [ "$CYP_COUNT" -gt 0 ]; then
        echo "Generating CYP heatmap ($CYP_COUNT genes)..."
        python3 "${CODE_DIR}/scripts/generate_family_heatmap.py" \
            --genes "$CYP_FILE" \
            --counts "$COUNT_MATRIX" \
            --metadata "$METADATA" \
            --output "${OUTPUT_DIR}/cyp_heatmap.pdf" \
            --title "${SPECIES} Cytochrome P450 (CYP) Expression"
    else
        echo "No CYP genes found, skipping heatmap"
    fi
fi

# OMT heatmap
OMT_FILE="${OUTPUT_DIR}/${SPECIES}_OMT_verified.tsv"
if [ -f "$OMT_FILE" ]; then
    OMT_COUNT=$(tail -n +2 "$OMT_FILE" | wc -l)
    if [ "$OMT_COUNT" -gt 0 ]; then
        echo ""
        echo "Generating OMT heatmap ($OMT_COUNT genes)..."
        python3 "${CODE_DIR}/scripts/generate_family_heatmap.py" \
            --genes "$OMT_FILE" \
            --counts "$COUNT_MATRIX" \
            --metadata "$METADATA" \
            --output "${OUTPUT_DIR}/omt_heatmap.pdf" \
            --title "${SPECIES} O-Methyltransferase (OMT) Expression"
    else
        echo "No OMT genes found, skipping heatmap"
    fi
fi

# Combined heatmap
COMBINED_FILE="${OUTPUT_DIR}/${SPECIES}_CYP_OMT_combined.tsv"
if [ -f "$COMBINED_FILE" ]; then
    COMBINED_COUNT=$(tail -n +2 "$COMBINED_FILE" | grep -cv "^#" || true)
    if [ "$COMBINED_COUNT" -gt 0 ]; then
        echo ""
        echo "Generating combined CYP+OMT heatmap ($COMBINED_COUNT genes)..."
        python3 "${CODE_DIR}/scripts/generate_family_heatmap.py" \
            --genes "$COMBINED_FILE" \
            --counts "$COUNT_MATRIX" \
            --metadata "$METADATA" \
            --output "${OUTPUT_DIR}/cyp_omt_combined_heatmap.pdf" \
            --title "${SPECIES} CYP + OMT Expression" \
            --color-by-family
    fi
fi

# Significant-only heatmaps
if [ -f "$CYP_FILE" ]; then
    echo ""
    echo "Generating CYP significant-only heatmap..."
    python3 "${CODE_DIR}/scripts/generate_family_heatmap.py" \
        --genes "$CYP_FILE" \
        --counts "$COUNT_MATRIX" \
        --metadata "$METADATA" \
        --output "${OUTPUT_DIR}/cyp_heatmap_sig_only.pdf" \
        --title "${SPECIES} CYP DE Genes Only" \
        --sig-only || echo "  (no significant CYP genes to plot)"
fi

# ============================================================================
# DONE
# ============================================================================
echo ""
echo "=========================================="
echo "GENE FAMILY ANALYSIS COMPLETE!"
echo "=========================================="
echo ""
echo "Output directory: $OUTPUT_DIR"
echo ""
echo "Files:"
ls -la "$OUTPUT_DIR"/ 2>/dev/null | grep -v "^total" | grep -v "^\." || true
echo ""
echo "------------------------------------------"
echo "NEXT STEPS"
echo "------------------------------------------"
echo ""
echo "1. View heatmaps:"
echo "   scp daisycortesj@tinkercliffs.arc.vt.edu:${OUTPUT_DIR}/*.pdf ~/Desktop/"
echo ""
echo "2. Re-run with full evidence (after HMMER finishes):"
echo "   sbatch scripts/run_gene_family_heatmap.sbatch $SPECIES $DBNAME $BLAST_MODE full"
echo ""
echo "3. Compare DC vs DG (after both species are done):"
echo "   sbatch scripts/run_compare_species.sbatch DC DG $DBNAME $BLAST_MODE"
echo ""
date
echo "=========================================="
