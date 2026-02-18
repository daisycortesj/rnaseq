#!/bin/bash
################################################################################
# Combined gene family extraction + heatmaps for TWO species
################################################################################
#
# Extracts CYP + OMT genes from both species in one run, using BLAST + HMMER
# evidence from each species. Produces a merged table and heatmaps showing
# both species side by side.
#
# Both species must share the same reference genome (same gene IDs).
# HMMER is optional -- uses whatever evidence is available.
#
# Requires (for each species):
#   - Combined annotated file (from run_combine_filter.sbatch)
#   - Count matrix (from build_count_matrix.sbatch)
#   - Optionally: HMMER results (from run_hmmer.sbatch)
#
# How to run:
#   # BLAST only:
#   sbatch scripts/run_gene_family_combined.sbatch DC DG swissprot discovery
#
#   # Full evidence (HMMER included):
#   sbatch scripts/run_gene_family_combined.sbatch DC DG swissprot discovery full
#
# Output:
#   06_analysis/gene_families_DC_DG/
#     DC_DG_CYP_verified.tsv
#     DC_DG_OMT_verified.tsv
#     DC_DG_CYP_OMT_combined.tsv
#     DC_DG_CYP_OMT_shared_same_direction.tsv
#     DC_DG_CYP_OMT_DC_only.tsv
#     DC_DG_CYP_OMT_DG_only.tsv
#     DC_DG_gene_family_summary.txt
#     cyp_heatmap_DC.pdf / cyp_heatmap_DG.pdf
#     omt_heatmap_DC.pdf / omt_heatmap_DG.pdf
#
################################################################################

#SBATCH --job-name=gene_fam_combined
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=30:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu

#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis
#SBATCH -o 06_analysis/gene_family_combined_%j.out
#SBATCH -e 06_analysis/gene_family_combined_%j.err

set -euo pipefail

source ~/.bashrc
conda activate rnaseq

# ============================================================================
# ARGUMENTS
# ============================================================================
SP1="${1:-}"
SP2="${2:-}"
DBNAME="${3:-swissprot}"
BLAST_MODE="${4:-discovery}"
EVIDENCE="${5:-blast}"

if [ -z "$SP1" ] || [ -z "$SP2" ]; then
    echo "ERROR: Missing species arguments"
    echo ""
    echo "Usage: sbatch scripts/run_gene_family_combined.sbatch SP1 SP2 [DATABASE] [BLAST_MODE] [EVIDENCE]"
    echo ""
    echo "  SP1:        First species  (DC, DG, or MF)"
    echo "  SP2:        Second species (DC, DG, or MF)"
    echo "  DATABASE:   swissprot (default), nr, refseq_protein"
    echo "  BLAST_MODE: discovery (default) or strict"
    echo "  EVIDENCE:   blast (default) or full (includes HMMER)"
    echo ""
    echo "Examples:"
    echo "  sbatch scripts/run_gene_family_combined.sbatch DC DG swissprot discovery"
    echo "  sbatch scripts/run_gene_family_combined.sbatch DC DG swissprot discovery full"
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
    echo "ERROR: Species must be different. Use run_gene_family_heatmap.sbatch for single species."
    exit 1
fi

# ============================================================================
# PATHS
# ============================================================================
BASE_DIR="/projects/tholl_lab_1/daisy_analysis"
CODE_DIR="${BASE_DIR}/05_rnaseq-code"

SP1_BLAST="${BASE_DIR}/06_analysis/combined_${SP1}/${SP1}_${DBNAME}_${BLAST_MODE}_annotated.tsv"
SP2_BLAST="${BASE_DIR}/06_analysis/combined_${SP2}/${SP2}_${DBNAME}_${BLAST_MODE}_annotated.tsv"

OUTPUT_DIR="${BASE_DIR}/06_analysis/gene_families_${SP1}_${SP2}"

# Count matrices
get_count_folder() {
    case "$1" in
        DC) echo "00_1_DC" ;;
        DG) echo "00_2_DG" ;;
        MF) echo "00_3_MF" ;;
    esac
}

SP1_COUNT_FOLDER=$(get_count_folder "$SP1")
SP2_COUNT_FOLDER=$(get_count_folder "$SP2")

SP1_COUNTS="${BASE_DIR}/03_count_tables/${SP1_COUNT_FOLDER}/gene_count_matrix.tsv"
SP2_COUNTS="${BASE_DIR}/03_count_tables/${SP2_COUNT_FOLDER}/gene_count_matrix.tsv"
SP1_META="${BASE_DIR}/03_count_tables/${SP1_COUNT_FOLDER}/sample_metadata.tsv"
SP2_META="${BASE_DIR}/03_count_tables/${SP2_COUNT_FOLDER}/sample_metadata.tsv"

# HMMER (optional)
SP1_HMMER="${BASE_DIR}/06_analysis/hmmer_${SP1}/all_genes_protein_pfam_domains.txt"
SP2_HMMER="${BASE_DIR}/06_analysis/hmmer_${SP2}/all_genes_protein_pfam_domains.txt"

# ============================================================================
# VALIDATE
# ============================================================================
for F in "$SP1_BLAST" "$SP2_BLAST"; do
    if [ ! -f "$F" ]; then
        echo "ERROR: Annotated file not found: $F"
        echo "Run combine+filter first for both species."
        exit 1
    fi
done

mkdir -p "$OUTPUT_DIR"

echo "=========================================="
echo "Combined Gene Family: ${SP1} + ${SP2}"
echo "=========================================="
echo ""
echo "Species 1:    $SP1  ($SP1_BLAST)"
echo "Species 2:    $SP2  ($SP2_BLAST)"
echo "Database:     $DBNAME ($BLAST_MODE)"
echo "Evidence:     $EVIDENCE"
echo "Output:       $OUTPUT_DIR"
echo ""

# Build HMMER arguments
SP1_HMMER_ARG=""
SP2_HMMER_ARG=""

if [ "$EVIDENCE" = "full" ]; then
    [ -f "$SP1_HMMER" ]   && SP1_HMMER_ARG="--sp1-hmmer $SP1_HMMER"     && echo "${SP1} HMMER:   found"
    [ -f "$SP2_HMMER" ]   && SP2_HMMER_ARG="--sp2-hmmer $SP2_HMMER"     && echo "${SP2} HMMER:   found"
    [ ! -f "$SP1_HMMER" ] && echo "${SP1} HMMER:   not found (skipping)"
    [ ! -f "$SP2_HMMER" ] && echo "${SP2} HMMER:   not found (skipping)"
else
    echo "HMMER: skipped (use 'full' to include)"
fi

echo ""
date
echo "------------------------------------------"

# ============================================================================
# STEP 1: EXTRACT GENE FAMILIES (COMBINED)
# ============================================================================
echo ""
echo "=========================================="
echo "STEP 1/2: Extracting CYP + OMT from ${SP1} + ${SP2}"
echo "=========================================="
echo ""

python3 "${CODE_DIR}/scripts/extract_gene_families_combined.py" \
    --sp1-blast "$SP1_BLAST" \
    --sp2-blast "$SP2_BLAST" \
    $SP1_HMMER_ARG $SP2_HMMER_ARG \
    --sp1-name "$SP1" \
    --sp2-name "$SP2" \
    --output "$OUTPUT_DIR"

if [ $? -ne 0 ]; then
    echo "ERROR: Combined extraction failed"
    exit 1
fi

# ============================================================================
# STEP 2: GENERATE HEATMAPS (per species, using verified gene lists)
# ============================================================================
echo ""
echo "=========================================="
echo "STEP 2/2: Generating Heatmaps"
echo "=========================================="
echo ""

# CYP heatmap for each species
CYP_FILE="${OUTPUT_DIR}/${SP1}_${SP2}_CYP_verified.tsv"
if [ -f "$CYP_FILE" ]; then
    CYP_COUNT=$(tail -n +2 "$CYP_FILE" | wc -l)
    if [ "$CYP_COUNT" -gt 0 ]; then
        if [ -f "$SP1_COUNTS" ] && [ -f "$SP1_META" ]; then
            echo "Generating ${SP1} CYP heatmap..."
            python3 "${CODE_DIR}/scripts/generate_family_heatmap.py" \
                --genes "$CYP_FILE" \
                --counts "$SP1_COUNTS" \
                --metadata "$SP1_META" \
                --output "${OUTPUT_DIR}/cyp_heatmap_${SP1}.pdf" \
                --title "${SP1} Cytochrome P450 (CYP) Expression" || true
        fi

        if [ -f "$SP2_COUNTS" ] && [ -f "$SP2_META" ]; then
            echo "Generating ${SP2} CYP heatmap..."
            python3 "${CODE_DIR}/scripts/generate_family_heatmap.py" \
                --genes "$CYP_FILE" \
                --counts "$SP2_COUNTS" \
                --metadata "$SP2_META" \
                --output "${OUTPUT_DIR}/cyp_heatmap_${SP2}.pdf" \
                --title "${SP2} Cytochrome P450 (CYP) Expression" || true
        fi
    fi
fi

# OMT heatmap for each species
OMT_FILE="${OUTPUT_DIR}/${SP1}_${SP2}_OMT_verified.tsv"
if [ -f "$OMT_FILE" ]; then
    OMT_COUNT=$(tail -n +2 "$OMT_FILE" | wc -l)
    if [ "$OMT_COUNT" -gt 0 ]; then
        if [ -f "$SP1_COUNTS" ] && [ -f "$SP1_META" ]; then
            echo "Generating ${SP1} OMT heatmap..."
            python3 "${CODE_DIR}/scripts/generate_family_heatmap.py" \
                --genes "$OMT_FILE" \
                --counts "$SP1_COUNTS" \
                --metadata "$SP1_META" \
                --output "${OUTPUT_DIR}/omt_heatmap_${SP1}.pdf" \
                --title "${SP1} O-Methyltransferase (OMT) Expression" || true
        fi

        if [ -f "$SP2_COUNTS" ] && [ -f "$SP2_META" ]; then
            echo "Generating ${SP2} OMT heatmap..."
            python3 "${CODE_DIR}/scripts/generate_family_heatmap.py" \
                --genes "$OMT_FILE" \
                --counts "$SP2_COUNTS" \
                --metadata "$SP2_META" \
                --output "${OUTPUT_DIR}/omt_heatmap_${SP2}.pdf" \
                --title "${SP2} O-Methyltransferase (OMT) Expression" || true
        fi
    fi
fi

# ============================================================================
# DONE
# ============================================================================
echo ""
echo "=========================================="
echo "COMBINED GENE FAMILY ANALYSIS COMPLETE!"
echo "=========================================="
echo ""
echo "Output: $OUTPUT_DIR"
echo ""
echo "Files:"
ls -la "$OUTPUT_DIR"/ 2>/dev/null | grep -v "^total" | grep -v "^\." || true
echo ""
echo "------------------------------------------"
echo "NEXT STEPS"
echo "------------------------------------------"
echo ""
echo "1. Download results:"
echo "   scp -r daisycortesj@tinkercliffs.arc.vt.edu:${OUTPUT_DIR}/ ~/Desktop/gene_families_${SP1}_${SP2}/"
echo ""
echo "2. Re-run with full evidence (if HMMER finishes later):"
echo "   sbatch scripts/run_gene_family_combined.sbatch $SP1 $SP2 $DBNAME $BLAST_MODE full"
echo ""
date
echo "=========================================="
