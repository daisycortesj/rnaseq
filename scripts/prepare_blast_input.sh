#!/bin/bash
#
# prepare_blast_input.sh - Master script for BLAST input preparation
#
# This script runs the complete workflow:
# 1. Parse GTF for gene annotations
# 2. Inner join with gene IDs list
# 3. Extract CDS sequences for FASTA
#
# USAGE:
#   bash scripts/prepare_blast_input.sh DC
#   bash scripts/prepare_blast_input.sh DG
#   bash scripts/prepare_blast_input.sh MF
#
# Or with explicit paths:
#   bash scripts/prepare_blast_input.sh \
#     path/to/annotation.gtf \
#     path/to/genome.fna \
#     path/to/gene_ids.txt \
#     path/to/output_dir/
#
# Author: Daisy Cortes

set -e  # Exit on any error
set -u  # Exit on undefined variable
set -o pipefail  # Exit on pipe failure

# Color codes for pretty output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored messages
print_header() {
    echo -e "${BLUE}======================================================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}======================================================================${NC}"
}

print_step() {
    echo -e "${GREEN}✓ $1${NC}"
}

print_error() {
    echo -e "${RED}✗ ERROR: $1${NC}" >&2
}

print_warning() {
    echo -e "${YELLOW}⚠ WARNING: $1${NC}"
}

# Species configurations
declare -A SPECIES_GTF
declare -A SPECIES_GENOME
declare -A SPECIES_GENE_IDS
declare -A SPECIES_OUTPUT

SPECIES_GTF["DC"]="04_reference/dc_genomic.gtf"
SPECIES_GENOME["DC"]="04_reference/GCF_001625215.2_DH1_v3.0_genomic.fna"
SPECIES_GENE_IDS["DC"]="06_analysis/pydeseq2_DC_step1_unfiltered/all_gene_ids.txt"
SPECIES_OUTPUT["DC"]="06_analysis/blast_input_DC"

SPECIES_GTF["DG"]="04_reference/dc_genomic.gtf"
SPECIES_GENOME["DG"]="04_reference/GCF_001625215.2_DH1_v3.0_genomic.fna"
SPECIES_GENE_IDS["DG"]="06_analysis/pydeseq2_DG_step1_unfiltered/all_gene_ids.txt"
SPECIES_OUTPUT["DG"]="06_analysis/blast_input_DG"

SPECIES_GTF["MF"]="04_reference/mf_genomic.gtf"
SPECIES_GENOME["MF"]="04_reference/MYU_GWHGDHP00000000.1_genomic.fna"
SPECIES_GENE_IDS["MF"]="06_analysis/pydeseq2_MF_step1_unfiltered/all_gene_ids.txt"
SPECIES_OUTPUT["MF"]="06_analysis/blast_input_MF"

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Parse arguments
if [ $# -eq 0 ]; then
    echo "Usage: $0 <species_code>"
    echo "   or: $0 <gtf> <genome> <gene_ids> <output_dir>"
    echo ""
    echo "Species codes: DC, DG, MF"
    echo ""
    echo "Examples:"
    echo "  bash scripts/prepare_blast_input.sh DC"
    echo "  bash scripts/prepare_blast_input.sh DG"
    echo "  bash scripts/prepare_blast_input.sh path/to/gtf path/to/genome path/to/ids path/to/output"
    exit 1
fi

# Mode 1: Species code
if [ $# -eq 1 ]; then
    SPECIES=$(echo "$1" | tr '[:lower:]' '[:upper:]')
    
    if [[ ! -v SPECIES_GTF["$SPECIES"] ]]; then
        print_error "Unknown species code: $SPECIES"
        echo "Available species: DC, DG, MF"
        exit 1
    fi
    
    GTF_FILE="${PROJECT_ROOT}/${SPECIES_GTF[$SPECIES]}"
    GENOME_FILE="${PROJECT_ROOT}/${SPECIES_GENOME[$SPECIES]}"
    GENE_IDS_FILE="${PROJECT_ROOT}/${SPECIES_GENE_IDS[$SPECIES]}"
    OUTPUT_DIR="${PROJECT_ROOT}/${SPECIES_OUTPUT[$SPECIES]}"
    
    print_header "BLAST INPUT PREPARATION - Species: $SPECIES"

# Mode 2: Explicit paths
elif [ $# -eq 4 ]; then
    GTF_FILE="$1"
    GENOME_FILE="$2"
    GENE_IDS_FILE="$3"
    OUTPUT_DIR="$4"
    
    print_header "BLAST INPUT PREPARATION - Custom Paths"
else
    print_error "Invalid number of arguments"
    echo "Provide either: 1 argument (species code) OR 4 arguments (all paths)"
    exit 1
fi

# Print configuration
echo "Configuration:"
echo "  GTF file:     $GTF_FILE"
echo "  Genome file:  $GENOME_FILE"
echo "  Gene IDs:     $GENE_IDS_FILE"
echo "  Output dir:   $OUTPUT_DIR"
echo ""

# Check input files exist
for file in "$GTF_FILE" "$GENOME_FILE" "$GENE_IDS_FILE"; do
    if [ ! -f "$file" ]; then
        print_error "File not found: $file"
        exit 1
    fi
done

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Define output files
GENES_PARSED="$OUTPUT_DIR/genes_parsed.tsv"
FILTERED_GENES="$OUTPUT_DIR/filtered_genes_with_annotations.tsv"
FASTA_OUTPUT="$OUTPUT_DIR/all_genes_cds.fasta"

echo ""
print_header "STEP 1/3: Parse GTF for Gene Annotations"
echo "Running: parse_refseq_gtf.py"
echo ""

python "$SCRIPT_DIR/parse_refseq_gtf.py" \
    --gtf "$GTF_FILE" \
    --feature gene \
    --out "$GENES_PARSED"

if [ $? -eq 0 ]; then
    print_step "Step 1 completed successfully"
else
    print_error "Step 1 failed"
    exit 1
fi

echo ""
print_header "STEP 2/3: Filter Genes (Inner Join)"
echo "Running: join_gtf_with_gene_ids.py"
echo ""

python "$SCRIPT_DIR/join_gtf_with_gene_ids.py" \
    --parsed "$GENES_PARSED" \
    --gene-ids "$GENE_IDS_FILE" \
    --out "$FILTERED_GENES"

if [ $? -eq 0 ]; then
    print_step "Step 2 completed successfully"
else
    print_error "Step 2 failed"
    exit 1
fi

echo ""
print_header "STEP 3/3: Extract CDS Sequences"
echo "Running: extract_cds_from_gtf.py"
echo ""

python "$SCRIPT_DIR/extract_cds_from_gtf.py" \
    "$GTF_FILE" \
    "$GENOME_FILE" \
    "$GENE_IDS_FILE" \
    "$FASTA_OUTPUT"

if [ $? -eq 0 ]; then
    print_step "Step 3 completed successfully"
else
    print_error "Step 3 failed"
    exit 1
fi

# Print summary
echo ""
print_header "WORKFLOW COMPLETED SUCCESSFULLY!"
echo ""
echo "Output files created:"
echo "  1. $GENES_PARSED"
echo "     → All genes from GTF with annotations"
echo ""
echo "  2. $FILTERED_GENES"
echo "     → Only your genes with annotations (INNER JOIN)"
echo ""
echo "  3. $FASTA_OUTPUT"
echo "     → CDS sequences for BLAST"
echo ""
echo "Next steps:"
echo "  • Run BLAST on: $FASTA_OUTPUT"
echo "  • Use annotations from: $FILTERED_GENES"
echo ""
print_header "All Done! 🎉"
