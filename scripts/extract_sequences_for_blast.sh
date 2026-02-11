#!/bin/bash
#
# Extract protein sequences for ALL genes from PyDESeq2 results
# Requires a protein FASTA file (from NCBI).
#
# Usage:
#   bash scripts/extract_sequences_for_blast.sh DC    (D. carota)
#   bash scripts/extract_sequences_for_blast.sh DG    (Daucus glaber)
#   bash scripts/extract_sequences_for_blast.sh MF    (Myristica fragrans)
#
# Override paths with env vars:
#   PYDESEQ_RESULTS=... PROTEIN_FASTA=... OUTPUT_DIR=... bash scripts/extract_sequences_for_blast.sh DC

set -e

BASE_DIR="${BASE_DIR:-.}"

# Get species from command line
SPECIES="${1:-}"

if [ -n "${SPECIES}" ]; then
    SPECIES=$(echo "${SPECIES}" | tr '[:lower:]' '[:upper:]')
    case "${SPECIES}" in
        DC)
            PYDESEQ_RESULTS="${PYDESEQ_RESULTS:-${BASE_DIR}/06_analysis/pydeseq2_DC_step1_unfiltered/pydeseq2_results_UNFILTERED.tsv}"
            PROTEIN_FASTA="${PROTEIN_FASTA:-${BASE_DIR}/04_reference/GCF_001625215.2_DH1_v3.0_protein.faa}"
            OUTPUT_DIR="${OUTPUT_DIR:-${BASE_DIR}/06_analysis/blast_input_DC}"
            SPECIES_NAME="D. carota ssp. maximus"
            ;;
        DG)
            PYDESEQ_RESULTS="${PYDESEQ_RESULTS:-${BASE_DIR}/06_analysis/pydeseq2_DG_step1_unfiltered/pydeseq2_results_UNFILTERED.tsv}"
            PROTEIN_FASTA="${PROTEIN_FASTA:-${BASE_DIR}/04_reference/GCF_001625215.2_DH1_v3.0_protein.faa}"
            OUTPUT_DIR="${OUTPUT_DIR:-${BASE_DIR}/06_analysis/blast_input_DG}"
            SPECIES_NAME="Daucus glaber"
            ;;
        MF)
            PYDESEQ_RESULTS="${PYDESEQ_RESULTS:-${BASE_DIR}/06_analysis/pydeseq2_MF_step1_unfiltered/pydeseq2_results_UNFILTERED.tsv}"
            PROTEIN_FASTA="${PROTEIN_FASTA:-${BASE_DIR}/04_reference/MYU_GWHGDHP00000000.1_protein.faa}"
            OUTPUT_DIR="${OUTPUT_DIR:-${BASE_DIR}/06_analysis/blast_input_MF}"
            SPECIES_NAME="Myristica fragrans"
            ;;
        *)
            echo "ERROR: Unknown species: ${SPECIES}"
            echo "Valid options: DC, DG, MF"
            echo ""
            echo "Usage: bash scripts/extract_sequences_for_blast.sh [DC|DG|MF]"
            exit 1
            ;;
    esac
else
    # Fallback to env vars / defaults
    PYDESEQ_RESULTS="${PYDESEQ_RESULTS:-${BASE_DIR}/06_analysis/pydeseq2_step1_unfiltered/pydeseq2_results_UNFILTERED.tsv}"
    PROTEIN_FASTA="${PROTEIN_FASTA:-${BASE_DIR}/04_reference/GCF_001625215.2_DH1_v3.0_protein.faa}"
    OUTPUT_DIR="${OUTPUT_DIR:-${BASE_DIR}/06_analysis/blast_input}"
    SPECIES_NAME="(default)"
fi

echo "=========================================="
echo "Extract Sequences for BLAST"
echo "=========================================="
echo ""
echo "Species: ${SPECIES_NAME}"
echo ""

# Check inputs
if [ ! -f "$PYDESEQ_RESULTS" ]; then
    echo "ERROR: PyDESeq2 results not found: $PYDESEQ_RESULTS"
    echo "Run PyDESeq2 Step 1 first!"
    exit 1
fi

if [ ! -f "$PROTEIN_FASTA" ]; then
    echo "ERROR: Protein FASTA not found: $PROTEIN_FASTA"
    echo ""
    echo "Options:"
    echo "  1. Download from NCBI and set PROTEIN_FASTA=/path/to/protein.faa"
    echo "  2. Use Option B from BLAST_ALL_GENES_WORKFLOW.md to extract from GTF"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "Input files:"
echo "  PyDESeq2 results: $PYDESEQ_RESULTS"
echo "  Protein FASTA:    $PROTEIN_FASTA"
echo "  Output directory: $OUTPUT_DIR"
echo ""

# Step 1: Extract gene IDs (exclude NA genes)
echo "Step 1: Extracting gene IDs from PyDESeq2 results..."
GENE_IDS="$OUTPUT_DIR/all_gene_ids.txt"

awk -F'\t' 'NR>1 && $6!="NA" {print $1}' "$PYDESEQ_RESULTS" > "$GENE_IDS"

TOTAL_GENES=$(wc -l < "$GENE_IDS")
echo "  ✓ Extracted $TOTAL_GENES gene IDs"
echo "  ✓ Saved to: $GENE_IDS"
echo ""

# Step 2: Extract protein sequences
echo "Step 2: Extracting protein sequences..."
SEQUENCES="$OUTPUT_DIR/all_genes.faa"

grep -A1 -Ff "$GENE_IDS" "$PROTEIN_FASTA" | grep -v "^--$" > "$SEQUENCES"

TOTAL_SEQS=$(grep -c "^>" "$SEQUENCES" || true)
echo "  ✓ Extracted $TOTAL_SEQS protein sequences"
echo "  ✓ Saved to: $SEQUENCES"
echo ""

# Check if counts match
if [ "$TOTAL_GENES" -ne "$TOTAL_SEQS" ]; then
    MISSING=$((TOTAL_GENES - TOTAL_SEQS))
    echo "⚠ WARNING: ID mismatch!"
    echo "  Gene IDs in PyDESeq2: $TOTAL_GENES"
    echo "  Sequences found:       $TOTAL_SEQS"
    echo "  Missing:               $MISSING"
    echo ""
    echo "This usually means:"
    echo "  - Gene IDs use different format (LOC123 vs LOC123.1)"
    echo "  - Some genes don't have protein sequences"
    echo ""
    echo "See BLAST_ALL_GENES_WORKFLOW.md 'TROUBLESHOOTING' section"
    echo ""
    
    # Find missing genes
    MISSING_IDS="$OUTPUT_DIR/missing_gene_ids.txt"
    grep "^>" "$SEQUENCES" | sed 's/^>//' | cut -d' ' -f1 | sort > "$OUTPUT_DIR/found_ids.txt"
    sort "$GENE_IDS" > "$OUTPUT_DIR/sorted_gene_ids.txt"
    comm -23 "$OUTPUT_DIR/sorted_gene_ids.txt" "$OUTPUT_DIR/found_ids.txt" > "$MISSING_IDS"
    
    echo "Missing gene IDs saved to: $MISSING_IDS"
    echo "First 10 missing genes:"
    head -10 "$MISSING_IDS"
    echo ""
else
    echo "✓ All genes found! Ready for BLAST"
    echo ""
fi

# Step 3: Create summary
SUMMARY="$OUTPUT_DIR/extraction_summary.txt"
cat > "$SUMMARY" << EOF
Sequence Extraction Summary
===========================

Date: $(date)

Input Files:
  PyDESeq2 results: $PYDESEQ_RESULTS
  Protein FASTA:    $PROTEIN_FASTA

Output Files:
  Gene IDs:         $GENE_IDS
  Sequences:        $SEQUENCES

Statistics:
  Total gene IDs:   $TOTAL_GENES
  Sequences found:  $TOTAL_SEQS
  Match rate:       $(awk "BEGIN {printf \"%.1f%%\", 100.0*$TOTAL_SEQS/$TOTAL_GENES}")

Next Steps:
==========

1. Run BLAST:
   sbatch scripts/run_cyp_blast.sbatch

2. Or manually:
   blastp \\
       -query $SEQUENCES \\
       -db /path/to/database \\
       -out $OUTPUT_DIR/blast_results.txt \\
       -outfmt "6 qseqid sseqid pident length evalue bitscore stitle" \\
       -evalue 1e-10 \\
       -num_threads 16

3. Combine BLAST results with PyDESeq2:
   python scripts/combine_blast_deseq.py

EOF

echo "=========================================="
echo "✓ Extraction complete!"
echo "=========================================="
echo ""
echo "Output files:"
echo "  Gene IDs:    $GENE_IDS"
echo "  Sequences:   $SEQUENCES"
echo "  Summary:     $SUMMARY"
echo ""
echo "Next step: Run BLAST"
echo "  sbatch scripts/run_cyp_blast.sbatch"
echo ""
