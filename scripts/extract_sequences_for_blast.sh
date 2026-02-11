#!/bin/bash
#
# =============================================================================
#  EXTRACT PROTEIN SEQUENCES FOR BLAST 
# =============================================================================
#
#  WHAT THIS SCRIPT DOES:
#  ----------------------------------------
#  1. You already ran PyDESeq2 and have a spreadsheet of genes (a .tsv file).
#  2. You also have a big file with ALL protein sequences for your species.
#  3. This script takes the LIST of genes from the spreadsheet and pulls out
#     ONLY those proteins from the big file. So you get one smaller file with
#     just the proteins you care about — ready to run BLAST on.
#
#  WHAT YOU NEED BEFORE RUNNING:
#  ----------------------------
#  - PyDESeq2 Step 1 must be done for your species (you have a ..._UNFILTERED.tsv)
#  - The protein file for that species must be in 04_reference/ (e.g. ..._protein.faa)
#
#  HOW TO RUN (pick ONE species code):
#  -----------------------------------
#    bash scripts/extract_sequences_for_blast.sh DC    (carrot - D. carota)
#    bash scripts/extract_sequences_for_blast.sh DG    (Daucus glaber)
#    bash scripts/extract_sequences_for_blast.sh MF    (nutmeg - M. fragrans)
#
#  You type DC, DG, or MF so the script knows which spreadsheet and which
#  protein file to use. Each species has its own folder and output file.
#
# =============================================================================

# If any command fails, stop the script (safer for beginners)
set -e

# Where is your project? (Default: current folder. You can set BASE_DIR if needed.)
BASE_DIR="${BASE_DIR:-.}"

# -----------------------------------------------------------------------------
#  PART 1: Which species are you working with?
# -----------------------------------------------------------------------------
# The first thing you type when you run the script (e.g. DC) is stored here.
# We use it to choose the right spreadsheet and the right protein file.

SPECIES="${1:-}"

if [ -n "${SPECIES}" ]; then
    # Make it uppercase so "dc" and "DC" both work
    SPECIES=$(echo "${SPECIES}" | tr '[:lower:]' '[:upper:]')

    # Set file paths based on species
    case "${SPECIES}" in
        DC)
            # Carrot: paths to carrot PyDESeq2 results and carrot protein file
            PYDESEQ_RESULTS="${PYDESEQ_RESULTS:-${BASE_DIR}/06_analysis/pydeseq2_DC_step1_unfiltered/pydeseq2_results_UNFILTERED.tsv}"
            PROTEIN_FASTA="${PROTEIN_FASTA:-${BASE_DIR}/04_reference/GCF_001625215.2_DH1_v3.0_protein.faa}"
            OUTPUT_DIR="${OUTPUT_DIR:-${BASE_DIR}/06_analysis/blast_input_DC}"
            SPECIES_NAME="D. carota ssp. maximus"
            ;;
        DG)
            # Daucus glaber: its own PyDESeq2 folder and same carrot protein ref if applicable
            PYDESEQ_RESULTS="${PYDESEQ_RESULTS:-${BASE_DIR}/06_analysis/pydeseq2_DG_step1_unfiltered/pydeseq2_results_UNFILTERED.tsv}"
            PROTEIN_FASTA="${PROTEIN_FASTA:-${BASE_DIR}/04_reference/GCF_001625215.2_DH1_v3.0_protein.faa}"
            OUTPUT_DIR="${OUTPUT_DIR:-${BASE_DIR}/06_analysis/blast_input_DG}"
            SPECIES_NAME="Daucus glaber"
            ;;
        MF)
            # Nutmeg: its own PyDESeq2 folder and nutmeg protein file
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
    # No species given: use default paths (one generic folder)
    PYDESEQ_RESULTS="${PYDESEQ_RESULTS:-${BASE_DIR}/06_analysis/pydeseq2_step1_unfiltered/pydeseq2_results_UNFILTERED.tsv}"
    PROTEIN_FASTA="${PROTEIN_FASTA:-${BASE_DIR}/04_reference/GCF_001625215.2_DH1_v3.0_protein.faa}"
    OUTPUT_DIR="${OUTPUT_DIR:-${BASE_DIR}/06_analysis/blast_input}"
    SPECIES_NAME="(default)"
fi

# -----------------------------------------------------------------------------
#  SHOW WHAT WE'RE DOING
# -----------------------------------------------------------------------------
echo ""
echo "=========================================="
echo "  Extract Sequences for BLAST"
echo "=========================================="
echo ""
echo "  Species: ${SPECIES_NAME}"
echo ""

# -----------------------------------------------------------------------------
#  PART 2: Check that your input files exist
# -----------------------------------------------------------------------------
# If the spreadsheet or protein file is missing, we stop and tell you instead
# of failing later with a confusing error.

if [ ! -f "$PYDESEQ_RESULTS" ]; then
    echo "  ERROR: PyDESeq2 results file not found."
    echo "  Expected: $PYDESEQ_RESULTS"
    echo ""
    echo "  Run PyDESeq2 Step 1 first for this species!"
    exit 1
fi

if [ ! -f "$PROTEIN_FASTA" ]; then
    echo "  ERROR: Protein file not found."
    echo "  Expected: $PROTEIN_FASTA"
    echo ""
    echo "  You need to download the protein FASTA for this species from NCBI"
    echo "  and put it in 04_reference/. See the guide for the download link."
    exit 1
fi

# Create the folder where we'll save results (no error if it already exists)
mkdir -p "$OUTPUT_DIR"

echo "  Input files we're using:"
echo "    Spreadsheet (gene list): $PYDESEQ_RESULTS"
echo "    Protein file:            $PROTEIN_FASTA"
echo "    Output folder:           $OUTPUT_DIR"
echo ""

# -----------------------------------------------------------------------------
#  PART 3: Step 1 — Get the list of gene IDs from the spreadsheet
# -----------------------------------------------------------------------------
# The spreadsheet has many columns. We only need the first column (gene ID)
# and we skip rows where the stats are "NA" (no number). We save that list
# to a simple text file: one gene ID per line.

echo "  Step 1: Reading the spreadsheet and making a list of gene IDs..."
echo "          (We skip the header row and rows where padj is NA.)"
echo ""

GENE_IDS="$OUTPUT_DIR/all_gene_ids.txt"

# awk: read the TSV; skip line 1 (header); keep rows where column 6 is not "NA"; print column 1
awk -F'\t' 'NR>1 && $6!="NA" {print $1}' "$PYDESEQ_RESULTS" > "$GENE_IDS"

TOTAL_GENES=$(wc -l < "$GENE_IDS")
echo "    Done. We found $TOTAL_GENES gene IDs."
echo "    Saved list to: $GENE_IDS"
echo ""

# -----------------------------------------------------------------------------
#  PART 4: Step 2 — Pull out only those proteins from the big protein file
# -----------------------------------------------------------------------------
# We take the list of gene IDs and search the big protein file. For every ID
# that matches a line starting with ">", we keep that line and the next line
# (the sequence). That gives us a small FASTA file with only our genes.

echo "  Step 2: Pulling out protein sequences for those genes..."
echo "          (Searching the big protein file for each gene ID.)"
echo ""

if [ -n "${SPECIES}" ]; then
    SEQUENCES="$OUTPUT_DIR/${SPECIES}_gene.faa"
else
    SEQUENCES="$OUTPUT_DIR/all_genes.faa"
fi

# grep -f: use each line of the gene list as a search string
# grep -A1: when we find a match, also print the next line (the sequence)
grep -A1 -Ff "$GENE_IDS" "$PROTEIN_FASTA" | grep -v "^--$" > "$SEQUENCES"

TOTAL_SEQS=$(grep -c "^>" "$SEQUENCES" || true)
echo "    Done. We found $TOTAL_SEQS protein sequences."
echo "    Saved to: $SEQUENCES"
echo ""

# -----------------------------------------------------------------------------
#  PART 5: Check if the numbers make sense
# -----------------------------------------------------------------------------
# Ideally the number of gene IDs and the number of sequences should match.
# If not, we write which genes were missing and show a few.

if [ "$TOTAL_GENES" -ne "$TOTAL_SEQS" ]; then
    MISSING=$((TOTAL_GENES - TOTAL_SEQS))
    echo "  WARNING: The numbers don't match."
    echo "    Gene IDs in spreadsheet: $TOTAL_GENES"
    echo "    Sequences we found:      $TOTAL_SEQS"
    echo "    Missing:                 $MISSING"
    echo ""
    echo "  This can happen if gene IDs in the spreadsheet don't exactly match"
    echo "  the headers in the protein file (e.g. different format)."
    echo ""

    # Save the list of missing gene IDs to a file so you can look at them
    MISSING_IDS="$OUTPUT_DIR/missing_gene_ids.txt"
    grep "^>" "$SEQUENCES" | sed 's/^>//' | cut -d' ' -f1 | sort > "$OUTPUT_DIR/found_ids.txt"
    sort "$GENE_IDS" > "$OUTPUT_DIR/sorted_gene_ids.txt"
    comm -23 "$OUTPUT_DIR/sorted_gene_ids.txt" "$OUTPUT_DIR/found_ids.txt" > "$MISSING_IDS"

    echo "  Missing gene IDs were saved to: $MISSING_IDS"
    echo "  First 10 missing:"
    head -10 "$MISSING_IDS"
    echo ""
else
    echo "  All genes had a matching sequence. You're ready for BLAST!"
    echo ""
fi

# -----------------------------------------------------------------------------
#  PART 6: Write a short summary file
# -----------------------------------------------------------------------------
# So you have a record of what was run and where the files are.

SUMMARY="$OUTPUT_DIR/extraction_summary.txt"
cat > "$SUMMARY" << EOF
Sequence Extraction Summary
===========================

Date: $(date)

What we used:
  - Spreadsheet (PyDESeq2): $PYDESEQ_RESULTS
  - Protein file:          $PROTEIN_FASTA

What we created:
  - Gene ID list:    $GENE_IDS
  - Protein FASTA:   $SEQUENCES

Counts:
  - Gene IDs:        $TOTAL_GENES
  - Sequences found: $TOTAL_SEQS

Next steps:
  1. Run BLAST using the file: $SEQUENCES
  2. Then combine BLAST results with PyDESeq2: python scripts/combine_blast_deseq.py
EOF

# -----------------------------------------------------------------------------
#  DONE — Tell the user clearly
# -----------------------------------------------------------------------------
echo "=========================================="
echo "  Extraction complete!"
echo "=========================================="
echo ""
echo "  Your output files:"
echo "    Gene list:   $GENE_IDS"
echo "    Proteins:    $SEQUENCES   <-- use this file for BLAST"
echo "    Summary:     $SUMMARY"
echo ""
echo "  Next step: Run BLAST (e.g. blastp) with $SEQUENCES as the query."
echo ""
