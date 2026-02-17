#!/bin/bash
################################################################################
# Complete BLAST Input Preparation Workflow (SLURM)
################################################################################
#
# This script prepares the input for BLAST analysis by:
# What this does:
#   1. Parse GTF → extract gene annotations WITH transcript IDs (XM numbers)
#   2. Inner join → filter for only your genes from all_gene_ids.txt
#   3. Extract CDS → create FASTA with enhanced headers for BLAST
# 
# Uses:
#   - parse_refseq_gtf.py - extracts gene annotations with transcript IDs (XM numbers)
#   - join_gtf_with_gene_ids.py - filters for only your genes from all_gene_ids.txt
#   - extract_cds_from_gtf.py - extracts CDS sequences from the GTF file and creates a FASTA file with enhanced headers for BLAST
#   
#
# How to run:
#   sbatch scripts/blast_input.sbatch DC
#   sbatch scripts/blast_input.sbatch DG
#   sbatch scripts/blast_input.sbatch MF
#
# What you get:
#   - filtered_genes_complete.tsv (annotations with XM numbers)
#   - all_genes_cds.fasta (sequences with enhanced headers)
#
################################################################################

#SBATCH --job-name=extract_cds           # Job name
#SBATCH --account=tholl_lab_1            # Your lab account
#SBATCH --partition=normal_q             # Which queue to use
#SBATCH --cpus-per-task=4                # Number of CPUs
#SBATCH --mem=16G                        # Memory (genome loading)
#SBATCH --time=1:00:00                   # Max time (1 hour)
#SBATCH --mail-type=END,FAIL             # Email when done or failed
#SBATCH --mail-user=daisycortesj@vt.edu  # Your email

# Where to save output and error messages
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis
#SBATCH -o 06_analysis/extract_cds_%j.out
#SBATCH -e 06_analysis/extract_cds_%j.err


# ══════════════════════════════════════════════════════════════════════════════
# PATHS
# ══════════════════════════════════════════════════════════════════════════════

BASE_DIR="/projects/tholl_lab_1/daisy_analysis"
CODE_DIR="${BASE_DIR}/05_rnaseq-code"


# ══════════════════════════════════════════════════════════════════════════════
# GET SPECIES FROM COMMAND LINE
# ══════════════════════════════════════════════════════════════════════════════

SPECIES=$1

if [ -z "${SPECIES}" ]; then
    echo "ERROR: You need to tell me which species!"
    echo ""
    echo "Run like this:"
    echo "  sbatch scripts/extract_cds_for_blast.sbatch DC"
    echo "  sbatch scripts/extract_cds_for_blast.sbatch DG"
    echo "  sbatch scripts/extract_cds_for_blast.sbatch MF"
    exit 1
fi

# Convert to uppercase
SPECIES=$(echo "${SPECIES}" | tr '[:lower:]' '[:upper:]')

# Set the species name
if [ "${SPECIES}" = "DC" ]; then
    SPECIES_NAME="D. carota ssp. maximus"
elif [ "${SPECIES}" = "DG" ]; then
    SPECIES_NAME="Daucus glaber"
elif [ "${SPECIES}" = "MF" ]; then
    SPECIES_NAME="Myristica fragrans"
else
    echo "ERROR: Unknown species: ${SPECIES}"
    exit 1
fi


# ══════════════════════════════════════════════════════════════════════════════
# SHOW WHAT WE'RE DOING
# ══════════════════════════════════════════════════════════════════════════════

echo "============================================================"
echo "Extract CDS for BLAST"
echo "============================================================"
echo ""
echo "Species: ${SPECIES_NAME} (${SPECIES})"
echo "Start time: $(date)"
echo ""


# ══════════════════════════════════════════════════════════════════════════════
# ACTIVATE CONDA ENVIRONMENT
# ══════════════════════════════════════════════════════════════════════════════

source ~/.bashrc
conda activate rnaseq

if [ -z "${CONDA_DEFAULT_ENV}" ]; then
    echo "ERROR: Conda environment didn't activate"
    exit 1
fi

echo "Using conda environment: ${CONDA_DEFAULT_ENV}"
echo ""


# ══════════════════════════════════════════════════════════════════════════════
# RUN COMPLETE BLAST INPUT PREPARATION WORKFLOW
# ══════════════════════════════════════════════════════════════════════════════

echo "============================================================"
echo "STEP 1/3: Parse GTF with Auto Transcript Lookup"
echo "============================================================"
echo ""
echo "Extracting gene annotations WITH transcript IDs (XM numbers)..."
echo ""

cd "${BASE_DIR}"

# Determine input/output paths based on species
if [ "${SPECIES}" = "DC" ]; then
    GTF="04_reference/dc_genomic.gtf"
    GENOME="04_reference/GCF_001625215.2_DH1_v3.0_genomic.fna"
    GENE_IDS="06_analysis/pydeseq2_DC_step1_unfiltered/all_gene_ids.txt"
    OUTPUT_DIR="06_analysis/blast_input_DC"
elif [ "${SPECIES}" = "DG" ]; then
    GTF="04_reference/dc_genomic.gtf"
    GENOME="04_reference/GCF_001625215.2_DH1_v3.0_genomic.fna"
    GENE_IDS="06_analysis/pydeseq2_DG_step1_unfiltered/all_gene_ids.txt"
    OUTPUT_DIR="06_analysis/blast_input_DG"
elif [ "${SPECIES}" = "MF" ]; then
    GTF="04_reference/mf_genomic.gtf"
    GENOME="04_reference/MYU_GWHGDHP00000000.1_genomic.fna"
    GENE_IDS="06_analysis/pydeseq2_MF_step1_unfiltered/all_gene_ids.txt"
    OUTPUT_DIR="06_analysis/blast_input_MF"
fi

GENES_PARSED="${OUTPUT_DIR}/genes_with_XM.tsv"
FILTERED_GENES="${OUTPUT_DIR}/filtered_genes_complete.tsv"
FILTERED_IDS="${OUTPUT_DIR}/filtered_gene_ids.txt"
FASTA_OUTPUT="${OUTPUT_DIR}/all_genes_cds.fasta"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Step 1: Parse GTF with smart parser
python3 "${CODE_DIR}/scripts/parse_refseq_gtf.py" \
    --gtf "${GTF}" \
    --out "${GENES_PARSED}"

if [ $? -ne 0 ]; then
    echo "ERROR: GTF parsing failed"
    exit 1
fi

echo ""
echo "============================================================"
echo "STEP 2/3: Filter Genes (Inner Join)"
echo "============================================================"
echo ""
echo "Keeping only genes in your analysis list..."
echo ""

python3 "${CODE_DIR}/scripts/join_gtf_with_gene_ids.py" \
    --parsed "${GENES_PARSED}" \
    --gene-ids "${GENE_IDS}" \
    --out "${FILTERED_GENES}"

if [ $? -ne 0 ]; then
    echo "ERROR: Inner join failed"
    exit 1
fi

# Extract gene IDs from filtered TSV for sequence extraction
echo ""
echo "Extracting gene IDs from filtered results..."
tail -n +2 "${FILTERED_GENES}" | cut -f2 > "${FILTERED_IDS}"
echo "  Created: ${FILTERED_IDS}"

echo ""
echo "============================================================"
echo "STEP 3/3: Extract CDS Sequences"
echo "============================================================"
echo ""
echo "Extracting sequences with enhanced FASTA headers..."
echo ""

python3 "${CODE_DIR}/scripts/extract_cds_from_gtf.py" \
    "${GTF}" \
    "${GENOME}" \
    "${FILTERED_IDS}" \
    "${FASTA_OUTPUT}"

if [ $? -ne 0 ]; then
    echo "ERROR: CDS extraction failed"
    exit 1
fi


# ══════════════════════════════════════════════════════════════════════════════
# WORKFLOW COMPLETE!
# ══════════════════════════════════════════════════════════════════════════════

echo ""
echo "============================================================"
echo "BLAST INPUT PREPARATION COMPLETE!"
echo "============================================================"
echo ""
echo "Species: ${SPECIES_NAME} (${SPECIES})"
echo ""
echo "Output files created:"
echo "  1. ${GENES_PARSED}"
echo "     → All genes with transcript IDs (XM numbers)"
echo ""
echo "  2. ${FILTERED_GENES}"
echo "     → Only your genes with complete annotations"
echo ""
echo "  3. ${FASTA_OUTPUT}"
echo "     → CDS sequences with enhanced headers"
echo "     Format: >gene_id|transcript_id description"
echo ""
echo "------------------------------------------------------------"
echo "NEXT STEPS"
echo "------------------------------------------------------------"
echo ""
echo "1. Run BLAST (command-line):"
echo "   blastx -query ${FASTA_OUTPUT} \\"
echo "          -db nr \\"
echo "          -out 06_analysis/blast_results_${SPECIES}.txt \\"
echo "          -outfmt 6 \\"
echo "          -evalue 1e-5 \\"
echo "          -num_threads 16"
echo ""
echo "2. Or download FASTA to use Geneious:"
echo "   scp daisycortesj@tinkercliffs.arc.vt.edu:${BASE_DIR}/${FASTA_OUTPUT} ~/Desktop/"
echo ""
echo "3. Combine BLAST results with annotations:"
echo "   python scripts/combine_blast_deseq.py ..."
echo ""
echo "End time: $(date)"
echo "============================================================"

exit 0
