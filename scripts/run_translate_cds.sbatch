#!/bin/bash
################################################################################
# Translate CDS sequences to protein sequences
################################################################################
#
# Converts all_genes_cds.fasta to all_genes_protein.fasta for each species.
# Uses Biopython to translate nucleotide CDS to amino acid sequences.
#
# How to run:
#   sbatch scripts/run_translate_cds.sbatch DC
#   sbatch scripts/run_translate_cds.sbatch DG
#   sbatch scripts/run_translate_cds.sbatch MF
#   sbatch scripts/run_translate_cds.sbatch all  # Translate all species
#
################################################################################

#SBATCH --job-name=translate_cds
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu

#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/06_analysis
#SBATCH -o translate_cds_%j.out
#SBATCH -e translate_cds_%j.err

# Strict error handling
set -euo pipefail


# install biopython:
pip install biopython



# Get species argument (default to "all")
SPECIES="${1:-all}"

echo "=========================================="
echo "CDS → Protein Translation"
echo "=========================================="
echo "Date: $(date)"
echo "Hostname: $(hostname)"
echo "Species: $SPECIES"
echo ""

# Check if biopython is installed
echo "Checking Biopython installation..."
python -c "from Bio import SeqIO; print('✓ Biopython is installed')" || {
    echo "ERROR: Biopython not found in rnaseq environment"
    echo "Install with: conda install -c conda-forge biopython"
    exit 1
}
echo ""

# Set paths
BASE_DIR=/projects/tholl_lab_1/daisy_analysis
CODE_DIR=${BASE_DIR}/05_rnaseq-code

# Function to translate one species
translate_species() {
    local sp=$1
    echo "=========================================="
    echo "Translating $sp"
    echo "=========================================="
    python ${CODE_DIR}/scripts/translate_cds_to_protein.py "$sp"
    echo ""
}

# Translate based on argument
if [ "$SPECIES" == "all" ]; then
    echo "Translating all species (DC, DG, MF)..."
    echo ""
    translate_species DC
    translate_species DG
    translate_species MF
else
    # Validate species
    if [[ ! "$SPECIES" =~ ^(DC|DG|MF)$ ]]; then
        echo "ERROR: Invalid species '$SPECIES'. Must be DC, DG, MF, or all"
        echo "Usage: sbatch run_translate_cds.sbatch [SPECIES]"
        echo "  SPECIES: DC, DG, MF, or all"
        exit 1
    fi
    translate_species "$SPECIES"
fi

echo "=========================================="
echo "Translation Complete"
echo "Date: $(date)"
echo "=========================================="
echo ""
echo "Next steps:"
echo "  1. Run BLASTp (protein vs protein) for faster searches:"
echo "     sbatch scripts/blastp_strictfilter.sbatch DC nr"
echo "     sbatch scripts/blastp_strictfilter.sbatch DG swissprot"
echo "     sbatch scripts/blastp_strictfilter.sbatch MF refseq_protein"
echo ""
