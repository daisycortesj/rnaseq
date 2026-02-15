#!/bin/bash
################################################################################
# Extract CDS sequences for BLAST
################################################################################
#
# What this does:
#   Extracts CDS (coding sequence) from GTF + genome for genes in PyDESeq2 results
#   Creates nucleotide FASTA for BLAST (blastx vs NR)
#
# How to run:
#   sbatch scripts/extract_cds_for_blast.sbatch DC    (for D. carota)
#   sbatch scripts/extract_cds_for_blast.sbatch DG    (for D. glaber)
#   sbatch scripts/extract_cds_for_blast.sbatch MF    (for M. fragrans)
#
# What you get:
#   06_analysis/blast_input_DC/all_genes_cds.fasta
#   (Ready for BLAST or import into Geneious)
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
# RUN CDS EXTRACTION
# ══════════════════════════════════════════════════════════════════════════════

echo "============================================================"
echo "EXTRACTING CDS SEQUENCES"
echo "============================================================"
echo ""
echo "This will:"
echo "  1. Read gene IDs from PyDESeq2 results"
echo "  2. Find CDS coordinates in GTF"
echo "  3. Extract sequences from genome"
echo "  4. Write FASTA file for BLAST"
echo ""

# Run the Python script with species code
cd "${BASE_DIR}"

python3 "${CODE_DIR}/scripts/extract_cds_from_gtf.py" "${SPECIES}"

# Check if it worked
if [ $? -ne 0 ]; then
    echo ""
    echo "ERROR: CDS extraction failed"
    exit 1
fi


# ══════════════════════════════════════════════════════════════════════════════
# DONE!
# ══════════════════════════════════════════════════════════════════════════════

echo ""
echo "============================================================"
echo "CDS EXTRACTION COMPLETE!"
echo "============================================================"
echo ""
echo "Species: ${SPECIES_NAME}"
echo ""
echo "Output: 06_analysis/blast_input_${SPECIES}/all_genes_cds.fasta"
echo ""
echo "------------------------------------------------------------"
echo "NEXT STEPS"
echo "------------------------------------------------------------"
echo ""
echo "1. Download FASTA to your computer:"
echo "   scp daisycortesj@tinkercliffs.arc.vt.edu:${BASE_DIR}/06_analysis/blast_input_${SPECIES}/all_genes_cds.fasta ~/Desktop/"
echo ""
echo "2. Import into Geneious and run BLAST"
echo "   See: GENEOUS_BEGINNER_GUIDE.md"
echo ""
echo "3. Or run command-line BLAST:"
echo "   blastx -query 06_analysis/blast_input_${SPECIES}/all_genes_cds.fasta -db nr ..."
echo ""
echo "End time: $(date)"
echo "============================================================"

exit 0
