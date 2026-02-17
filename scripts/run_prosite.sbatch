#!/bin/bash
################################################################################
# Run PROSITE motif scanning for protein pattern identification
################################################################################
#
# Uses PROSITE patterns to identify functional motifs in protein sequences.
# PROSITE is complementary to HMMER (Pfam):
#   - PROSITE: Short, highly specific motifs/signatures (e.g., active sites)
#   - Pfam: Longer protein domains and families
#
# Use this AFTER:
#   1. Translating CDS to protein (run_translate_cds.sbatch)
#   2. Running HMMER (run_hmmer.sbatch) - for domain information
#   3. Optionally filtering to candidate genes
#
# How to run:
#   # Scan all proteins:
#   sbatch scripts/run_prosite.sbatch DC all_genes_protein.fasta
#
#   # Scan filtered proteins:
#   sbatch scripts/run_prosite.sbatch DC filtered_proteins.fasta
#
# Output:
#   - {filename}_prosite_motifs.txt (motif table)
#
################################################################################

#SBATCH --job-name=prosite_scan
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu

#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/06_analysis
#SBATCH -o prosite_%j.out
#SBATCH -e prosite_%j.err

# Strict error handling
set -euo pipefail

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq

# ============================================================================
# USER-EDITABLE VARIABLES
# ============================================================================
SPECIES="${1:-DC}"
QUERY_FILE="${2:-all_genes_protein.fasta}"

# Validate species
if [[ ! "$SPECIES" =~ ^(DC|DG|MF)$ ]]; then
    echo "ERROR: Invalid species '$SPECIES'. Must be DC, DG, or MF"
    echo "Usage: sbatch run_prosite.sbatch [SPECIES] [protein.fasta]"
    exit 1
fi

# Set paths
BASE_DIR=/projects/tholl_lab_1/daisy_analysis
INPUT_DIR="${BASE_DIR}/06_analysis/blast_input_${SPECIES}"
OUTPUT_DIR="${BASE_DIR}/06_analysis/prosite_${SPECIES}"
PROSITE_DIR="${BASE_DIR}/07_NRdatabase/prosite"

# Determine input file path
if [[ "$QUERY_FILE" == /* ]]; then
    # Absolute path provided
    QUERY_FASTA="$QUERY_FILE"
else
    # Relative path - look in blast_input directory
    QUERY_FASTA="${INPUT_DIR}/${QUERY_FILE}"
fi

# Validate input file exists
if [ ! -f "$QUERY_FASTA" ]; then
    echo "ERROR: Query file not found: $QUERY_FASTA"
    echo "Make sure to run: sbatch scripts/run_translate_cds.sbatch $SPECIES"
    exit 1
fi

# Validate PROSITE database exists
if [ ! -f "${PROSITE_DIR}/prosite.dat" ]; then
    echo "ERROR: PROSITE database not found: ${PROSITE_DIR}/prosite.dat"
    echo ""
    echo "Download PROSITE database:"
    echo "  sbatch scripts/download_protein_databases.sbatch"
    echo ""
    echo "Or manually:"
    echo "  mkdir -p ${PROSITE_DIR}"
    echo "  cd ${PROSITE_DIR}"
    echo "  wget ftp://ftp.expasy.org/databases/prosite/prosite.dat"
    echo "  wget ftp://ftp.expasy.org/databases/prosite/prosite.psa"
    exit 1
fi

# ============================================================================

# Create output directory
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

# Get base filename for output
BASENAME=$(basename "$QUERY_FASTA" .fasta)

# Print info
echo "=========================================="
echo "PROSITE Motif Scan Started"
echo "=========================================="
echo "Species:      $SPECIES"
echo "Query:        $QUERY_FASTA"
echo "PROSITE DB:   $PROSITE_DIR"
echo "Output dir:   $OUTPUT_DIR"
echo "CPUs:         ${SLURM_CPUS_PER_TASK:-1}"
echo "------------------------------------------"
date
echo "Hostname: $(hostname)"
echo "------------------------------------------"

# Check for ps_scan.pl (from EMBOSS or Perl Bio::Tools::Prosite)
if command -v ps_scan.pl &> /dev/null; then
    echo "Using ps_scan.pl: $(which ps_scan.pl)"
elif command -v prositecan &> /dev/null; then
    echo "Using prositecan: $(which prositecan)"
else
    echo "ERROR: No PROSITE scanning tool found"
    echo ""
    echo "Install with one of:"
    echo "  conda install -c bioconda emboss"
    echo "  conda install -c bioconda perl-bio-prosite"
    exit 1
fi

echo "=========================================="

# Run PROSITE scan
echo ""
echo "Running PROSITE scan..."
echo "This may take 30 min - 4 hours depending on input size"
echo ""

# Try ps_scan.pl first (most common)
if command -v ps_scan.pl &> /dev/null; then
    echo "Using ps_scan.pl for scanning..."
    
    # ps_scan.pl output format:
    # >sequence_id
    # PS00001 motif_name START-END score
    
    ps_scan.pl \
        -d "${PROSITE_DIR}" \
        -o gff \
        "$QUERY_FASTA" \
        > "${BASENAME}_prosite_motifs_raw.gff"
    
    # Convert to more readable table format
    echo "Converting to table format..."
    grep -v "^#" "${BASENAME}_prosite_motifs_raw.gff" | \
        awk -F'\t' 'BEGIN {
            print "sequence_id\tprosite_id\tmotif_name\tstart\tend\tscore\tstrand"
        } {
            # Extract PROSITE ID and name from column 9 (attributes)
            split($9, attr, ";")
            for (i in attr) {
                if (attr[i] ~ /Name=/) {
                    split(attr[i], name, "=")
                    prosite_name = name[2]
                }
                if (attr[i] ~ /ID=/) {
                    split(attr[i], id, "=")
                    prosite_id = id[2]
                }
            }
            print $1"\t"prosite_id"\t"prosite_name"\t"$4"\t"$5"\t"$6"\t"$7
        }' > "${BASENAME}_prosite_motifs.txt"
    
    # Clean up raw file
    rm -f "${BASENAME}_prosite_motifs_raw.gff"

# Alternative: try prositecan
elif command -v prositecan &> /dev/null; then
    echo "Using prositecan for scanning..."
    
    prositecan \
        -i "$QUERY_FASTA" \
        -d "${PROSITE_DIR}/prosite.dat" \
        -o "${BASENAME}_prosite_motifs.txt"

else
    echo "ERROR: No compatible PROSITE tool found"
    exit 1
fi

# Check if successful
if [ $? -ne 0 ]; then
    echo "ERROR: PROSITE scan failed"
    exit 1
fi

echo ""
echo "=========================================="
echo "PROSITE Scan Completed"
echo "=========================================="
echo "Output file:"
echo "  ${BASENAME}_prosite_motifs.txt"
echo ""

# Quick summary
if [ -f "${BASENAME}_prosite_motifs.txt" ]; then
    echo "Quick Summary:"
    echo "------------------------------------------"
    
    # Count proteins with motifs
    n_with_motifs=$(tail -n +2 "${BASENAME}_prosite_motifs.txt" | cut -f1 | sort -u | wc -l)
    echo "Proteins with motifs: $n_with_motifs"
    
    # Count total motif hits
    n_motifs=$(tail -n +2 "${BASENAME}_prosite_motifs.txt" | wc -l)
    echo "Total motif hits: $n_motifs"
    
    # Top 10 most common motifs
    echo ""
    echo "Top 10 most common motifs:"
    tail -n +2 "${BASENAME}_prosite_motifs.txt" | cut -f2,3 | sort | uniq -c | sort -rn | head -10
    echo "------------------------------------------"
fi

echo ""
date
echo "=========================================="
