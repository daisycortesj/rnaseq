#!/bin/bash
# =============================================================================
# SCRIPT:  run_interproscan.sbatch
#
# WHAT IT DOES:
#   Scans nutmeg proteins for domains, InterPro IDs, and GO terms.
#
# BEFORE YOU RUN THIS:
#   1. Download InterProScan once (skips if already there):
#        sbatch scripts/setup/download_interproscan.sbatch
#   2. Have TransDecoder proteins (not CD-HIT DNA):
#        sbatch scripts/03_assembly/run_transdecoder.sbatch MF
#
# HOW TO SUBMIT:
#   cd /projects/tholl_lab_1/daisy_analysis/05_rnaseq-code
#   sbatch scripts/07_domains/run_interproscan.sbatch
#
# Chain download → scan (scan waits until download finishes):
#   DL=$(sbatch --parsable scripts/setup/download_interproscan.sbatch)
#   sbatch --dependency=afterok:$DL scripts/07_domains/run_interproscan.sbatch
# =============================================================================

#SBATCH --job-name=InterProScan
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/06_analysis
#SBATCH --output=InterProScan_%j.out
#SBATCH --error=InterProScan_%j.err

set -euo pipefail


############################################################
## 1. Folders and files
############################################################

WORK_DIR="/projects/tholl_lab_1/daisy_analysis/06_analysis/interproscan_MF"
INSTALL_DIR="/projects/tholl_lab_1/daisy_analysis/04_reference/software"
VERSION="5.75-106.0"
INTERPRO_DIR="$INSTALL_DIR/interproscan-$VERSION"
INTERPRO_SH="$INTERPRO_DIR/interproscan.sh"

# Proteins from TransDecoder (NOT the CD-HIT .fasta DNA file)
INPUT_FASTA="/projects/tholl_lab_1/daisy_analysis/01_processed/00_6_cdhit/MF_trinity_cdhit95.fasta.transdecoder.pep"

OUTPUT_PREFIX="$WORK_DIR/MF_trinity_cdhit95_interproscan"
TEMP_DIR="$WORK_DIR/interproscan_temp_${SLURM_JOB_ID:-manual}"
CLEAN_FASTA="$WORK_DIR/MF_proteins_clean.faa"

mkdir -p "$WORK_DIR"
mkdir -p "$TEMP_DIR"


############################################################
## 2. Java 11 (InterProScan needs this version)
############################################################
# If this fails, on a login node run:  module spider Java
# Then change 11.0.20 to the Java 11 name your cluster shows.

module load Java/11.0.20 || module load Java/11 || true

echo "Java version:"
java -version


############################################################
## 3. Check InterProScan is already downloaded
############################################################

if [[ ! -f "$INTERPRO_SH" ]]; then
    echo "ERROR: InterProScan is not installed."
    echo "Download it once (the script skips if it is already there):"
    echo "  sbatch scripts/setup/download_interproscan.sbatch"
    exit 1
fi

echo "Using InterProScan: $INTERPRO_SH"
chmod +x "$INTERPRO_SH"


############################################################
## 4. Check the protein FASTA exists
############################################################

if [[ ! -f "$INPUT_FASTA" ]]; then
    echo "ERROR: Cannot find $INPUT_FASTA"
    echo "Run TransDecoder first:"
    echo "  sbatch scripts/03_assembly/run_transdecoder.sbatch MF"
    exit 1
fi


############################################################
## 5. Make a clean copy of the protein FASTA
############################################################
# TransDecoder files have two things InterProScan does not like:
#   - extra words after the protein ID on the >header line
#   - a * character (stop codon) at the end of the sequence
# We write a cleaned copy. Your original .pep file is not changed.

echo "Cleaning FASTA (strip * and extra header text)..."

awk '{
    if ($0 ~ /^>/) {
        print $1
    } else {
        gsub(/\*/, "")
        print
    }
}' "$INPUT_FASTA" > "$CLEAN_FASTA"

echo "Proteins to scan:"
grep -c "^>" "$CLEAN_FASTA"


############################################################
## 6. Remove old result files with the same names
############################################################

rm -f "${OUTPUT_PREFIX}.tsv"
rm -f "${OUTPUT_PREFIX}.gff3"


############################################################
## 7. Run InterProScan
############################################################
# -i  input proteins
# -b  output name (no extension)
# -f  make a TSV table and a GFF3 file
# -cpu  use the CPUs requested above
# -t p  input is protein (not DNA)
# -goterms / -iprlookup / -pa  add GO terms, InterPro IDs, pathways
# -dp  do not use UniProt precomputed hits (your proteins are new)
# -T  temp folder

cd "$WORK_DIR"

echo "Running InterProScan..."
echo "Started: $(date)"

"$INTERPRO_SH" \
    -i "$CLEAN_FASTA" \
    -b "$OUTPUT_PREFIX" \
    -f TSV,GFF3 \
    -cpu "${SLURM_CPUS_PER_TASK:-8}" \
    -t p \
    -goterms \
    -iprlookup \
    -pa \
    -dp \
    -T "$TEMP_DIR"

echo "Finished scan: $(date)"


############################################################
## 8. Make a simpler table (keep the useful columns only)
############################################################
# InterProScan TSV columns we keep:
#   $1 protein ID   $4 database   $5 signature   $6 description
#   $12 InterPro ID   $13 InterPro description   $14 GO terms

awk -F'\t' '
BEGIN {
    OFS = "\t"
    print "Protein", "Database", "Signature", "Description", \
          "InterPro_ID", "InterPro_description", "GO_terms"
}
{
    print $1, $4, $5, $6, $12, $13, $14
}
' "${OUTPUT_PREFIX}.tsv" \
> "$WORK_DIR/MF_interproscan_summary.tsv"


############################################################
## 9. Keep only oxidative-enzyme rows (CYP / P450 / oxidase)
############################################################

awk -F'\t' '
BEGIN { OFS = "\t" }
NR == 1 { print; next }
tolower($0) ~ /cytochrome|p450|cyp|monooxygenase|oxidase|oxygenase|hydroxylase|demethylase|epoxidase|dioxygenase|peroxidase|oxidoreductase|flavin|fadh|fmn|ferredoxin|reductase/ {
    print
}
' "$WORK_DIR/MF_interproscan_summary.tsv" \
> "$WORK_DIR/MF_interproscan_oxidative.tsv"


############################################################
## 10. Clean up and print what was created
############################################################

rm -rf "$TEMP_DIR"

echo "DONE"
echo "Results:"
echo "${OUTPUT_PREFIX}.tsv"
echo "${OUTPUT_PREFIX}.gff3"
echo "$WORK_DIR/MF_interproscan_summary.tsv"
echo "$WORK_DIR/MF_interproscan_oxidative.tsv"
