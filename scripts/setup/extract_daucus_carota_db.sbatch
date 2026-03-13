#!/bin/bash
################################################################################
# Extract Daucus carota Proteins from NR and Build Custom BLAST Database
################################################################################
#
# Extracts all Daucus carota (taxid 79200) protein sequences from your
# existing nr database and builds a small, fast custom BLAST database.
#
# Prerequisites:
#   - nr database must already be downloaded (run download_blastdb.sbatch first)
#   - taxdb must be present in the BLASTDB directory
#
# The resulting database is placed in the same BLASTDB directory so all
# existing BLAST scripts can use it by passing "daucus_carota" as the database.
#
# Usage:
#   sbatch scripts/extract_daucus_carota_db.sbatch
#
# Runtime: ~5-15 minutes
#
# After this completes, run BLAST against it:
#   sbatch scripts/blastp_discoveryfilter.sbatch DC DCdb
#   sbatch scripts/blastx_discoveryfilter.sbatch DC DCdb
#   sbatch scripts/blastp_strictfilter.sbatch DC DCdb
#   sbatch scripts/blastx_strictfilter.sbatch DC DCdb
#
################################################################################

# ══════════════════════════════════════════════════════════════════════════════
# SLURM SETTINGS (lightweight - small database extraction)
# ══════════════════════════════════════════════════════════════════════════════

#SBATCH --job-name=extract_daucus_db
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daisycortesj@vt.edu

#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/07_NRdatabase
#SBATCH -o extract_daucus_db_%j.out
#SBATCH -e extract_daucus_db_%j.err

set -euo pipefail

echo "=========================================="
echo "Daucus carota Database Extraction Started"
echo "=========================================="
echo "Date: $(date)"
echo "Hostname: $(hostname)"
echo ""

# Initialize conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq

# Database locations
BLASTDB_DIR=/projects/tholl_lab_1/daisy_analysis/07_NRdatabase/blastdb
export BLASTDB="$BLASTDB_DIR"

TAXID=79200
DB_NAME=DCdb
FASTA_FILE="${BLASTDB_DIR}/${DB_NAME}.fasta"

cd "$BLASTDB_DIR"

echo "BLASTDB directory: $BLASTDB_DIR"
echo "Daucus carota NCBI taxid: $TAXID"
echo "Output database name: $DB_NAME"
echo ""

# ══════════════════════════════════════════════════════════════════════════════
# VERIFY PREREQUISITES
# ══════════════════════════════════════════════════════════════════════════════

echo "=========================================="
echo "Checking Prerequisites"
echo "=========================================="

echo "Tools:"
echo "  blastdbcmd: $(which blastdbcmd)"
echo "  makeblastdb: $(which makeblastdb)"
echo ""

if ! blastdbcmd -db nr -info &>/dev/null; then
    echo "ERROR: nr database not found or not ready in $BLASTDB_DIR"
    echo "Run download_blastdb.sbatch first to download nr."
    exit 1
fi
echo "  ✓ nr database found"

if [ ! -f "${BLASTDB_DIR}/taxdb.btd" ] || [ ! -f "${BLASTDB_DIR}/taxdb.bti" ]; then
    echo "ERROR: taxdb not found in $BLASTDB_DIR"
    echo "Run download_blastdb.sbatch first (it downloads taxdb)."
    exit 1
fi
echo "  ✓ taxdb found"
echo ""

# ══════════════════════════════════════════════════════════════════════════════
# CHECK IF ALREADY BUILT
# ══════════════════════════════════════════════════════════════════════════════

if blastdbcmd -db "$DB_NAME" -info &>/dev/null; then
    echo "=========================================="
    echo "Database already exists!"
    echo "=========================================="
    blastdbcmd -db "$DB_NAME" -info
    echo ""
    echo "To rebuild, delete the existing files first:"
    echo "  rm -f ${BLASTDB_DIR}/${DB_NAME}.p*"
    echo "Then resubmit this script."
    echo ""
    echo "=========================================="
    echo "Done (nothing to do)"
    echo "Date: $(date)"
    echo "=========================================="
    exit 0
fi

# ══════════════════════════════════════════════════════════════════════════════
# STEP 1: EXTRACT DAUCUS CAROTA SEQUENCES FROM NR
# ══════════════════════════════════════════════════════════════════════════════

echo "=========================================="
echo "Step 1: Extracting Daucus carota from nr"
echo "=========================================="
echo "Start time: $(date)"
echo ""

echo "Running: blastdbcmd -db nr -taxids $TAXID -out $FASTA_FILE"
blastdbcmd -db nr -taxids "$TAXID" -out "$FASTA_FILE"

if [ ! -f "$FASTA_FILE" ] || [ ! -s "$FASTA_FILE" ]; then
    echo "ERROR: Extraction produced no sequences"
    echo "This could mean:"
    echo "  - taxdb is missing or outdated"
    echo "  - nr doesn't contain Daucus carota sequences (unlikely)"
    exit 1
fi

SEQ_COUNT=$(grep -c '^>' "$FASTA_FILE")
FILE_SIZE=$(du -h "$FASTA_FILE" | cut -f1)
echo ""
echo "✓ Extracted $SEQ_COUNT Daucus carota protein sequences ($FILE_SIZE)"
echo "End time: $(date)"
echo ""

# ══════════════════════════════════════════════════════════════════════════════
# STEP 2: BUILD BLAST DATABASE
# ══════════════════════════════════════════════════════════════════════════════

echo "=========================================="
echo "Step 2: Building BLAST database"
echo "=========================================="
echo "Start time: $(date)"
echo ""

makeblastdb \
    -in "$FASTA_FILE" \
    -dbtype prot \
    -out "${BLASTDB_DIR}/${DB_NAME}" \
    -title "Daucus carota proteins from nr (taxid ${TAXID})" \
    -parse_seqids \
    -taxid "$TAXID"

echo ""
echo "✓ Database built"
echo "End time: $(date)"
echo ""

# ══════════════════════════════════════════════════════════════════════════════
# STEP 3: VERIFY
# ══════════════════════════════════════════════════════════════════════════════

echo "=========================================="
echo "Step 3: Verification"
echo "=========================================="

if blastdbcmd -db "$DB_NAME" -info &>/dev/null; then
    echo "✓ Database is valid and ready to use"
    echo ""
    blastdbcmd -db "$DB_NAME" -info
else
    echo "ERROR: Database verification failed"
    exit 1
fi

echo ""

# Clean up the intermediate FASTA (the .p* database files are what BLAST uses)
rm -f "$FASTA_FILE"
echo "Cleaned up intermediate FASTA file"
echo ""

# ══════════════════════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════════════════════

echo "=========================================="
echo "SUMMARY"
echo "=========================================="
echo ""
echo "Database:   $DB_NAME"
echo "Location:   $BLASTDB_DIR"
echo "Sequences:  $SEQ_COUNT Daucus carota proteins"
echo "Taxid:      $TAXID"
echo ""
echo "Database files:"
ls -lh "${BLASTDB_DIR}/${DB_NAME}".p* 2>/dev/null
echo ""
echo "Disk usage:"
du -sh "${BLASTDB_DIR}/${DB_NAME}".p* 2>/dev/null | tail -1
echo ""

echo "=========================================="
echo "Next Steps"
echo "=========================================="
echo ""
echo "Run BLAST against the Daucus carota database:"
echo ""
echo "  # BLASTp (protein vs protein) - fast, use with translated proteins"
echo "  sbatch scripts/blastp_discoveryfilter.sbatch DC DCdb"
echo "  sbatch scripts/blastp_strictfilter.sbatch DC DCdb"
echo ""
echo "  # BLASTX (nucleotide vs protein) - use with CDS nucleotide sequences"
echo "  sbatch scripts/blastx_discoveryfilter.sbatch DC DCdb"
echo "  sbatch scripts/blastx_strictfilter.sbatch DC DCdb"
echo ""
echo "  # TIP: DCdb is tiny, so you can override SLURM resources:"
echo "  sbatch --mem=8G --cpus-per-task=4 --time=1:00:00 scripts/blastp_discoveryfilter.sbatch DC DCdb"
echo ""

echo "=========================================="
echo "Daucus carota Database Extraction Complete"
echo "Date: $(date)"
echo "=========================================="
