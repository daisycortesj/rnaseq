#!/bin/bash
################################################################################
# Download and Prepare Protein Domain/Motif Databases
################################################################################
#
# This script downloads and prepares protein domain and motif databases:
#   - Pfam-A (for HMMER, ~600 MB compressed, ~2.3 GB uncompressed)
#   - PROSITE (for motif scanning, ~30 MB)
#
# Features:
#   - Automatically extracts and indexes databases
#   - Skips already-prepared databases
#   - Verifies databases after processing
#   - Uses the same directory as BLAST databases for consistency
#
# Usage:
#   sbatch scripts/download_protein_databases.sbatch
#
# Runtime: 10-60 minutes depending on download speed
#
################################################################################

# ══════════════════════════════════════════════════════════════════════════════
# SLURM SETTINGS 
# ══════════════════════════════════════════════════════════════════════════════

#SBATCH --job-name=protein_db_download        # Job name
#SBATCH --account=tholl_lab_1                 # Your lab account
#SBATCH --partition=normal_q                  # Which queue to use
#SBATCH --cpus-per-task=4                     # Number of CPUs
#SBATCH --mem=16G                             # Memory
#SBATCH --time=4:00:00                        # Max time (4 hours)
#SBATCH --mail-type=ALL                       # Email when done
#SBATCH --mail-user=daisycortesj@vt.edu       # Your email

# Where to save output and error messages
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/07_NRdatabase
#SBATCH -o protein_db_download_%j.out
#SBATCH -e protein_db_download_%j.err

# Exit on error, undefined variables, and pipe failures
set -euo pipefail

echo "=========================================="
echo "Protein Database Download Job Started"
echo "=========================================="
echo "Date: $(date)"
echo "Hostname: $(hostname)"
echo ""

# Initialize conda
echo "Initializing conda..."
source ~/miniconda3/etc/profile.d/conda.sh

# Activate rnaseq environment
echo "Activating rnaseq conda environment..."
conda activate rnaseq

# Verify environment is activated
echo "Active environment: $CONDA_DEFAULT_ENV"
echo "Conda prefix: $CONDA_PREFIX"
echo "PATH: $PATH" | grep -o "$CONDA_PREFIX[^:]*"

# Set database directory (same as BLAST databases)
DB_DIR=/projects/tholl_lab_1/daisy_analysis/07_NRdatabase
echo "Database directory: $DB_DIR"
echo ""

# Create and navigate to database directory
mkdir -p "$DB_DIR"
cd "$DB_DIR"
echo "Current directory: $(pwd)"
echo ""

# Print diagnostic information
echo "=========================================="
echo "Environment Check"
echo "=========================================="

# Check for required tools
tools_missing=false

echo "Checking for required tools:"
if command -v hmmscan &> /dev/null; then
    echo "  ✓ hmmscan: $(which hmmscan)"
else
    echo "  ✗ hmmscan: NOT FOUND (required for HMMER)"
    tools_missing=true
fi

if command -v hmmpress &> /dev/null; then
    echo "  ✓ hmmpress: $(which hmmpress)"
else
    echo "  ✗ hmmpress: NOT FOUND (required for indexing Pfam)"
    tools_missing=true
fi

if command -v ps_scan.pl &> /dev/null; then
    echo "  ✓ ps_scan.pl: $(which ps_scan.pl)"
else
    echo "  ⚠ ps_scan.pl: NOT FOUND (optional for PROSITE)"
fi

echo ""

if $tools_missing; then
    echo "=========================================="
    echo "ERROR: Required tools not installed"
    echo "=========================================="
    echo ""
    echo "Please update your conda environment to include HMMER:"
    echo ""
    echo "  conda activate rnaseq"
    echo "  conda install -c bioconda hmmer emboss"
    echo ""
    echo "Or update from environment.yml:"
    echo "  conda env update -f environment.yml --prune"
    echo ""
    echo "=========================================="
    exit 1
fi

# ══════════════════════════════════════════════════════════════════════════════
# DOWNLOAD AND PREPARE PFAM DATABASE
# ══════════════════════════════════════════════════════════════════════════════

echo "=========================================="
echo "Processing Pfam Database for HMMER"
echo "=========================================="
echo ""

PFAM_HMM="Pfam-A.hmm"

# Check if Pfam is already prepared (indexed)
if [ -f "${PFAM_HMM}" ] && [ -f "${PFAM_HMM}.h3m" ] && [ -f "${PFAM_HMM}.h3i" ] && [ -f "${PFAM_HMM}.h3f" ] && [ -f "${PFAM_HMM}.h3p" ]; then
    echo "✓ Pfam database already prepared and indexed"
    echo "  Location: ${DB_DIR}/${PFAM_HMM}"
else
    echo "Preparing Pfam database..."
    echo "Start time: $(date)"
    echo ""
    
    # Download if not present
    if [ ! -f "${PFAM_HMM}" ]; then
        echo "Step 1: Downloading Pfam-A.hmm.gz (~600 MB)..."
        echo "Source: http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/"
        
        # Try wget first
        download_success=false
        
        if command -v wget &> /dev/null; then
            echo "Trying wget..."
            if wget --timeout=60 --tries=3 http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz; then
                echo "✓ Download successful with wget"
                download_success=true
            else
                echo "✗ wget failed"
            fi
        fi
        
        # Try curl if wget failed
        if ! $download_success && command -v curl &> /dev/null; then
            echo "Trying curl..."
            if curl --connect-timeout 60 --max-time 3600 -o Pfam-A.hmm.gz http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz; then
                echo "✓ Download successful with curl"
                download_success=true
            else
                echo "✗ curl failed"
            fi
        fi
        
        if ! $download_success; then
            echo "✗ Failed to download Pfam database"
            echo "Please download manually:"
            echo "  wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
            exit 1
        fi
        
        # Verify download
        filesize=$(stat -f%z Pfam-A.hmm.gz 2>/dev/null || stat -c%s Pfam-A.hmm.gz 2>/dev/null)
        echo "Downloaded file size: $((filesize / 1024 / 1024)) MB"
        
        if [ "$filesize" -lt 100000000 ]; then
            echo "✗ Downloaded file too small, likely corrupted"
            rm -f Pfam-A.hmm.gz
            exit 1
        fi
        
        echo ""
        echo "Step 2: Extracting Pfam-A.hmm.gz..."
        if gunzip -v Pfam-A.hmm.gz; then
            echo "✓ Extraction successful"
        else
            echo "✗ Extraction failed"
            exit 1
        fi
    else
        echo "Pfam-A.hmm already exists, skipping download"
    fi
    
    echo ""
    echo "Step 3: Indexing Pfam database with hmmpress..."
    echo "This creates .h3m, .h3i, .h3f, .h3p index files (takes ~5-10 minutes)"
    
    if command -v hmmpress &> /dev/null; then
        if hmmpress -f "$PFAM_HMM"; then
            echo "✓ Pfam database indexed successfully"
        else
            echo "✗ hmmpress failed"
            exit 1
        fi
    else
        echo "✗ hmmpress not found in PATH"
        echo "Make sure HMMER is installed in your conda environment:"
        echo "  conda install -c bioconda hmmer"
        exit 1
    fi
    
    echo ""
    echo "Pfam preparation completed at: $(date)"
fi

echo ""

# ══════════════════════════════════════════════════════════════════════════════
# DOWNLOAD AND PREPARE PROSITE DATABASE
# ══════════════════════════════════════════════════════════════════════════════

echo "=========================================="
echo "Processing PROSITE Database"
echo "=========================================="
echo ""

PROSITE_DIR="${DB_DIR}/prosite"
PROSITE_DAT="prosite.dat"
PROSITE_PSA="prosite.psa"

mkdir -p "$PROSITE_DIR"

# Check if PROSITE is already prepared
if [ -f "${PROSITE_DIR}/${PROSITE_DAT}" ] && [ -f "${PROSITE_DIR}/${PROSITE_PSA}" ]; then
    echo "✓ PROSITE database already prepared"
    echo "  Location: ${PROSITE_DIR}"
else
    echo "Preparing PROSITE database..."
    echo "Start time: $(date)"
    echo ""
    
    cd "$PROSITE_DIR"
    
    # Download prosite.dat
    if [ ! -f "${PROSITE_DAT}" ]; then
        echo "Step 1: Downloading prosite.dat (~30 MB)..."
        echo "Source: ftp://ftp.expasy.org/databases/prosite/"
        
        download_success=false
        
        if command -v wget &> /dev/null; then
            echo "Trying wget..."
            if wget --timeout=60 --tries=3 ftp://ftp.expasy.org/databases/prosite/prosite.dat; then
                echo "✓ Download successful with wget"
                download_success=true
            else
                echo "⚠ wget failed, trying alternative URL..."
            fi
        fi
        
        # Try alternative URL with curl
        if ! $download_success && command -v curl &> /dev/null; then
            echo "Trying curl with alternative source..."
            if curl --connect-timeout 60 --max-time 600 -o prosite.dat https://ftp.expasy.org/databases/prosite/prosite.dat; then
                echo "✓ Download successful with curl"
                download_success=true
            else
                echo "✗ curl failed"
            fi
        fi
        
        if ! $download_success; then
            echo "⚠ Failed to download PROSITE database"
            echo "You can download manually later:"
            echo "  cd ${PROSITE_DIR}"
            echo "  wget ftp://ftp.expasy.org/databases/prosite/prosite.dat"
            echo "  wget ftp://ftp.expasy.org/databases/prosite/prosite.psa"
            cd "$DB_DIR"
            echo ""
            echo "Continuing without PROSITE (Pfam is more commonly used anyway)..."
        else
            # Verify download
            filesize=$(stat -f%z prosite.dat 2>/dev/null || stat -c%s prosite.dat 2>/dev/null)
            echo "Downloaded file size: $((filesize / 1024 / 1024)) MB"
            
            echo ""
            echo "Step 2: Downloading prosite.psa (patterns file)..."
            if command -v wget &> /dev/null; then
                wget --timeout=60 --tries=3 ftp://ftp.expasy.org/databases/prosite/prosite.psa || \
                    curl -o prosite.psa https://ftp.expasy.org/databases/prosite/prosite.psa
            else
                curl -o prosite.psa https://ftp.expasy.org/databases/prosite/prosite.psa
            fi
            
            echo "✓ PROSITE database downloaded"
            echo ""
            echo "PROSITE preparation completed at: $(date)"
            
            cd "$DB_DIR"
        fi
    else
        echo "PROSITE database already exists"
        cd "$DB_DIR"
    fi
fi

echo ""

# ══════════════════════════════════════════════════════════════════════════════
# VERIFY DATABASES
# ══════════════════════════════════════════════════════════════════════════════

echo "=========================================="
echo "Verifying Databases"
echo "=========================================="
echo ""

cd "$DB_DIR"

# Verify Pfam
echo "Checking Pfam database:"
if [ -f "${PFAM_HMM}" ]; then
    filesize=$(stat -f%z "${PFAM_HMM}" 2>/dev/null || stat -c%s "${PFAM_HMM}" 2>/dev/null)
    echo "  ✓ ${PFAM_HMM} exists ($((filesize / 1024 / 1024)) MB)"
    
    # Check index files
    if [ -f "${PFAM_HMM}.h3m" ] && [ -f "${PFAM_HMM}.h3i" ] && [ -f "${PFAM_HMM}.h3f" ] && [ -f "${PFAM_HMM}.h3p" ]; then
        echo "  ✓ Index files (.h3m, .h3i, .h3f, .h3p) present"
        
        # Try to get database stats
        if command -v hmmstat &> /dev/null; then
            echo ""
            echo "  Pfam database statistics:"
            hmmstat "${PFAM_HMM}" | head -20
        fi
    else
        echo "  ✗ Index files missing - run: hmmpress ${PFAM_HMM}"
    fi
else
    echo "  ✗ ${PFAM_HMM} not found"
fi

echo ""

# Verify PROSITE
echo "Checking PROSITE database:"
if [ -d "${PROSITE_DIR}" ] && [ -f "${PROSITE_DIR}/${PROSITE_DAT}" ]; then
    filesize=$(stat -f%z "${PROSITE_DIR}/${PROSITE_DAT}" 2>/dev/null || stat -c%s "${PROSITE_DIR}/${PROSITE_DAT}" 2>/dev/null)
    echo "  ✓ ${PROSITE_DAT} exists ($((filesize / 1024 / 1024)) MB)"
    
    if [ -f "${PROSITE_DIR}/${PROSITE_PSA}" ]; then
        echo "  ✓ ${PROSITE_PSA} (patterns) exists"
    else
        echo "  ⚠ ${PROSITE_PSA} missing"
    fi
    
    # Count entries
    n_entries=$(grep -c "^ID   " "${PROSITE_DIR}/${PROSITE_DAT}" || echo "0")
    echo "  Total PROSITE entries: $n_entries"
else
    echo "  ✗ PROSITE database not found in ${PROSITE_DIR}"
    echo "  (This is optional - Pfam is more commonly used)"
fi

echo ""

# ══════════════════════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════════════════════

echo "=========================================="
echo "SUMMARY"
echo "=========================================="
echo ""
echo "Database directory: ${DB_DIR}"
echo ""
echo "Available databases:"
echo ""

if [ -f "${DB_DIR}/${PFAM_HMM}" ] && [ -f "${DB_DIR}/${PFAM_HMM}.h3m" ]; then
    echo "  ✓ Pfam-A (for HMMER)"
    echo "    Location: ${DB_DIR}/${PFAM_HMM}"
    echo "    Ready for: hmmscan --cpu 16 ${DB_DIR}/${PFAM_HMM} protein.fasta"
else
    echo "  ✗ Pfam-A - NOT READY"
fi

echo ""

if [ -f "${PROSITE_DIR}/${PROSITE_DAT}" ]; then
    echo "  ✓ PROSITE"
    echo "    Location: ${PROSITE_DIR}"
    echo "    Ready for: ps_scan.pl -d ${PROSITE_DIR} protein.fasta"
else
    echo "  ⚠ PROSITE - Not available (optional)"
fi

echo ""
echo "Disk usage:"
du -sh "${DB_DIR}/Pfam-A.hmm"* 2>/dev/null || echo "  Pfam: Not found"
du -sh "${PROSITE_DIR}" 2>/dev/null || echo "  PROSITE: Not found"
echo ""

echo "=========================================="
echo "Next Steps"
echo "=========================================="
echo ""
echo "1. Run HMMER domain scan on your proteins:"
echo "   sbatch scripts/run_hmmer.sbatch DC all_genes_protein.fasta"
echo ""
echo "2. (Optional) For PROSITE motif scanning:"
echo "   # Install PROSITE tools if needed:"
echo "   # conda install -c bioconda emboss"
echo "   # Then create a run_prosite.sbatch script similar to run_hmmer.sbatch"
echo ""

echo "=========================================="
echo "Protein Database Download Job Completed"
echo "Date: $(date)"
echo "=========================================="
