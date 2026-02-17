#!/bin/bash
################################################################################
# Download and Extract BLAST Databases
################################################################################
#
# This script downloads and prepares BLAST protein databases:
#   - swissprot (curated, fast, ~1-2 GB)
#   - nr (comprehensive, slow, ~128 GB)
#   - refseq_protein (balanced, ~30 GB)
#   - taxdb (taxonomy for species filtering)
#
# Features:
#   - Automatically extracts existing .tar.gz files before downloading
#   - Skips databases that are already ready
#   - Fallback FTP download for swissprot if update_blastdb.pl fails
#   - Verifies all databases after processing
#
# Usage:
#   sbatch scripts/download_blastdb.sbatch
#
# Runtime: 1-48 hours depending on which databases need downloading
#
################################################################################

# ══════════════════════════════════════════════════════════════════════════════
# SLURM SETTINGS 
# ══════════════════════════════════════════════════════════════════════════════

#SBATCH --job-name=blastdb_download         # Job name
#SBATCH --account=tholl_lab_1             # Your lab account
#SBATCH --partition=normal_q              # Which queue to use
#SBATCH --cpus-per-task=8                 # Number of CPUs (8 for speed)
#SBATCH --mem=32G                         # Memory (16GB for PyDESeq2)
#SBATCH --time=48:00:00                    # Max time (48 hours)
#SBATCH --mail-type=ALL                   # Email when done
#SBATCH --mail-user=daisycortesj@vt.edu   # Your email

# Where to save output and error messages
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/07_NRdatabase
#SBATCH -o blastdb_download_%j.out
#SBATCH -e blastdb_download_%j.err

# Exit on error, undefined variables, and pipe failures
set -euo pipefail

echo "=========================================="
echo "BLAST Database Download Job Started"
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

# Export BLASTDB path
export BLASTDB=/projects/tholl_lab_1/daisy_analysis/07_NRdatabase/blastdb
echo "BLASTDB set to: $BLASTDB"
echo ""

# Create and navigate to BLASTDB directory
echo "Creating and navigating to BLASTDB directory..."
mkdir -p "$BLASTDB"
cd "$BLASTDB"
echo "Current directory: $(pwd)"
echo ""

# Print diagnostic information
echo "=========================================="
echo "Environment Check"
echo "=========================================="
echo "PATH: $PATH"
echo ""
echo "which update_blastdb.pl: $(which update_blastdb.pl)"
echo "which blastdbcmd: $(which blastdbcmd)"
echo "which blastx: $(which blastx)"
echo ""
echo "BLAST version:"
blastx -version
echo ""

# Download and decompress databases
echo "=========================================="
echo "Downloading BLAST Databases"
echo "=========================================="

# Function to extract existing tar.gz files
extract_existing_archives() {
    local db_pattern=$1
    echo "Checking for existing ${db_pattern} archives to extract..."
    
    if ls ${db_pattern}*.tar.gz 1> /dev/null 2>&1; then
        echo "Found ${db_pattern} archives. Extracting..."
        for archive in ${db_pattern}*.tar.gz; do
            if [ -f "$archive" ]; then
                echo "  Extracting: $archive"
                tar -xzf "$archive" || echo "  Warning: Failed to extract $archive"
            fi
        done
        echo "Extraction complete for ${db_pattern}"
    else
        echo "No ${db_pattern} archives found to extract."
    fi
    echo ""
}

# Check and download/extract swissprot
echo "=========================================="
echo "Processing SwissProt Database"
echo "=========================================="
if blastdbcmd -db swissprot -info &>/dev/null; then
    echo "✓ swissprot database already ready"
else
    # Try extracting existing archives first
    extract_existing_archives "swissprot"
    
    # If still not ready, try downloading
    if ! blastdbcmd -db swissprot -info &>/dev/null; then
        echo "Attempting to download swissprot..."
        echo "Start time: $(date)"
        
        # Try update_blastdb.pl first
        if update_blastdb.pl --passive --decompress swissprot 2>/dev/null; then
            echo "✓ swissprot downloaded successfully via update_blastdb.pl"
        else
            echo "⚠ update_blastdb.pl failed. Trying direct FTP download..."
            
            # Fallback: direct download from NCBI FTP
            if command -v wget &> /dev/null; then
                wget -q ftp://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz && \
                tar -xzf swissprot.tar.gz && \
                echo "✓ swissprot downloaded and extracted via FTP" || \
                echo "✗ FTP download also failed"
            elif command -v curl &> /dev/null; then
                curl -s ftp://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz -o swissprot.tar.gz && \
                tar -xzf swissprot.tar.gz && \
                echo "✓ swissprot downloaded and extracted via FTP" || \
                echo "✗ FTP download also failed"
            else
                echo "✗ Neither wget nor curl available. Cannot download swissprot."
            fi
        fi
        echo "swissprot processing completed at: $(date)"
    fi
fi
echo ""

# Check and download nr
echo "=========================================="
echo "Processing NR Database"
echo "=========================================="
if blastdbcmd -db nr -info &>/dev/null; then
    echo "✓ nr database already ready"
else
    # Try extracting existing archives first
    extract_existing_archives "nr"
    
    # If still not ready, download
    if ! blastdbcmd -db nr -info &>/dev/null; then
        echo "Downloading nr..."
        echo "Start time: $(date)"
        update_blastdb.pl --passive --decompress nr
        echo "nr download completed at: $(date)"
    fi
fi
echo ""

# Check and download/extract refseq_protein
echo "=========================================="
echo "Processing RefSeq Protein Database"
echo "=========================================="
if blastdbcmd -db refseq_protein -info &>/dev/null; then
    echo "✓ refseq_protein database already ready"
else
    # Try extracting existing archives first
    extract_existing_archives "refseq_protein"
    
    # If still not ready, download
    if ! blastdbcmd -db refseq_protein -info &>/dev/null; then
        echo "Downloading refseq_protein..."
        echo "Start time: $(date)"
        update_blastdb.pl --passive --decompress refseq_protein
        echo "refseq_protein download completed at: $(date)"
    fi
fi
echo ""

# Download taxdb (taxonomy database - useful for filtering by species)
echo "=========================================="
echo "Processing TaxDB (Taxonomy Database)"
echo "=========================================="
if [ -f "taxdb.btd" ] && [ -f "taxdb.bti" ]; then
    echo "✓ taxdb already exists"
else
    echo "Downloading taxdb..."
    echo "Start time: $(date)"
    update_blastdb.pl --passive --decompress taxdb
    echo "taxdb download completed at: $(date)"
fi
echo ""

# Verify all databases
echo "=========================================="
echo "Verifying Downloaded Databases"
echo "=========================================="

echo "Verifying swissprot database:"
if blastdbcmd -db swissprot -info &>/dev/null; then
    blastdbcmd -db swissprot -info | head -n 10
    echo "✓ swissprot is ready to use"
else
    echo "✗ swissprot verification failed - database may not be usable"
fi
echo ""

echo "Verifying nr database:"
if blastdbcmd -db nr -info &>/dev/null; then
    blastdbcmd -db nr -info | head -n 10
    echo "✓ nr is ready to use"
else
    echo "✗ nr verification failed - database may not be usable"
fi
echo ""

echo "Verifying refseq_protein database:"
if blastdbcmd -db refseq_protein -info &>/dev/null; then
    blastdbcmd -db refseq_protein -info | head -n 10
    echo "✓ refseq_protein is ready to use"
else
    echo "✗ refseq_protein verification failed - database may not be usable"
fi
echo ""

# Summary
echo "=========================================="
echo "SUMMARY"
echo "=========================================="
echo "Available databases:"
ls -lh *.pin 2>/dev/null | awk '{print "  " $9}' | sed 's/.pin$//' | sort -u || echo "  None found"
echo ""
echo "Disk usage:"
du -sh . 2>/dev/null || echo "  Unable to calculate"
echo ""

echo "=========================================="
echo "BLAST Database Download Job Completed"
echo "Date: $(date)"
echo "=========================================="
