#!/bin/bash

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

# Check and download swissprot
if blastdbcmd -db swissprot -info &>/dev/null; then
    echo "swissprot database already exists, skipping download..."
else
    echo "Downloading swissprot..."
    echo "Start time: $(date)"
    update_blastdb.pl --passive --decompress swissprot
    echo "swissprot download completed at: $(date)"
fi
echo ""

# Check and download nr
if blastdbcmd -db nr -info &>/dev/null; then
    echo "nr database already exists, skipping download..."
else
    echo "Downloading nr..."
    echo "Start time: $(date)"
    update_blastdb.pl --passive --decompress nr
    echo "nr download completed at: $(date)"
fi
echo ""

# Uncomment these later when ready to download additional databases:
# echo "Downloading refseq_protein..."
# echo "Start time: $(date)"
# update_blastdb.pl --passive --decompress refseq_protein
# echo "refseq_protein download completed at: $(date)"
# echo ""
# 
# echo "Downloading taxdb..."
# echo "Start time: $(date)"
# update_blastdb.pl --passive --decompress taxdb
# echo "taxdb download completed at: $(date)"
# echo ""

# Verify downloads
echo "=========================================="
echo "Verifying Downloaded Databases"
echo "=========================================="

echo "Verifying swissprot database:"
blastdbcmd -db swissprot -info | head -n 40
echo ""

echo "Verifying nr database:"
blastdbcmd -db nr -info | head -n 40
echo ""

echo "=========================================="
echo "BLAST Database Download Job Completed"
echo "Date: $(date)"
echo "=========================================="
