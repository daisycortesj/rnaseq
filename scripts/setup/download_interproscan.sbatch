#!/bin/bash
# =============================================================================
# SCRIPT:  download_interproscan.sbatch
#
# WHAT IT DOES:
#   Downloads InterProScan ONCE (software + domain databases, ~14 GB).
#   If it is already installed, this script does nothing and exits.
#
#   Same idea as download_blastdb.sbatch: one setup job, then every
#   analysis job reuses the files.
#
# HOW TO SUBMIT (from the repo folder):
#   cd /projects/tholl_lab_1/daisy_analysis/05_rnaseq-code
#   sbatch scripts/setup/download_interproscan.sbatch
#
# Safe to submit again — it will print SKIP and stop.
#
# AFTER THIS SUCCEEDS:
#   sbatch scripts/07_domains/run_interproscan.sbatch
# =============================================================================

#SBATCH --job-name=download_interproscan
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/04_reference
#SBATCH --output=download_interproscan_%j.out
#SBATCH --error=download_interproscan_%j.err

set -euo pipefail


############################################################
## 1. Where to install (must match run_interproscan.sbatch)
############################################################

INSTALL_DIR="/projects/tholl_lab_1/daisy_analysis/04_reference/software"
VERSION="5.75-106.0"
INTERPRO_DIR="$INSTALL_DIR/interproscan-$VERSION"

mkdir -p "$INSTALL_DIR"

echo "========================================================"
echo "InterProScan download"
echo "Install folder: $INTERPRO_DIR"
echo "Host:           $(hostname)"
echo "Started:        $(date)"
echo "========================================================"


############################################################
## 2. Skip if already downloaded
############################################################

if [[ -f "$INTERPRO_DIR/interproscan.sh" ]]; then
    echo "SKIP — InterProScan is already installed:"
    echo "  $INTERPRO_DIR/interproscan.sh"
    echo "Nothing was downloaded."
    echo "Next: sbatch scripts/07_domains/run_interproscan.sbatch"
    exit 0
fi


############################################################
## 3. Download, check, unzip, index
############################################################

cd "$INSTALL_DIR"

echo "Downloading InterProScan $VERSION (this is large, be patient)..."

wget -c \
    "https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/$VERSION/interproscan-$VERSION-64-bit.tar.gz"

wget -c \
    "https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/$VERSION/interproscan-$VERSION-64-bit.tar.gz.md5"

echo "Checking the download is not corrupted..."
md5sum -c "interproscan-$VERSION-64-bit.tar.gz.md5"

echo "Unzipping..."
tar -pxzf "interproscan-$VERSION-64-bit.tar.gz"

echo "Indexing databases (required once after unzipping)..."
cd "$INTERPRO_DIR"
python3 setup.py -f interproscan.properties

chmod +x "$INTERPRO_DIR/interproscan.sh"


############################################################
## 4. Confirm it worked
############################################################

if [[ ! -f "$INTERPRO_DIR/interproscan.sh" ]]; then
    echo "ERROR: Download finished but interproscan.sh is missing."
    echo "Looked in: $INTERPRO_DIR"
    exit 1
fi

echo "========================================================"
echo "DONE — InterProScan is installed"
echo "Path:     $INTERPRO_DIR"
echo "Finished: $(date)"
echo "Next:     sbatch scripts/07_domains/run_interproscan.sbatch"
echo "========================================================"
