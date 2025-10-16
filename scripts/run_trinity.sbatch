#!/bin/bash
# ---------------------------------------------------------------------
# PURPOSE: Run Trinity de novo transcriptome assembly for ONE sample
# USAGE:   sbatch run_trinity.sbatch <sample_dir> <sample_prefix>
#          where sample_dir is: 00_1_DC, 00_2_DG, 00_3_MF
#          and sample_prefix is: DC1L1, DC1L2, DC1R1, etc.
# EXAMPLE: sbatch run_trinity.sbatch 00_1_DC DC1L1
# ---------------------------------------------------------------------

# ===== SLURM directives =====
#SBATCH --job-name=trinity
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daisycortesj@vt.edu
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/01_processed
#SBATCH -o trinity_%x_%j.out
#SBATCH -e trinity_%x_%j.err

set -euo pipefail

# ===== Configuration =====
SAMPLE_DIR="${1:-}"
SAMPLE_PREFIX="${2:-}"
RAWDATA_DIR="/projects/tholl_lab_1/daisy_analysis/00_rawdata"
CODE_DIR="/projects/tholl_lab_1/daisy_analysis/05_rnaseq-code"
TRINITY_OUTDIR="/projects/tholl_lab_1/daisy_analysis/01_processed"
THREADS="${SLURM_CPUS_PER_TASK:-16}"
MEM_GB="128"

echo "============================================================"
echo "Trinity De Novo Transcriptome Assembly - Individual Sample"
echo "============================================================"
echo "SLURM job running on host: $(hostname)"
echo "CPUs: ${THREADS}  MEM: ${SLURM_MEM_PER_NODE:-128G}"
echo "Sample directory: ${SAMPLE_DIR}"
echo "Sample prefix: ${SAMPLE_PREFIX}"
echo "Start time: $(date)"
echo "============================================================"

# ===== Validate Arguments =====
if [[ -z "${SAMPLE_DIR}" ]] || [[ -z "${SAMPLE_PREFIX}" ]]; then
    echo "ERROR: Missing required arguments"
    echo "Usage: sbatch run_trinity.sbatch <sample_dir> <sample_prefix>"
    echo "Example: sbatch run_trinity.sbatch 00_1_DC DC1L1"
    echo ""
    echo "Sample directories: 00_1_DC, 00_2_DG, 00_3_MF"
    echo "Sample prefixes in 00_1_DC: DC1L1, DC1L2, DC1L3, DC1R1, DC1R2, DC1R3, DC2L1, DC2L2, DC2L3, DC2R1, DC2R2, DC2R3"
    exit 1
fi

# ===== Validate Sample Directory =====
FULL_SAMPLE_PATH="${RAWDATA_DIR}/${SAMPLE_DIR}"
if [[ ! -d "${FULL_SAMPLE_PATH}" ]]; then
    echo "ERROR: Sample directory not found: ${FULL_SAMPLE_PATH}"
    exit 1
fi

echo "✓ Sample directory found: ${FULL_SAMPLE_PATH}"

# ===== Activate Conda Environment =====
if ! bash -c "source ~/.bashrc" 2>/dev/null; then
    echo "Warning: .bashrc has issues, using alternative conda activation"
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate rnaseq
else
    echo "Loading .bashrc successfully"
    source ~/.bashrc
    conda activate rnaseq
fi

# Verify conda activation
conda info --envs | grep rnaseq || echo "ERROR: conda environment not activated"

# Check tools
echo "Checking tool versions:"
Trinity --version || echo "Trinity not found."
echo "------------------------------------------------------------"

# ===== Collect Read Files for This Sample =====
echo "Searching for ${SAMPLE_PREFIX} reads in ${FULL_SAMPLE_PATH}..."

# Find the specific sample files
LEFT_FILE=$(find "${FULL_SAMPLE_PATH}" -name "${SAMPLE_PREFIX}_1.fq.gz" | head -1)
RIGHT_FILE=$(find "${FULL_SAMPLE_PATH}" -name "${SAMPLE_PREFIX}_2.fq.gz" | head -1)

if [[ -z "${LEFT_FILE}" ]] || [[ -z "${RIGHT_FILE}" ]]; then
    echo "ERROR: Could not find read files for sample ${SAMPLE_PREFIX}"
    echo "Expected: ${SAMPLE_PREFIX}_1.fq.gz and ${SAMPLE_PREFIX}_2.fq.gz"
    echo ""
    echo "Available samples in ${SAMPLE_DIR}:"
    find "${FULL_SAMPLE_PATH}" -name "*_1.fq.gz" | xargs -n1 basename | sed 's/_1.fq.gz$//' | sort
    exit 1
fi

# Verify files are readable
if [[ ! -r "${LEFT_FILE}" ]]; then
    echo "ERROR: Cannot read file: ${LEFT_FILE}"
    exit 1
fi

if [[ ! -r "${RIGHT_FILE}" ]]; then
    echo "ERROR: Cannot read file: ${RIGHT_FILE}"
    exit 1
fi

echo "✓ Found read files:"
echo "  Left:  ${LEFT_FILE}"
echo "  Right: ${RIGHT_FILE}"

# Check if symlinks and show actual file location
if [[ -L "${LEFT_FILE}" ]]; then
    echo "  (symlink to: $(readlink -f "${LEFT_FILE}"))"
fi

echo "------------------------------------------------------------"

# ===== Set Output Directory =====
OUTPUT_DIR="${TRINITY_OUTDIR}/${SAMPLE_PREFIX}_trinity_assembly"
echo "✓ Output directory: ${OUTPUT_DIR}"

# Create output directory if it doesn't exist
mkdir -p "${TRINITY_OUTDIR}"

# ===== Go to Pipeline Directory =====
cd "${CODE_DIR}"

# ===== Build Trinity Command =====
CMD="python -m rna_pipeline.cli \
  --mode trinity \
  --reads-left \"${LEFT_FILE}\" \
  --reads-right \"${RIGHT_FILE}\" \
  --outdir \"${OUTPUT_DIR}\" \
  --threads ${THREADS} \
  --mem-gb ${MEM_GB}"

echo "Running command:"
echo "${CMD}"
echo "------------------------------------------------------------"

# ===== Check if output directory exists =====
if [[ -d "${OUTPUT_DIR}" ]]; then
    echo "WARNING: Output directory already exists: ${OUTPUT_DIR}"
    echo "Trinity will fail if this directory contains a previous run."
    echo "To restart, remove the directory first: rm -rf ${OUTPUT_DIR}"
    echo "Attempting to run anyway..."
fi

# ===== Run Pipeline =====
eval ${CMD}

# ===== Summary =====
echo "============================================================"
echo "Trinity assembly completed successfully!"
echo "End time: $(date)"
echo "Output directory: ${OUTPUT_DIR}"
echo "============================================================"

# Check for Trinity output
if [[ -f "${OUTPUT_DIR}/Trinity.fasta" ]]; then
    echo "✓ Trinity assembly file: ${OUTPUT_DIR}/Trinity.fasta"
    echo "Assembly stats:"
    grep -c ">" "${OUTPUT_DIR}/Trinity.fasta" || echo "  Contigs: (count not available)"
else
    echo "Warning: Trinity.fasta not found in expected location"
    echo "Checking output directory contents:"
    ls -lh "${OUTPUT_DIR}/" 2>/dev/null || echo "Directory not accessible"
fi

