#!/bin/bash
# ---------------------------------------------------------------------
# PURPOSE: Run Trinity de novo transcriptome assembly for ONE sample
# USAGE:   sbatch run_trinity.sbatch <sample_dir> <sample_prefix>
#          where sample_dir is: 00_1_DC, 00_2_DG, 00_3_MF
#          and sample_prefix is: DC1L1, DC1L2, DC1R1, etc.
# EXAMPLE: sbatch run_trinity.sbatch 00_1_DC DC1L1
#
# RESUME:  If a previous job timed out, simply resubmit - the script
#          will automatically detect the incomplete assembly and resume
#          from the checkpoint. Completed assemblies are skipped.
# ---------------------------------------------------------------------

# ===== SLURM directives =====
#SBATCH --job-name=trinity
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G
#SBATCH --time=72:00:00
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
THREADS="${SLURM_CPUS_PER_TASK:-24}"
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

# ===== Activate Conda Environment with Auto-Fix =====
# Source the automated setup script that checks and fixes dependencies
source "${CODE_DIR}/scripts/setup_conda_env.sh"
setup_rnaseq_env || {
    echo "ERROR: Failed to set up conda environment"
    exit 1
}

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

# ===== Check for existing/incomplete assembly =====
RESUME_FLAG=""
if [[ -d "${OUTPUT_DIR}" ]]; then
    if [[ -f "${OUTPUT_DIR}/Trinity.fasta" ]]; then
        echo "============================================================"
        echo "✓ Assembly already complete (Trinity.fasta exists)"
        echo "  Output: ${OUTPUT_DIR}/Trinity.fasta"
        echo "  Exiting successfully - no work needed"
        echo "============================================================"
        echo "To force reassembly, remove the directory first:"
        echo "  rm -rf ${OUTPUT_DIR}"
        exit 0
    else
        echo "============================================================"
        echo "⚠ INCOMPLETE ASSEMBLY DETECTED"
        echo "  Output directory exists but Trinity.fasta is missing"
        echo "  This likely means a previous job timed out"
        echo "  Will RESUME from checkpoint..."
        echo "============================================================"
        RESUME_FLAG="--resume"
    fi
fi

# ===== Verify Critical Dependencies =====
echo "Verifying Trinity dependencies..."

# Check for bowtie2-build-s (required for RSEM abundance estimation)
if ! command -v bowtie2-build-s &> /dev/null; then
    echo "⚠ WARNING: bowtie2-build-s not found!"
    echo "  This is required for abundance estimation via RSEM"
    echo "  Expected location: ${CONDA_PREFIX}/bin/bowtie2-build-s"
    echo ""
    echo "  To fix this, run on login node:"
    echo "    bash scripts/fix_trinity_dependencies.sh"
    echo ""
    echo "  Attempting to continue anyway..."
    echo "  (Trinity may fail during the abundance estimation step)"
fi

# ===== Build Trinity Command =====
CMD="python -m rna_pipeline.cli \
  --mode trinity \
  --reads-left \"${LEFT_FILE}\" \
  --reads-right \"${RIGHT_FILE}\" \
  --outdir \"${OUTPUT_DIR}\" \
  --threads ${THREADS} \
  --mem-gb ${MEM_GB} \
  ${RESUME_FLAG}"

echo "Running command:"
echo "${CMD}"
echo "------------------------------------------------------------"

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

