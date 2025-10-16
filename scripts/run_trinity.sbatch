#!/bin/bash
# ---------------------------------------------------------------------
# PURPOSE: Run Trinity de novo transcriptome assembly
# USAGE:   sbatch run_trinity.sbatch <sample_dir>
#          where sample_dir is: 00_1_DC, 00_2_DG, 00_3_MF
# EXAMPLE: sbatch run_trinity.sbatch 00_1_DC
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
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis
#SBATCH -o trinity_%x_%j.out
#SBATCH -e trinity_%x_%j.err

set -euo pipefail

# ===== Configuration =====
SAMPLE_DIR="${1:-}"
RAWDATA_DIR="/projects/tholl_lab_1/daisy_analysis/00_rawdata"
CODE_DIR="/projects/tholl_lab_1/daisy_analysis/05_rnaseq-code"
TRINITY_OUTDIR="/projects/tholl_lab_1/daisy_analysis/03_trinity"
THREADS="${SLURM_CPUS_PER_TASK:-16}"
MEM_GB="128"

echo "============================================================"
echo "Trinity De Novo Transcriptome Assembly"
echo "============================================================"
echo "SLURM job running on host: $(hostname)"
echo "CPUs: ${THREADS}  MEM: ${SLURM_MEM_PER_NODE:-128G}"
echo "Sample directory: ${SAMPLE_DIR}"
echo "Start time: $(date)"
echo "============================================================"

# ===== Validate Arguments =====
if [[ -z "${SAMPLE_DIR}" ]]; then
    echo "ERROR: No sample directory provided"
    echo "Usage: sbatch run_trinity.sbatch <sample_dir>"
    echo "Sample directories: 00_1_DC, 00_2_DG, 00_3_MF"
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

# ===== Collect Read Files =====
echo "Collecting read files from ${SAMPLE_DIR}..."

# Find all _1.fq.gz files (left reads)
LEFT_FILES=$(find "${FULL_SAMPLE_PATH}" -name "*_1.fq.gz" | sort | tr '\n' ',' | sed 's/,$//')
# Find all _2.fq.gz files (right reads)
RIGHT_FILES=$(find "${FULL_SAMPLE_PATH}" -name "*_2.fq.gz" | sort | tr '\n' ',' | sed 's/,$//')

if [[ -z "${LEFT_FILES}" ]] || [[ -z "${RIGHT_FILES}" ]]; then
    echo "ERROR: No paired-end read files found in ${FULL_SAMPLE_PATH}"
    exit 1
fi

echo "Left reads files:"
echo "${LEFT_FILES}" | tr ',' '\n'
echo ""
echo "Right reads files:"
echo "${RIGHT_FILES}" | tr ',' '\n'
echo "------------------------------------------------------------"

# ===== Set Output Directory =====
OUTPUT_DIR="${TRINITY_OUTDIR}/${SAMPLE_DIR}_assembly"
echo "✓ Output directory: ${OUTPUT_DIR}"

# Create output directory if it doesn't exist
mkdir -p "${TRINITY_OUTDIR}"

# ===== Go to Pipeline Directory =====
cd "${CODE_DIR}"

# ===== Build Trinity Command =====
CMD="python -m rna_pipeline.cli \
  --mode trinity \
  --reads-left ${LEFT_FILES} \
  --reads-right ${RIGHT_FILES} \
  --outdir ${OUTPUT_DIR} \
  --threads ${THREADS} \
  --mem-gb ${MEM_GB}"

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

