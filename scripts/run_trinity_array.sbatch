#!/bin/bash
# ---------------------------------------------------------------------
# PURPOSE: Run Trinity de novo assembly as a SLURM job array
#          Each sample runs as an independent job with its own time limit
# USAGE:   
#   1. First, generate the sample list:
#      ls /projects/tholl_lab_1/daisy_analysis/00_rawdata/00_*/*_1.fq.gz > sample_list.txt
#   2. Count samples:
#      wc -l sample_list.txt
#   3. Submit the array (replace N with sample count):
#      sbatch --array=1-N run_trinity_array.sbatch
#
# RESUME:  If jobs time out, simply resubmit - incomplete assemblies
#          will automatically resume from checkpoint.
#
# ADVANTAGES over run_trinity_all.sbatch:
#   - Each sample has its own 72-hour time limit
#   - Samples can run in parallel (faster overall)
#   - If one sample fails/times out, others are unaffected
# ---------------------------------------------------------------------

# ===== SLURM directives =====
#SBATCH --job-name=trinity_arr
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/01_processed
#SBATCH -o trinity_array_%A_%a.out
#SBATCH -e trinity_array_%A_%a.err

set -euo pipefail

# ===== Configuration =====
SAMPLE_LIST="/projects/tholl_lab_1/daisy_analysis/01_processed/sample_list.txt"
CODE_DIR="/projects/tholl_lab_1/daisy_analysis/05_rnaseq-code"
TRINITY_OUTDIR="/projects/tholl_lab_1/daisy_analysis/01_processed"
THREADS="${SLURM_CPUS_PER_TASK:-24}"
MEM_GB="128"

echo "============================================================"
echo "Trinity De Novo Assembly - Job Array"
echo "============================================================"
echo "Array Job ID: ${SLURM_ARRAY_JOB_ID:-N/A}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID:-N/A}"
echo "SLURM job running on host: $(hostname)"
echo "CPUs: ${THREADS}  MEM: ${SLURM_MEM_PER_NODE:-128G}"
echo "Start time: $(date)"
echo "============================================================"

# ===== Validate sample list =====
if [[ ! -f "${SAMPLE_LIST}" ]]; then
    echo "ERROR: Sample list not found: ${SAMPLE_LIST}"
    echo ""
    echo "Generate it with:"
    echo "  ls /projects/tholl_lab_1/daisy_analysis/00_rawdata/00_*/*_1.fq.gz > ${SAMPLE_LIST}"
    exit 1
fi

TOTAL_SAMPLES=$(wc -l < "${SAMPLE_LIST}")
echo "Total samples in list: ${TOTAL_SAMPLES}"

if [[ ${SLURM_ARRAY_TASK_ID:-1} -gt ${TOTAL_SAMPLES} ]]; then
    echo "ERROR: Task ID ${SLURM_ARRAY_TASK_ID} exceeds sample count ${TOTAL_SAMPLES}"
    exit 1
fi

# ===== Get sample for this task =====
LEFT_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID:-1}p" "${SAMPLE_LIST}")
base=$(basename "${LEFT_FILE}" _1.fq.gz)
sample_dir=$(dirname "${LEFT_FILE}")
RIGHT_FILE="${sample_dir}/${base}_2.fq.gz"

echo ""
echo "Processing sample ${SLURM_ARRAY_TASK_ID:-1} of ${TOTAL_SAMPLES}: ${base}"
echo "  Left:  ${LEFT_FILE}"
echo "  Right: ${RIGHT_FILE}"

# ===== Validate files =====
if [[ ! -r "${LEFT_FILE}" ]]; then
    echo "ERROR: Cannot read left file: ${LEFT_FILE}"
    exit 1
fi

if [[ ! -r "${RIGHT_FILE}" ]]; then
    echo "ERROR: Cannot read right file: ${RIGHT_FILE}"
    exit 1
fi

echo "  ✓ Both read files found and readable"

# ===== Activate Conda Environment =====
if ! bash -c "source ~/.bashrc" 2>/dev/null; then
    echo "Warning: .bashrc has issues, using alternative conda activation"
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate rnaseq
else
    source ~/.bashrc
    conda activate rnaseq
fi

# Verify conda activation
conda info --envs | grep rnaseq || { echo "ERROR: conda environment not activated"; exit 1; }

# ===== Go to Pipeline Directory =====
cd "${CODE_DIR}"

# ===== Set output directory =====
OUTPUT_DIR="${TRINITY_OUTDIR}/${base}_trinity_assembly"
echo "  Output: ${OUTPUT_DIR}"

# ===== Check for existing/incomplete assembly =====
RESUME_FLAG=""
if [[ -d "${OUTPUT_DIR}" ]]; then
    if [[ -f "${OUTPUT_DIR}/Trinity.fasta" ]]; then
        echo "  ✓ Assembly already complete (Trinity.fasta exists)"
        echo "  Exiting successfully - no work needed"
        exit 0
    else
        echo "  ⚠ Incomplete assembly detected - will RESUME from checkpoint"
        RESUME_FLAG="--resume"
    fi
fi

# ===== Run Trinity =====
echo ""
echo "Running Trinity assembly..."
CMD="python -m rna_pipeline.cli \
  --mode trinity \
  --reads-left \"${LEFT_FILE}\" \
  --reads-right \"${RIGHT_FILE}\" \
  --outdir \"${OUTPUT_DIR}\" \
  --threads ${THREADS} \
  --mem-gb ${MEM_GB} \
  ${RESUME_FLAG}"

echo "Command: ${CMD}"
echo "------------------------------------------------------------"

if eval ${CMD} 2>&1; then
    echo ""
    echo "✓ ${base} completed successfully"
    
    # Check for Trinity output
    if [[ -f "${OUTPUT_DIR}/Trinity.fasta" ]]; then
        CONTIG_COUNT=$(grep -c ">" "${OUTPUT_DIR}/Trinity.fasta" || echo "0")
        echo "✓ Trinity.fasta: ${CONTIG_COUNT} contigs"
    else
        echo "⚠ Warning: Trinity.fasta not found"
    fi
else
    EXIT_CODE=$?
    echo ""
    echo "✗ ${base} failed (exit code: ${EXIT_CODE})"
    exit ${EXIT_CODE}
fi

echo ""
echo "============================================================"
echo "End time: $(date)"
echo "============================================================"

exit 0
