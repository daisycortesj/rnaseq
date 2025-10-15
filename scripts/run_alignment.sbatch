#!/bin/bash
# ---------------------------------------------------------------------
# PURPOSE: Align RNA-seq reads using STAR
# USAGE:   sbatch run_alignment.sbatch <genome_type> <reads_R1> <reads_R2> [sample_name]
#          where genome_type is: nutmeg, carrot
# EXAMPLES:
#          sbatch run_alignment.sbatch carrot sample_R1.fastq.gz sample_R2.fastq.gz
#          sbatch run_alignment.sbatch nutmeg reads_R1.fq.gz reads_R2.fq.gz mysample
# ---------------------------------------------------------------------

# ===== SLURM directives =====
#SBATCH --job-name=star_align
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=6:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daisycortesj@vt.edu
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/02_mapped
#SBATCH -o alignment_%x_%j.out
#SBATCH -e alignment_%x_%j.err

set -euo pipefail

# ===== Configuration =====
GENOME_TYPE="${1:-}"
READS_R1="${2:-}"
READS_R2="${3:-}"
SAMPLE_NAME="${4:-}"

REFERENCE_DIR="/projects/tholl_lab_1/daisy_analysis/04_reference"
CODE_DIR="/projects/tholl_lab_1/daisy_analysis/05_rnaseq-code"
ALIGNMENT_DIR="/projects/tholl_lab_1/daisy_analysis/02_mapped"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

echo "============================================================"
echo "STAR Read Alignment Pipeline"
echo "============================================================"
echo "SLURM job running on host: $(hostname)"
echo "CPUs: ${THREADS}  MEM: ${SLURM_MEM_PER_NODE:-64G}"
echo "Genome type: ${GENOME_TYPE}"
echo "Start time: $(date)"
echo "============================================================"

# ===== Validate Arguments =====
if [[ -z "${GENOME_TYPE}" ]]; then
    echo "ERROR: No genome type provided"
    echo "Usage: sbatch run_alignment.sbatch <genome_type> <reads_R1> <reads_R2> [sample_name]"
    echo "Genome types: nutmeg, carrot"
    exit 1
fi

if [[ -z "${READS_R1}" ]] || [[ -z "${READS_R2}" ]]; then
    echo "ERROR: Read files not provided"
    echo "Usage: sbatch run_alignment.sbatch <genome_type> <reads_R1> <reads_R2> [sample_name]"
    exit 1
fi

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
STAR --version || echo "STAR not found."
echo "------------------------------------------------------------"

# ===== Auto-detect Sample Directory Name =====
# Extract sample identifier from reads path (e.g., 00_1_DC, 00_2_DG, 00_3_MF)
READS_DIR=$(dirname "${READS_R1}")
SAMPLE_DIR=$(basename "${READS_DIR}")

echo "Detected sample directory: ${SAMPLE_DIR}"

# ===== Set Genome-Specific Parameters =====
case "${GENOME_TYPE}" in
    nutmeg)
        echo "Aligning to NUTMEG genome (MYU assembly)"
        GENOME_INDEX="${REFERENCE_DIR}/nutmeg_star_index"
        QUANT_MODE="--no-quant-mode"  # No GTF for nutmeg
        OUTDIR="${ALIGNMENT_DIR}/${SAMPLE_DIR}"
        ;;
    
    carrot)
        echo "Aligning to CARROT genome (with gene quantification)"
        GENOME_INDEX="${REFERENCE_DIR}/carrot_star_index"
        QUANT_MODE="--quant-mode"  # Has GTF
        OUTDIR="${ALIGNMENT_DIR}/${SAMPLE_DIR}"
        ;;
    
    *)
        echo "ERROR: Unknown genome type '${GENOME_TYPE}'"
        echo "Valid options: nutmeg, carrot"
        exit 1
        ;;
esac

# ===== Validate Input Files =====
echo "Validating input files..."

if [[ ! -d "${GENOME_INDEX}" ]]; then
    echo "ERROR: Genome index not found: ${GENOME_INDEX}"
    echo "Please build the index first using run_genome_index.sbatch"
    exit 1
fi
echo "✓ Genome index: ${GENOME_INDEX}"

if [[ ! -f "${READS_R1}" ]]; then
    echo "ERROR: R1 reads file not found: ${READS_R1}"
    exit 1
fi
echo "✓ R1 reads: ${READS_R1}"

if [[ ! -f "${READS_R2}" ]]; then
    echo "ERROR: R2 reads file not found: ${READS_R2}"
    exit 1
fi
echo "✓ R2 reads: ${READS_R2}"

echo "✓ Output directory: ${OUTDIR}"
echo "------------------------------------------------------------"

# ===== Go to Pipeline Directory =====
cd "${CODE_DIR}"

# ===== Build STAR Alignment Command =====
if [[ -n "${SAMPLE_NAME}" ]]; then
    SAMPLE_ARG="--sample-name ${SAMPLE_NAME}"
else
    SAMPLE_ARG=""
fi

CMD="python -m rna_pipeline.cli \
  --mode align \
  --genome-index ${GENOME_INDEX} \
  --reads-left ${READS_R1} \
  --reads-right ${READS_R2} \
  --outdir ${OUTDIR} \
  --threads ${THREADS} \
  ${QUANT_MODE} \
  ${SAMPLE_ARG}"

echo "Running command:"
echo "${CMD}"
echo "------------------------------------------------------------"

# ===== Run Pipeline =====
eval ${CMD}

# ===== Summary =====
echo "============================================================"
echo "Alignment completed successfully!"
echo "End time: $(date)"
echo "Output directory: ${OUTDIR}"
echo "============================================================"

# List output files
if [[ -d "${OUTDIR}" ]]; then
    echo "Generated files:"
    ls -lh "${OUTDIR}"
fi

