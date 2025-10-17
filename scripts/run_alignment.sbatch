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
#SBATCH -o bldstridx_output.txt
#SBATCH -e bldstridx_output.txt

set -euo pipefail

# ===== Configuration =====
GENOME_TYPE="${1:-}"
READS_R1="${2:-}"
READS_R2="${3:-}"
SAMPLE_NAME="${4:-}"

#**The syntax `${1:-}` means:** "Get the first argument, or use an empty string if nothing was provided"


REFERENCE_DIR="/projects/tholl_lab_1/daisy_analysis/04_reference" # where genome lives
CODE_DIR="/projects/tholl_lab_1/daisy_analysis/05_rnaseq-code" # where python code live
ALIGNMENT_DIR="/projects/tholl_lab_1/daisy_analysis/02_mapped" # where results are saved
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
# if genome is not given, don't continue running
if [[ -z "${GENOME_TYPE}" ]]; then
    echo "ERROR: No genome type provided"
    echo "Usage: sbatch run_alignment.sbatch <genome_type> <sample_dir>"
    echo "Genome types: nutmeg, carrot"
    echo "Sample dirs: 00_1_DC, 00_2_DG, 00_3_MF"
    exit 1
fi

if [[ -z "${SAMPLE_DIR}" ]]; then
    echo "ERROR: Sample directory not provided"
    echo "Usage: sbatch run_alignment.sbatch <genome_type> <sample_dir>"
    echo "Example: sbatch run_alignment.sbatch carrot 00_1_DC"
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

# ===== Validate Sample Directory Exists =====
SAMPLE_PATH="${RAWDATA_DIR}/${SAMPLE_DIR}"
if [[ ! -d "${SAMPLE_PATH}" ]]; then
    echo "ERROR: Sample directory not found: ${SAMPLE_PATH}"
    echo "Available directories in ${RAWDATA_DIR}:"
    ls -1 "${RAWDATA_DIR}" 2>/dev/null || echo "(none)"
    exit 1
fi

echo "✓ Sample directory: ${SAMPLE_PATH}"

# ===== Set Genome-Specific Parameters =====
# Depending on which genome was chosen, use different settings
# nutmeg found nutmeg_star_index
# no GTF
case "${GENOME_TYPE}" in
    nutmeg)
        echo "Aligning to NUTMEG genome (MYU assembly)"
        GENOME_INDEX="${REFERENCE_DIR}/nutmeg_star_index"
        QUANT_MODE="--no-quant-mode"  # No GTF for nutmeg
        OUTDIR="${ALIGNMENT_DIR}/${SAMPLE_DIR}"
        ;;
# carrot found in carrot_star_index
# with GTF  
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

# ===== Validate Genome Index =====
echo "Validating genome index..."

if [[ ! -d "${GENOME_INDEX}" ]]; then
    echo "ERROR: Genome index not found: ${GENOME_INDEX}"
    echo "Please build the index first using run_genome_index.sbatch"
    exit 1
fi
echo "✓ Genome index: ${GENOME_INDEX}"
echo "✓ Output directory: ${OUTDIR}"
echo "------------------------------------------------------------"

# ===== Go to Pipeline Directory =====
cd "${CODE_DIR}"

# ===== Process All Samples in Directory =====
echo "Searching for sample files in ${SAMPLE_PATH}..."
echo "------------------------------------------------------------"

# Count samples first
SAMPLE_COUNT=$(find "${SAMPLE_PATH}" -name "*_1.fq.gz" 2>/dev/null | wc -l)

if [[ ${SAMPLE_COUNT} -eq 0 ]]; then
    echo "ERROR: No samples found matching pattern *_1.fq.gz in ${SAMPLE_PATH}"
    exit 1
fi

echo "Found ${SAMPLE_COUNT} sample(s) to process"
echo "============================================================"

# Loop through all R1 files in the sample directory
for file in "${SAMPLE_PATH}"/*_1.fq.gz ; do
    # Extract file information
    withpath="${file}"
    filename="${withpath##*/}"
    base="${filename%_1.fq.gz}"
    
    # Define paired-end files
    READS_R1="${SAMPLE_PATH}/${base}_1.fq.gz"
    READS_R2="${SAMPLE_PATH}/${base}_2.fq.gz"
    
    echo ""
    echo "▶ Processing sample: ${base}"
    echo "  R1: ${READS_R1}"
    echo "  R2: ${READS_R2}"
    
    # Validate that both files exist
    if [[ ! -f "${READS_R1}" ]]; then
        echo "  ✗ ERROR: R1 file not found: ${READS_R1}"
        echo "  Skipping ${base}..."
        continue
    fi
    
    if [[ ! -f "${READS_R2}" ]]; then
        echo "  ✗ ERROR: R2 file not found: ${READS_R2}"
        echo "  Skipping ${base}..."
        continue
    fi
    
    echo "  ✓ Both read files found"
    
    # Build command for this sample
    CMD="python -m rna_pipeline.cli \
      --mode align \
      --genome-index ${GENOME_INDEX} \
      --reads-left ${READS_R1} \
      --reads-right ${READS_R2} \
      --outdir ${OUTDIR} \
      --threads ${THREADS} \
      ${QUANT_MODE} \
      --sample-name ${base}"
    
    echo "  Running alignment for ${base}..."
    echo "  ${CMD}"
    echo "  ------------------------------------------------------------"
    
    # Execute the command
    if eval ${CMD}; then
        echo "  ✓ ${base} completed successfully"
    else
        echo "  ✗ ${base} failed (exit code: $?)"
        echo "  Continuing with next sample..."
    fi
    
    echo "============================================================"
done

# ===== Summary =====
echo ""
echo "============================================================"
echo "Batch alignment completed!"
echo "End time: $(date)"
echo "Processed ${SAMPLE_COUNT} sample(s) from ${SAMPLE_DIR}"
echo "Output directory: ${OUTDIR}"
echo "============================================================"

# List output files
if [[ -d "${OUTDIR}" ]]; then
    echo "Generated files:"
    ls -lh "${OUTDIR}"
fi

