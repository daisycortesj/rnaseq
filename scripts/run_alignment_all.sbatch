#!/bin/bash
# ---------------------------------------------------------------------
# PURPOSE: Align RNA-seq reads for ALL sample directories using STAR
# USAGE:   sbatch run_alignment_all.sbatch <genome_type>
#          where genome_type is: nutmeg, carrot
# EXAMPLE: sbatch run_alignment_all.sbatch carrot
# NOTE:    This processes ALL sample directories (00_1_DC, 00_2_DG, 00_3_MF)
#          For a single directory, use run_alignment.sbatch
# ---------------------------------------------------------------------

# ===== SLURM directives =====
#SBATCH --job-name=star_align_all
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daisycortesj@vt.edu
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/02_mapped
#SBATCH -o bldstridx_output.txt
#SBATCH -e bldstridx_output.txt

set -euo pipefail

# ===== Configuration =====
GENOME_TYPE="${1:-}"
RAWDATA_DIR="/projects/tholl_lab_1/daisy_analysis/00_rawdata"
REFERENCE_DIR="/projects/tholl_lab_1/daisy_analysis/04_reference"
CODE_DIR="/projects/tholl_lab_1/daisy_analysis/05_rnaseq-code"
ALIGNMENT_DIR="/projects/tholl_lab_1/daisy_analysis/02_mapped"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

echo "============================================================"
echo "STAR Read Alignment Pipeline - Process All Samples"
echo "============================================================"
echo "SLURM job running on host: $(hostname)"
echo "CPUs: ${THREADS}  MEM: ${SLURM_MEM_PER_NODE:-64G}"
echo "Genome type: ${GENOME_TYPE}"
echo "Start time: $(date)"
echo "============================================================"

# ===== Validate Arguments =====
if [[ -z "${GENOME_TYPE}" ]]; then
    echo "ERROR: No genome type provided"
    echo "Usage: sbatch run_alignment_all.sbatch <genome_type>"
    echo "Genome types: nutmeg, carrot"
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

# ===== Set Genome-Specific Parameters =====
case "${GENOME_TYPE}" in
    nutmeg)
        echo "Aligning to NUTMEG genome (MYU assembly)"
        GENOME_INDEX="${REFERENCE_DIR}/nutmeg_star_index"
        QUANT_MODE="--no-quant-mode"
        ;;
    carrot)
        echo "Aligning to CARROT genome (with gene quantification)"
        GENOME_INDEX="${REFERENCE_DIR}/carrot_star_index"
        QUANT_MODE="--quant-mode"
        ;;
    *)
        echo "ERROR: Unknown genome type '${GENOME_TYPE}'"
        echo "Valid options: nutmeg, carrot"
        exit 1
        ;;
esac

# ===== Validate Genome Index =====
if [[ ! -d "${GENOME_INDEX}" ]]; then
    echo "ERROR: Genome index not found: ${GENOME_INDEX}"
    echo "Please build the index first using run_genome_index.sbatch"
    exit 1
fi
echo "✓ Genome index: ${GENOME_INDEX}"
echo "------------------------------------------------------------"

# ===== Go to Pipeline Directory =====
cd "${CODE_DIR}"

# ===== Process All Sample Directories =====
echo "Searching for sample directories in ${RAWDATA_DIR}..."
echo "============================================================"

TOTAL_SAMPLES=0
SUCCESSFUL_SAMPLES=0
FAILED_SAMPLES=0

# Loop through all sample directories (00_1_DC, 00_2_DG, 00_3_MF)
for sample_dir_path in "${RAWDATA_DIR}"/00_*/ ; do
    # Extract directory name
    SAMPLE_DIR=$(basename "${sample_dir_path}")
    
    echo ""
    echo "╔════════════════════════════════════════════════════════════"
    echo "║ Processing directory: ${SAMPLE_DIR}"
    echo "╚════════════════════════════════════════════════════════════"
    
    # Set output directory for this sample directory
    OUTDIR="${ALIGNMENT_DIR}/${SAMPLE_DIR}"
    echo "Output directory: ${OUTDIR}"
    
    # Create output directory
    mkdir -p "${OUTDIR}"
    
    # Count samples in this directory
    SAMPLE_COUNT=$(find "${sample_dir_path}" -name "*_1.fq.gz" 2>/dev/null | wc -l)
    
    if [[ ${SAMPLE_COUNT} -eq 0 ]]; then
        echo "⚠ No samples found in ${SAMPLE_DIR}, skipping..."
        continue
    fi
    
    echo "Found ${SAMPLE_COUNT} sample(s) in ${SAMPLE_DIR}"
    echo "------------------------------------------------------------"
    
    # Loop through all samples in this directory
    for file in "${sample_dir_path}"*_1.fq.gz ; do
        # Extract file information
        withpath="${file}"
        filename="${withpath##*/}"
        base="${filename%_1.fq.gz}"
        
        # Define paired-end files
        READS_R1="${sample_dir_path}${base}_1.fq.gz"
        READS_R2="${sample_dir_path}${base}_2.fq.gz"
        
        TOTAL_SAMPLES=$((TOTAL_SAMPLES + 1))
        
        echo ""
        echo "▶ Sample ${TOTAL_SAMPLES}: ${base}"
        echo "  R1: ${READS_R1}"
        echo "  R2: ${READS_R2}"
        
        # Validate that both files exist
        if [[ ! -f "${READS_R1}" ]]; then
            echo "  ✗ ERROR: R1 file not found"
            FAILED_SAMPLES=$((FAILED_SAMPLES + 1))
            continue
        fi
        
        if [[ ! -f "${READS_R2}" ]]; then
            echo "  ✗ ERROR: R2 file not found"
            FAILED_SAMPLES=$((FAILED_SAMPLES + 1))
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
        
        echo "  Running alignment..."
        
        # Execute the command
        if eval ${CMD}; then
            echo "  ✓ ${base} completed successfully"
            SUCCESSFUL_SAMPLES=$((SUCCESSFUL_SAMPLES + 1))
        else
            echo "  ✗ ${base} failed (exit code: $?)"
            FAILED_SAMPLES=$((FAILED_SAMPLES + 1))
            echo "  Continuing with next sample..."
        fi
        
        echo "  ------------------------------------------------------------"
    done
    
    echo ""
    echo "Completed directory: ${SAMPLE_DIR}"
    echo "============================================================"
done

# ===== Final Summary =====
echo ""
echo "╔════════════════════════════════════════════════════════════"
echo "║ FINAL SUMMARY"
echo "╚════════════════════════════════════════════════════════════"
echo "Total samples processed: ${TOTAL_SAMPLES}"
echo "Successful: ${SUCCESSFUL_SAMPLES}"
echo "Failed: ${FAILED_SAMPLES}"
echo "Genome: ${GENOME_TYPE}"
echo "End time: $(date)"
echo "============================================================"

# List all output directories
echo ""
echo "Output directories created:"
ls -ld "${ALIGNMENT_DIR}"/00_*/ 2>/dev/null || echo "(none)"

exit 0

