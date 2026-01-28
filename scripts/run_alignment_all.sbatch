#!/bin/bash
# ---------------------------------------------------------------------
# PURPOSE: Align RNA-seq reads for ALL sample directories using STAR
# USAGE:   sbatch run_alignment_all.sbatch
# NOTE:    This processes ALL sample directories with auto-mapped genomes:
#          - 00_1_DC → carrot genome
#          - 00_2_DG → carrot genome
#          - 00_3_MF → nutmeg genome
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
#SBATCH -o alignment_all_%j.out  
#SBATCH -e alignment_all_%j.err       

set -euo pipefail

# ===== Configuration =====
# ┌─────────────────────────────────────────────────────────────┐
# │ IMPORTANT: Choose your input data source                   │
# │ Set INPUT_DATA to:                                         │
# │   "raw"   → Use original FASTQ files (if QC is good)      │
# │   "fastp" → Use cleaned FASTQ files (if QC shows problems)│
# └─────────────────────────────────────────────────────────────┘

INPUT_DATA="raw"           # Change to "fastp" after cleaning

# Directory paths
RAWDATA_DIR="/projects/tholl_lab_1/daisy_analysis/00_rawdata"
FASTP_DIR="/projects/tholl_lab_1/daisy_analysis/01_processed/fastp"
REFERENCE_DIR="/projects/tholl_lab_1/daisy_analysis/04_reference"
CODE_DIR="/projects/tholl_lab_1/daisy_analysis/05_rnaseq-code"
ALIGNMENT_DIR="/projects/tholl_lab_1/daisy_analysis/02_mapped"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

# Determine input directory based on choice
if [[ "${INPUT_DATA}" == "fastp" ]]; then
    INPUT_DIR="${FASTP_DIR}"
    INPUT_SOURCE="FASTP-CLEANED"
    FILE_SUFFIX="_clean"
else
    INPUT_DIR="${RAWDATA_DIR}"
    INPUT_SOURCE="RAW"
    FILE_SUFFIX=""
fi

echo "============================================================"
echo "STAR Read Alignment Pipeline - Process All Samples"
echo "============================================================"
echo "SLURM job running on host: $(hostname)"
echo "CPUs: ${THREADS}  MEM: ${SLURM_MEM_PER_NODE:-64G}"
echo ""
echo "┌─────────────────────────────────────────────────────────┐"
echo "│ DATA SOURCE: ${INPUT_SOURCE}"
echo "└─────────────────────────────────────────────────────────┘"
echo "  Input: ${INPUT_DIR}"
echo ""
echo "Genome mapping:"
echo "  - 00_1_DC → carrot"
echo "  - 00_2_DG → carrot"
echo "  - 00_3_MF → nutmeg"
echo "Start time: $(date)"
echo "============================================================"

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

# ===== Go to Pipeline Directory =====
cd "${CODE_DIR}"

# ===== Process All Sample Directories =====
echo "Searching for sample directories in ${INPUT_DIR}..."
echo "============================================================"

TOTAL_SAMPLES=0
SUCCESSFUL_SAMPLES=0
FAILED_SAMPLES=0

# Loop through all sample directories (00_1_DC, 00_2_DG, 00_3_MF)
for sample_dir_path in "${INPUT_DIR}"/00_*/ ; do
    # Extract directory name
    SAMPLE_DIR=$(basename "${sample_dir_path}")
    
    echo ""
    echo "╔════════════════════════════════════════════════════════════"
    echo "║ Processing directory: ${SAMPLE_DIR}"
    echo "╚════════════════════════════════════════════════════════════"
    
    # Determine genome based on directory name
    case "${SAMPLE_DIR}" in
        00_1_DC|00_2_DG)
            GENOME_TYPE="carrot"
            GENOME_INDEX="${REFERENCE_DIR}/carrot_star_index"
            QUANT_MODE="--quant-mode"
            echo "Using CARROT genome (with gene quantification)"
            ;;
        00_3_MF)
            GENOME_TYPE="nutmeg"
            GENOME_INDEX="${REFERENCE_DIR}/nutmeg_star_index"
            QUANT_MODE="--no-quant-mode"
            echo "Using NUTMEG genome (MYU assembly)"
            ;;
        *)
            echo "⚠ WARNING: Unknown directory '${SAMPLE_DIR}', skipping..."
            continue
            ;;
    esac
    
    # Validate genome index exists
    if [[ ! -d "${GENOME_INDEX}" ]]; then
        echo "✗ ERROR: Genome index not found: ${GENOME_INDEX}"
        echo "Please build the index first using run_genome_index.sbatch"
        continue
    fi
    echo "✓ Genome index: ${GENOME_INDEX}"
    
    # Set output directory for this sample directory
    OUTDIR="${ALIGNMENT_DIR}/${SAMPLE_DIR}"
    echo "Output directory: ${OUTDIR}"
    
    # Create output directory
    mkdir -p "${OUTDIR}"
    
    # Count samples in this directory (adjust pattern based on data source)
    if [[ "${INPUT_DATA}" == "fastp" ]]; then
        SAMPLE_COUNT=$(find "${sample_dir_path}" -name "*_1_clean.fq.gz" 2>/dev/null | wc -l)
    else
        SAMPLE_COUNT=$(find "${sample_dir_path}" -name "*_1.fq.gz" 2>/dev/null | wc -l)
    fi
    
    if [[ ${SAMPLE_COUNT} -eq 0 ]]; then
        echo "⚠ No samples found in ${SAMPLE_DIR}, skipping..."
        continue
    fi
    
    echo "Found ${SAMPLE_COUNT} sample(s) in ${SAMPLE_DIR}"
    echo "------------------------------------------------------------"
    
    # Determine file pattern based on data source
    if [[ "${INPUT_DATA}" == "fastp" ]]; then
        FILE_PATTERN="*_1_clean.fq.gz"
    else
        FILE_PATTERN="*_1.fq.gz"
    fi
    
    # Loop through all samples in this directory
    for file in "${sample_dir_path}"${FILE_PATTERN} ; do
        # Extract file information
        withpath="${file}"
        filename="${withpath##*/}"
        base="${filename%_1${FILE_SUFFIX}.fq.gz}"
        
        # Define paired-end files
        READS_R1="${sample_dir_path}${base}_1${FILE_SUFFIX}.fq.gz"
        READS_R2="${sample_dir_path}${base}_2${FILE_SUFFIX}.fq.gz"
        
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
echo ""
echo "DATA USED: ${INPUT_SOURCE}"
if [[ "${INPUT_SOURCE}" == "RAW" ]]; then
    echo "  (If alignment rates are low, consider using fastp-cleaned data)"
fi
echo ""
echo "Genome mapping used:"
echo "  - 00_1_DC → carrot"
echo "  - 00_2_DG → carrot"
echo "  - 00_3_MF → nutmeg"
echo "End time: $(date)"
echo "============================================================"

# List all output directories
echo ""
echo "Output directories created:"
ls -ld "${ALIGNMENT_DIR}"/00_*/ 2>/dev/null || echo "(none)"

exit 0
