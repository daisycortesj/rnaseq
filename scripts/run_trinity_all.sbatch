#!/bin/bash
# ---------------------------------------------------------------------
# PURPOSE: Run Trinity de novo assembly for ALL sample directories
# USAGE:   sbatch run_trinity_all.sbatch
# NOTE:    This processes ALL samples in ALL directories (00_1_DC, 00_2_DG, 00_3_MF)
#          For a single directory, use run_trinity_batch.sbatch
#          For a single sample, use run_trinity.sbatch
# ---------------------------------------------------------------------

# ===== SLURM directives =====
#SBATCH --job-name=trinity_all
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G
#SBATCH --time=120:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daisycortesj@vt.edu
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/01_processed
#SBATCH -o trinity_all_%j.out
#SBATCH -e trinity_all_%j.err

set -euo pipefail

# ===== Configuration =====
RAWDATA_DIR="/projects/tholl_lab_1/daisy_analysis/00_rawdata"
CODE_DIR="/projects/tholl_lab_1/daisy_analysis/05_rnaseq-code"
TRINITY_OUTDIR="/projects/tholl_lab_1/daisy_analysis/01_processed"
THREADS="${SLURM_CPUS_PER_TASK:-16}"
MEM_GB="128"

echo "============================================================"
echo "Trinity De Novo Assembly - Process All Samples"
echo "============================================================"
echo "SLURM job running on host: $(hostname)"
echo "CPUs: ${THREADS}  MEM: ${SLURM_MEM_PER_NODE:-128G}"
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
Trinity --version || echo "Trinity not found."
echo "------------------------------------------------------------"

# ===== Go to Pipeline Directory =====
cd "${CODE_DIR}"

# ===== Process All Sample Directories =====
echo "Searching for sample directories in ${RAWDATA_DIR}..."
echo "============================================================"

TOTAL_SAMPLES=0
SUCCESSFUL_SAMPLES=0
FAILED_SAMPLES=0
SKIPPED_SAMPLES=0

# Loop through all sample directories (00_1_DC, 00_2_DG, 00_3_MF)
for sample_dir_path in "${RAWDATA_DIR}"/00_*/ ; do
    # Extract directory name
    SAMPLE_DIR=$(basename "${sample_dir_path}")
    
    echo ""
    echo "╔════════════════════════════════════════════════════════════"
    echo "║ Processing directory: ${SAMPLE_DIR}"
    echo "╚════════════════════════════════════════════════════════════"
    
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
        LEFT_FILE="${sample_dir_path}${base}_1.fq.gz"
        RIGHT_FILE="${sample_dir_path}${base}_2.fq.gz"
        
        TOTAL_SAMPLES=$((TOTAL_SAMPLES + 1))
        
        echo ""
        echo "▶ Sample ${TOTAL_SAMPLES}: ${base}"
        echo "  Left:  ${LEFT_FILE}"
        echo "  Right: ${RIGHT_FILE}"
        
        # Validate that both files exist and are readable
        if [[ ! -r "${LEFT_FILE}" ]]; then
            echo "  ✗ ERROR: Cannot read left file"
            FAILED_SAMPLES=$((FAILED_SAMPLES + 1))
            continue
        fi
        
        if [[ ! -r "${RIGHT_FILE}" ]]; then
            echo "  ✗ ERROR: Cannot read right file"
            FAILED_SAMPLES=$((FAILED_SAMPLES + 1))
            continue
        fi
        
        echo "  ✓ Both read files found and readable"
        
        # Check if symlinks
        if [[ -L "${LEFT_FILE}" ]]; then
            echo "  (symlink to: $(readlink -f "${LEFT_FILE}"))"
        fi
        
        # Set output directory for this sample
        OUTPUT_DIR="${TRINITY_OUTDIR}/${base}_trinity_assembly"
        echo "  Output: ${OUTPUT_DIR}"
        
        # Check if output directory already exists
        if [[ -d "${OUTPUT_DIR}" ]]; then
            echo "  ⚠ WARNING: Output directory already exists"
            echo "  Skipping ${base}... (remove directory first to reprocess)"
            SKIPPED_SAMPLES=$((SKIPPED_SAMPLES + 1))
            continue
        fi
        
        # Build command for this sample
        CMD="python -m rna_pipeline.cli \
          --mode trinity \
          --reads-left \"${LEFT_FILE}\" \
          --reads-right \"${RIGHT_FILE}\" \
          --outdir \"${OUTPUT_DIR}\" \
          --threads ${THREADS} \
          --mem-gb ${MEM_GB}"
        
        echo "  Running Trinity assembly..."
        
        # Execute the command
        echo "  Command: ${CMD}"
        if eval ${CMD} 2>&1; then
            echo "  ✓ ${base} completed successfully"
            SUCCESSFUL_SAMPLES=$((SUCCESSFUL_SAMPLES + 1))
            
            # Check for Trinity output
            if [[ -f "${OUTPUT_DIR}/Trinity.fasta" ]]; then
                CONTIG_COUNT=$(grep -c ">" "${OUTPUT_DIR}/Trinity.fasta" || echo "0")
                echo "  ✓ Trinity.fasta: ${CONTIG_COUNT} contigs"
            else
                echo "  ⚠ Warning: Trinity.fasta not found"
            fi
        else
            EXIT_CODE=$?
            echo "  ✗ ${base} failed (exit code: ${EXIT_CODE})"
            echo "  Error details logged to: trinity_all_${SLURM_JOB_ID}.err"
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
echo "Total samples found: ${TOTAL_SAMPLES}"
echo "Successful: ${SUCCESSFUL_SAMPLES}"
echo "Failed: ${FAILED_SAMPLES}"
echo "Skipped (already exists): ${SKIPPED_SAMPLES}"
echo "End time: $(date)"
echo "============================================================"

# List all output assemblies
echo ""
echo "Generated assemblies:"
ls -1d "${TRINITY_OUTDIR}"/*_trinity_assembly 2>/dev/null || echo "(none)"

exit 0

