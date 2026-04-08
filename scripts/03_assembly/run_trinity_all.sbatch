#!/bin/bash
# ---------------------------------------------------------------------
# PURPOSE: Run Trinity de novo assembly (one species or all)
# USAGE:   sbatch run_trinity_all.sbatch [SPECIES_CODE]
#
# EXAMPLES:
#   sbatch run_trinity_all.sbatch MF    # only nutmeg
#   sbatch run_trinity_all.sbatch DC    # only D. carota
#   sbatch run_trinity_all.sbatch       # all species
#
# AUTO-DETECTION (no manual settings needed):
#   For each species, the script checks if fastp-cleaned reads exist
#   in 01_processed/. If they do, it uses those. If not, it falls
#   back to raw reads from 00_rawdata/.
#
#   Example:
#     00_1_DC → no *_clean.fq.gz in 01_processed/ → uses RAW reads
#     00_3_MF → has *_clean.fq.gz in 01_processed/ → uses CLEAN reads
#
# PATHS: All directory paths come from scripts/config.sh
# ---------------------------------------------------------------------

# ===== SLURM directives =====
#SBATCH --job-name=trinity_all
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G
#SBATCH --time=164:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daisycortesj@vt.edu
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/01_processed
#SBATCH -o trinity_all_%j.out
#SBATCH -e trinity_all_%j.err

set -euo pipefail

# ===== Load shared config =====
source "${SLURM_SUBMIT_DIR:-.}/scripts/config.sh"

# --- Optional: filter to a single species ---
# If a species code is given (e.g., MF), only process that species.
# If no code is given, process all species directories.
SPECIES_FILTER="${1:-}"
if [[ -n "${SPECIES_FILTER}" ]]; then
    get_sample_info "${SPECIES_FILTER}"
    FILTER_DIR="${SPECIES_DIR}"
fi

# ===== Configuration =====
TRINITY_OUTDIR="${PROCESSED_DIR}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"
MEM_GB="128"

echo "============================================================"
if [[ -n "${SPECIES_FILTER}" ]]; then
    echo "Trinity De Novo Assembly — ${SPECIES_NAME} (${SPECIES_FILTER})"
else
    echo "Trinity De Novo Assembly - Process All Samples"
fi
echo "============================================================"
echo "SLURM job running on host: $(hostname)"
echo "CPUs: ${THREADS}  MEM: ${SLURM_MEM_PER_NODE:-128G}"
echo "Start time: $(date)"
echo "============================================================"

# ===== Activate Conda Environment =====
if ! bash -c "source ~/.bashrc" 2>/dev/null; then
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate rnaseq
else
    source ~/.bashrc
    conda activate rnaseq
fi

# Make sure Trinity is available
if ! command -v Trinity &> /dev/null; then
    echo "ERROR: Trinity not found. Install with: conda install -c bioconda trinity"
    exit 1
fi

echo "Trinity version: $(Trinity --version 2>&1 | head -1)"
echo "------------------------------------------------------------"

# ===== Go to Pipeline Directory =====
cd "${CODE_DIR}"

# ===== Process All Sample Directories =====
echo "Scanning species directories in ${RAWDATA_DIR}/ ..."
echo "============================================================"

TOTAL_SAMPLES=0
SUCCESSFUL_SAMPLES=0
FAILED_SAMPLES=0
SKIPPED_SAMPLES=0

# ====================================================================
# MAIN PROCESSING LOOP
# ====================================================================
# We discover species from 00_rawdata/ (the canonical list).
# For each species, we auto-detect whether fastp-cleaned reads exist
# in 01_processed/. If they do, we use those. Otherwise we fall back
# to raw reads. This means:
#   DC/DG (no fastp) → reads from 00_rawdata/00_1_DC/
#   MF    (has fastp) → reads from 01_processed/00_3_MF/
# ====================================================================
for sample_dir_path in "${RAWDATA_DIR}"/00_*/ ; do
    SAMPLE_DIR=$(basename "${sample_dir_path}")
    
    # If a species filter was given, skip directories that don't match
    if [[ -n "${SPECIES_FILTER}" ]] && [[ "${SAMPLE_DIR}" != "${FILTER_DIR}" ]]; then
        continue
    fi
    
    echo ""
    echo "╔════════════════════════════════════════════════════════════"
    echo "║ Processing directory: ${SAMPLE_DIR}"
    echo "╚════════════════════════════════════════════════════════════"
    
    # --- Auto-detect: use clean reads if available, otherwise raw ---
    # Check if fastp-cleaned reads exist in 01_processed/<SPECIES_DIR>/
    CLEAN_DIR="${PROCESSED_DIR}/${SAMPLE_DIR}"
    CLEAN_COUNT=$(find "${CLEAN_DIR}" -maxdepth 1 -name "*_1_clean.fq.gz" 2>/dev/null | wc -l)
    
    if [[ "${CLEAN_COUNT}" -gt 0 ]]; then
        # Fastp-cleaned reads found → use them
        SPECIES_INPUT_DIR="${CLEAN_DIR}"
        FILE_PATTERN="*_1_clean.fq.gz"
        FILE_SUFFIX="_clean"
        INPUT_SOURCE="FASTP-CLEANED"
    else
        # No clean reads → fall back to raw reads
        SPECIES_INPUT_DIR="${RAWDATA_DIR}/${SAMPLE_DIR}"
        FILE_PATTERN="*_1.fq.gz"
        FILE_SUFFIX=""
        INPUT_SOURCE="RAW"
    fi
    
    echo "  Data source: ${INPUT_SOURCE}"
    echo "  Input dir:   ${SPECIES_INPUT_DIR}"
    
    # Count samples in this directory
    SAMPLE_COUNT=$(find "${SPECIES_INPUT_DIR}" -maxdepth 1 -name "${FILE_PATTERN}" 2>/dev/null | wc -l)
    
    if [[ ${SAMPLE_COUNT} -eq 0 ]]; then
        echo "  ⚠ No samples found in ${SAMPLE_DIR}, skipping..."
        continue
    fi
    
    echo "  Samples found: ${SAMPLE_COUNT}"
    echo "------------------------------------------------------------"
    
    # Loop through all R1 files for this species
    for file in "${SPECIES_INPUT_DIR}"/${FILE_PATTERN} ; do
        # Extract sample base name
        # Raw:   MFF1_1.fq.gz       → base = MFF1
        # Clean: MFF1_1_clean.fq.gz → base = MFF1
        filename="${file##*/}"
        base="${filename%_1${FILE_SUFFIX}.fq.gz}"
        
        # Define paired-end files
        LEFT_FILE="${SPECIES_INPUT_DIR}/${base}_1${FILE_SUFFIX}.fq.gz"
        RIGHT_FILE="${SPECIES_INPUT_DIR}/${base}_2${FILE_SUFFIX}.fq.gz"
        
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
        RESUME_FLAG=""
        if [[ -d "${OUTPUT_DIR}" ]]; then
            # Check if Trinity.fasta exists (completed assembly)
            if [[ -f "${OUTPUT_DIR}/Trinity.fasta" ]]; then
                echo "  ✓ Assembly already complete (Trinity.fasta exists)"
                echo "  Skipping ${base}..."
                SKIPPED_SAMPLES=$((SKIPPED_SAMPLES + 1))
                continue
            else
                echo "  ⚠ Incomplete assembly detected - will RESUME"
                echo "  (Trinity.fasta not found, checkpoint files may exist)"
                RESUME_FLAG="--resume"
            fi
        fi
        
        # Build command for this sample
        CMD="python -m rna_pipeline.cli \
          --mode trinity \
          --reads-left \"${LEFT_FILE}\" \
          --reads-right \"${RIGHT_FILE}\" \
          --outdir \"${OUTPUT_DIR}\" \
          --threads ${THREADS} \
          --mem-gb ${MEM_GB} \
          ${RESUME_FLAG}"
        
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
echo ""
echo "Data source: auto-detected per species (clean if available, else raw)"
echo "End time: $(date)"
echo "============================================================"

# List all output assemblies
echo ""
echo "Generated assemblies:"
ls -1d "${TRINITY_OUTDIR}"/*_trinity_assembly 2>/dev/null || echo "(none)"

exit 0

