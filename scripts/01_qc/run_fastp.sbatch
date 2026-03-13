#!/bin/bash
# ---------------------------------------------------------------------
# PURPOSE: Clean/filter all raw FASTQ files using fastp
# USAGE:   sbatch fastp_all.sbatch
# OUTPUT:  Cleaned reads in /projects/tholl_lab_1/daisy_analysis/01_processed/fastp/
#          Quality reports: *_fastp.html
# ---------------------------------------------------------------------

# ===== SLURM directives =====
#SBATCH --job-name=fastp_clean
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daisycortesj@vt.edu
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis
#SBATCH -o 01_processed/fastp_%j.out
#SBATCH -e 01_processed/fastp_%j.err

set -euo pipefail

# ===== Configuration =====
RAWDATA_DIR="/projects/tholl_lab_1/daisy_analysis/00_rawdata"
FASTP_OUTPUT="/projects/tholl_lab_1/daisy_analysis/01_processed/fastp"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
MIN_QUALITY=20    # Minimum quality score (Q20 = 99% accurate)
MIN_LENGTH=50     # Minimum read length after trimming

echo "============================================================"
echo "fastp - Quality Trimming & Filtering (All Samples)"
echo "============================================================"
echo "SLURM job running on host: $(hostname)"
echo "CPUs: ${THREADS}  MEM: ${SLURM_MEM_PER_NODE:-16G}"
echo "Start time: $(date)"
echo "Parameters:"
echo "  Minimum quality: Q${MIN_QUALITY}"
echo "  Minimum length:  ${MIN_LENGTH}bp"
echo "============================================================"

# ===== Create output directory =====
mkdir -p "${FASTP_OUTPUT}"

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

# Check if fastp is installed
if ! command -v fastp &> /dev/null; then
    echo "fastp not found. Installing..."
    conda install -c bioconda fastp -y
fi

# Check tool version
echo "Checking tool version:"
fastp --version 2>&1 | head -1
echo "------------------------------------------------------------"

# ===== Find all R1 files (forward reads) =====
echo "Searching for FASTQ files in ${RAWDATA_DIR}..."

# Find all R1/forward files
SAMPLE_COUNT=0
CLEANED_COUNT=0
FAILED_COUNT=0

for R1 in $(find "${RAWDATA_DIR}" -name "*_1.fq.gz" -o -name "*_R1.fq.gz" | sort); do
    # Get corresponding R2 file
    R2="${R1/_1.fq.gz/_2.fq.gz}"
    R2="${R2/_R1.fq.gz/_R2.fq.gz}"
    
    # Extract sample name
    SAMPLE=$(basename "${R1}" | sed 's/_1\.fq\.gz$//' | sed 's/_R1\.fq\.gz$//')
    
    # Check if R2 exists
    if [[ ! -f "${R2}" ]]; then
        echo "ERROR: R2 file not found for ${SAMPLE}"
        echo "  Expected: ${R2}"
        ((FAILED_COUNT++))
        continue
    fi
    
    # Output filenames
    OUT_R1="${FASTP_OUTPUT}/${SAMPLE}_1_clean.fq.gz"
    OUT_R2="${FASTP_OUTPUT}/${SAMPLE}_2_clean.fq.gz"
    JSON_REPORT="${FASTP_OUTPUT}/${SAMPLE}_fastp.json"
    HTML_REPORT="${FASTP_OUTPUT}/${SAMPLE}_fastp.html"
    
    # Check if already processed
    if [[ -f "${OUT_R1}" ]] && [[ -f "${OUT_R2}" ]] && [[ -f "${HTML_REPORT}" ]]; then
        echo "⏭️  Skipping ${SAMPLE} (already processed)"
        ((SAMPLE_COUNT++))
        continue
    fi
    
    echo "------------------------------------------------------------"
    echo "Processing sample: ${SAMPLE}"
    echo "  Input R1:  ${R1}"
    echo "  Input R2:  ${R2}"
    echo "  Output R1: ${OUT_R1}"
    echo "  Output R2: ${OUT_R2}"
    
    # Run fastp
    fastp \
        -i "${R1}" \
        -I "${R2}" \
        -o "${OUT_R1}" \
        -O "${OUT_R2}" \
        --json "${JSON_REPORT}" \
        --html "${HTML_REPORT}" \
        --detect_adapter_for_pe \
        --qualified_quality_phred ${MIN_QUALITY} \
        --length_required ${MIN_LENGTH} \
        --thread ${THREADS} \
        2>&1 | grep -E "Read1 before filtering:|Read1 after filtering:|Read2 before filtering:|Read2 after filtering:"
    
    # Check if successful
    if [[ $? -eq 0 ]] && [[ -f "${OUT_R1}" ]] && [[ -f "${OUT_R2}" ]]; then
        echo "✓ ${SAMPLE} complete"
        ((CLEANED_COUNT++))
    else
        echo "✗ ${SAMPLE} failed"
        ((FAILED_COUNT++))
    fi
    
    ((SAMPLE_COUNT++))
done

echo "============================================================"
echo "fastp processing completed!"
echo "End time: $(date)"
echo ""
echo "SUMMARY:"
echo "  Total samples processed: ${SAMPLE_COUNT}"
echo "  Successfully cleaned:    ${CLEANED_COUNT}"
echo "  Failed:                  ${FAILED_COUNT}"
echo ""
echo "RESULTS:"
echo "  Cleaned reads: ${FASTP_OUTPUT}/*_clean.fq.gz"
echo "  Reports:       ${FASTP_OUTPUT}/*_fastp.html"
echo ""
echo "NEXT STEPS:"
echo "1. Review fastp reports:"
echo "   Download: ${FASTP_OUTPUT}/*_fastp.html"
echo "2. Re-run QC on cleaned data:"
echo "   python -m rna_pipeline.cli --mode qc \\"
echo "     --fastq-dir ${FASTP_OUTPUT} \\"
echo "     --outdir 01_processed/qc_fastp"
echo "3. Compare QC before/after:"
echo "   Before: 01_processed/qc/multiqc_report.html"
echo "   After:  01_processed/qc_fastp/multiqc_report.html"
echo "4. Use cleaned reads for Trinity/STAR:"
echo "   sbatch scripts/run_trinity.sbatch (with cleaned files)"
echo "============================================================"

# ===== Generate summary statistics =====
echo ""
echo "Quick statistics (from fastp reports):"
if command -v jq &> /dev/null; then
    # If jq is available, parse JSON reports nicely
    for JSON in "${FASTP_OUTPUT}"/*_fastp.json; do
        if [[ -f "${JSON}" ]]; then
            SAMPLE=$(basename "${JSON}" _fastp.json)
            BEFORE=$(jq -r '.summary.before_filtering.total_reads' "${JSON}" 2>/dev/null || echo "N/A")
            AFTER=$(jq -r '.summary.after_filtering.total_reads' "${JSON}" 2>/dev/null || echo "N/A")
            if [[ "${BEFORE}" != "N/A" ]] && [[ "${AFTER}" != "N/A" ]]; then
                PERCENT=$(awk "BEGIN {printf \"%.1f\", (${AFTER}/${BEFORE})*100}")
                echo "  ${SAMPLE}: ${AFTER}/${BEFORE} reads retained (${PERCENT}%)"
            fi
        fi
    done
else
    echo "  (Install 'jq' for detailed statistics: conda install -c conda-forge jq)"
fi
