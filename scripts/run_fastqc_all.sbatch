#!/bin/bash
# ---------------------------------------------------------------------
# PURPOSE: Run FastQC on all raw FASTQ files and generate MultiQC report
# USAGE:   sbatch run_fastqc_all.sbatch
# OUTPUT:  QC reports in /projects/tholl_lab_1/daisy_analysis/01_processed/qc/
#          Log files in /projects/tholl_lab_1/daisy_analysis/01_processed/
# ---------------------------------------------------------------------

# ===== SLURM directives =====
#SBATCH --job-name=fastqc
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daisycortesj@vt.edu
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis
#SBATCH -o 01_processed/qc_%j.out
#SBATCH -e 01_processed/qc_%j.err

set -euo pipefail

# ===== Configuration =====
RAWDATA_DIR="/projects/tholl_lab_1/daisy_analysis/00_rawdata"
QC_OUTPUT="/projects/tholl_lab_1/daisy_analysis/01_processed/qc"
CODE_DIR="/projects/tholl_lab_1/daisy_analysis/05_rnaseq-code"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

echo "============================================================"
echo "FastQC Quality Control - All Samples"
echo "============================================================"
echo "SLURM job running on host: $(hostname)"
echo "CPUs: ${THREADS}  MEM: ${SLURM_MEM_PER_NODE:-16G}"
echo "Start time: $(date)"
echo "============================================================"

# ===== Create output directory =====
mkdir -p "${QC_OUTPUT}"

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

echo "------------------------------------------------------------"

# ===== Go to Pipeline Directory =====
cd "${CODE_DIR}"

# ===== Run QC Pipeline =====
# The pipeline will automatically:
#   1. Check if FastQC/MultiQC are installed
#   2. Install them if missing (using conda)
#   3. Find all FASTQ files in the directory
#   4. Run FastQC on all files
#   5. Aggregate results with MultiQC

CMD="python -m rna_pipeline.cli \
  --mode qc \
  --fastq-dir \"${RAWDATA_DIR}\" \
  --outdir \"${QC_OUTPUT}\" \
  --threads ${THREADS}"

echo "Running QC pipeline:"
echo "${CMD}"
echo "------------------------------------------------------------"

eval ${CMD}

echo "------------------------------------------------------------"

# ===== Summary =====
echo "============================================================"
echo "Quality control completed successfully!"
echo "End time: $(date)"
echo ""
echo "RESULTS:"
echo "  Individual FastQC reports: ${QC_OUTPUT}/fastqc_results/"
echo "  MultiQC summary report:    ${QC_OUTPUT}/multiqc_report.html"
echo "  Job logs:                  01_processed/qc_${SLURM_JOB_ID}.out"
echo ""
echo "NEXT STEPS:"
echo "1. Download and open: ${QC_OUTPUT}/multiqc_report.html"
echo "2. Review using the guide: ReadME/QC_INTERPRETATION_GUIDE.md"
echo "3. Look for:"
echo "   - Per base quality scores (should be >28, ideally >30)"
echo "   - Adapter contamination (should be <5%)"
echo "   - Per sequence GC content (should match your organism)"
echo "4. If quality is good → proceed with confidence!"
echo "5. If quality issues found → consider trimming and re-running"
echo "============================================================"
