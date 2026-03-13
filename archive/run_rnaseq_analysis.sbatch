#!/bin/bash
# ===== SLURM directives =====
#SBATCH --job-name=rnaseq_analysis
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daisycortesj@vt.edu
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/05_rnaseq-code/scripts/
#SBATCH -o /projects/tholl_lab_1/daisy_analysis/06_analysis/rnaseq_analysis_output.txt
#SBATCH -e /projects/tholl_lab_1/daisy_analysis/06_analysis/rnaseq_analysis_error.txt

# ===== Job Script =====
# RNA-seq Count Matrix and Statistical Analysis Pipeline

set -e  # Exit on any error

# Activate conda environment
source ~/.bashrc
conda activate rnaseq

# Set paths (corrected for your directory structure)
STAR_COUNT_DIR="/projects/tholl_lab_1/daisy_analysis/03_count_tables/star"
OUTPUT_DIR="/projects/tholl_lab_1/daisy_analysis/06_analysis"
SCRIPT_DIR="/projects/tholl_lab_1/daisy_analysis/05_rnaseq-code/rna_pipeline/tools"

echo "Starting RNA-seq analysis pipeline..."
echo "Input directory: $STAR_COUNT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Script directory: $SCRIPT_DIR"

# Validate input directory
if [ ! -d "$STAR_COUNT_DIR" ]; then
    echo "ERROR: Directory does not exist: $STAR_COUNT_DIR"
    exit 1
fi

# Check for STAR count files
COUNT_FILES=$(find "$STAR_COUNT_DIR" -name "*ReadsPerGene.out.tab" | wc -l)
if [ "$COUNT_FILES" -eq 0 ]; then
    echo "ERROR: No ReadsPerGene.out.tab files found in $STAR_COUNT_DIR"
    exit 1
fi

echo "Found $COUNT_FILES STAR count files in $STAR_COUNT_DIR"

# Create output directory
mkdir -p "$OUTPUT_DIR"
COUNT_MATRIX_DIR="$OUTPUT_DIR/count_matrices"
mkdir -p "$COUNT_MATRIX_DIR"

echo "Output directory: $OUTPUT_DIR"

# Step 1: Build count matrix
echo "Step 1: Building gene count matrix..."

if [ ! -f "$SCRIPT_DIR/build_count_matrix.py" ]; then
    echo "ERROR: build_count_matrix.py not found in $SCRIPT_DIR"
    exit 1
fi

python3 "$SCRIPT_DIR/build_count_matrix.py" "$STAR_COUNT_DIR" -o "$COUNT_MATRIX_DIR"

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to build count matrix"
    exit 1
fi

echo "SUCCESS: Count matrix created successfully"

# Step 2: Check if R is available
echo "Step 2: Checking R installation..."

if ! command -v Rscript &> /dev/null; then
    echo "ERROR: Rscript not found. Please install R to run statistical analysis."
    echo "WARNING: Count matrix has been created. You can run the R analysis separately later."
    exit 1
fi

# Step 3: Run statistical analysis
echo "Step 3: Running differential expression analysis..."

if [ ! -f "$SCRIPT_DIR/rnaseq_analysis.R" ]; then
    echo "ERROR: rnaseq_analysis.R not found in $SCRIPT_DIR"
    exit 1
fi

Rscript "$SCRIPT_DIR/rnaseq_analysis.R" \
    "$COUNT_MATRIX_DIR/gene_count_matrix.tsv" \
    "$COUNT_MATRIX_DIR/sample_metadata.tsv" \
    "$OUTPUT_DIR/statistical_analysis"

if [ $? -ne 0 ]; then
    echo "ERROR: Statistical analysis failed"
    exit 1
fi

echo "SUCCESS: Statistical analysis completed successfully"
echo "Pipeline completed! Check $OUTPUT_DIR for all results."
