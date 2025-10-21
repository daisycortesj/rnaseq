#!/bin/bash
# RNA-seq Count Matrix and Statistical Analysis Pipeline
# 
# This script builds a gene count matrix from STAR output files
# and performs differential expression analysis with DESeq2 and edgeR.
#
# Usage: ./run_rnaseq_analysis.sh [star_count_dir] [output_dir]

set -e  # Exit on any error

# Default values
STAR_COUNT_DIR=""
OUTPUT_DIR="rnaseq_analysis_results"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Function to show usage
show_usage() {
    echo "Usage: $0 <star_count_dir> [output_dir]"
    echo ""
    echo "Arguments:"
    echo "  star_count_dir  Directory containing STAR ReadsPerGene.out.tab files"
    echo "  output_dir      Output directory (default: rnaseq_analysis_results)"
    echo ""
    echo "Example:"
    echo "  $0 /path/to/star/counts"
    echo "  $0 /path/to/star/counts my_analysis_results"
    echo ""
    echo "The script will:"
    echo "  1. Build gene count matrix from STAR output"
    echo "  2. Create sample metadata"
    echo "  3. Run DESeq2 and edgeR analysis"
    echo "  4. Generate quality control plots"
    echo "  5. Create summary reports"
}

# Parse command line arguments
if [ $# -eq 0 ]; then
    show_usage
    exit 1
fi

STAR_COUNT_DIR="$1"
if [ $# -ge 2 ]; then
    OUTPUT_DIR="$2"
fi

# Validate input directory
if [ ! -d "$STAR_COUNT_DIR" ]; then
    print_error "Directory does not exist: $STAR_COUNT_DIR"
    exit 1
fi

# Check for STAR count files
COUNT_FILES=$(find "$STAR_COUNT_DIR" -name "*ReadsPerGene.out.tab" | wc -l)
if [ "$COUNT_FILES" -eq 0 ]; then
    print_error "No ReadsPerGene.out.tab files found in $STAR_COUNT_DIR"
    exit 1
fi

print_status "Found $COUNT_FILES STAR count files in $STAR_COUNT_DIR"

# Create output directory
mkdir -p "$OUTPUT_DIR"
COUNT_MATRIX_DIR="$OUTPUT_DIR/count_matrices"
mkdir -p "$COUNT_MATRIX_DIR"

print_status "Output directory: $OUTPUT_DIR"

# Step 1: Build count matrix
print_status "Step 1: Building gene count matrix..."

if [ ! -f "$SCRIPT_DIR/build_count_matrix.py" ]; then
    print_error "build_count_matrix.py not found in $SCRIPT_DIR"
    exit 1
fi

python3 "$SCRIPT_DIR/build_count_matrix.py" "$STAR_COUNT_DIR" -o "$COUNT_MATRIX_DIR"

if [ $? -ne 0 ]; then
    print_error "Failed to build count matrix"
    exit 1
fi

print_success "Count matrix created successfully"

# Step 2: Check if R is available
print_status "Step 2: Checking R installation..."

if ! command -v Rscript &> /dev/null; then
    print_error "Rscript not found. Please install R to run statistical analysis."
    print_warning "Count matrix has been created. You can run the R analysis separately later."
    exit 1
fi

# Check for required R packages
print_status "Checking required R packages..."

Rscript -e "
required_packages <- c('DESeq2', 'edgeR', 'ggplot2', 'pheatmap', 'RColorBrewer', 'dplyr', 'EnhancedVolcano')
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
    cat('Missing R packages:', paste(missing_packages, collapse = ', '), '\n')
    cat('Install with: install.packages(c(', paste(paste0('\"', missing_packages, '\"'), collapse = ', '), '))\n')
    quit(status = 1)
} else {
    cat('All required R packages are available\n')
    quit(status = 0)
}
"

if [ $? -ne 0 ]; then
    print_error "Some required R packages are missing. Please install them first."
    print_warning "Count matrix has been created. Install R packages and run analysis separately."
    exit 1
fi

# Step 3: Run statistical analysis
print_status "Step 3: Running differential expression analysis..."

if [ ! -f "$SCRIPT_DIR/rnaseq_analysis.R" ]; then
    print_error "rnaseq_analysis.R not found in $SCRIPT_DIR"
    exit 1
fi

Rscript "$SCRIPT_DIR/rnaseq_analysis.R" \
    "$COUNT_MATRIX_DIR/gene_count_matrix.tsv" \
    "$COUNT_MATRIX_DIR/sample_metadata.tsv" \
    "$OUTPUT_DIR/statistical_analysis"

if [ $? -ne 0 ]; then
    print_error "Statistical analysis failed"
    exit 1
fi

print_success "Statistical analysis completed successfully"

# Step 4: Create final summary
print_status "Step 4: Creating final summary..."

cat > "$OUTPUT_DIR/README.md" << EOF
# RNA-seq Analysis Results

This directory contains the results of RNA-seq differential expression analysis.

## Files Generated

### Count Matrix and Metadata
- \`count_matrices/gene_count_matrix.tsv\` - Gene count matrix (genes × samples)
- \`count_matrices/sample_metadata.tsv\` - Sample information and grouping
- \`count_matrices/count_summary.txt\` - Summary statistics of count data

### Statistical Analysis Results
- \`statistical_analysis/deseq2_results.tsv\` - DESeq2 differential expression results
- \`statistical_analysis/edger_results.tsv\` - edgeR differential expression results
- \`statistical_analysis/analysis_summary.txt\` - Summary of analysis results

### Quality Control Plots
- \`statistical_analysis/qc_total_counts.pdf\` - Total read counts per sample
- \`statistical_analysis/qc_gene_counts.pdf\` - Gene count distributions
- \`statistical_analysis/heatmap_top_variable_genes.pdf\` - Heatmap of top variable genes

### Analysis Plots
- \`statistical_analysis/deseq2_ma_plot.pdf\` - DESeq2 MA plot
- \`statistical_analysis/deseq2_volcano_plot.pdf\` - DESeq2 volcano plot
- \`statistical_analysis/deseq2_vs_edger.pdf\` - Comparison of DESeq2 vs edgeR

## Next Steps

1. **Review Results**: Check the analysis summary and significant genes
2. **Functional Analysis**: Perform GO/KEGG pathway analysis on significant genes
3. **Visualization**: Create additional plots for publication
4. **Validation**: Validate key findings with qPCR or other methods

## Sample Information

Based on the sample names, your experiment appears to have:
- **DC samples**: Treatment group 1 (DC1, DC2)
- **DG samples**: Treatment group 2 (DG)
- **L/R conditions**: Different conditions or replicates
- **Replicates**: Biological replicates (1, 2, 3)

## Statistical Analysis Details

- **DESeq2**: Uses negative binomial model with shrinkage estimation
- **edgeR**: Uses negative binomial model with empirical Bayes methods
- **Filtering**: Genes with < 10 total counts were filtered out
- **Significance**: FDR < 0.05 considered significant

## Contact

For questions about this analysis, please refer to the pipeline documentation.
EOF

print_success "Analysis pipeline completed successfully!"
print_status "Results are available in: $OUTPUT_DIR"
print_status "Check README.md for detailed information about the results"

# Show quick summary
if [ -f "$OUTPUT_DIR/statistical_analysis/analysis_summary.txt" ]; then
    echo ""
    print_status "Quick Summary:"
    echo "==============="
    tail -n 20 "$OUTPUT_DIR/statistical_analysis/analysis_summary.txt"
fi

echo ""
print_success "Pipeline completed! Check $OUTPUT_DIR for all results."
