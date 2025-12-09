#!/bin/bash
# RNA-seq Count Matrix and Statistical Analysis Pipeline
# 
# This script builds a gene count matrix from STAR or RSEM output files
# and performs differential expression analysis with DESeq2/edgeR (R) or PyDESeq2 (Python).
#
# Usage: ./run_rnaseq_analysis.sh [count_dir] [output_dir] [--method METHOD]
#   METHOD: 'r' (DESeq2/edgeR), 'python' (PyDESeq2), or 'both' (default: 'both')

set -e  # Exit on any error

# Default values
COUNT_DIR=""
OUTPUT_DIR="rnaseq_analysis_results"
METHOD="both"  # 'r', 'python', or 'both'
COUNT_TYPE="auto"  # 'star', 'rsem', or 'auto'
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
    echo "Usage: $0 <count_dir> [output_dir] [--method METHOD] [--type TYPE]"
    echo ""
    echo "Arguments:"
    echo "  count_dir       Directory containing count files (STAR or RSEM)"
    echo "  output_dir      Output directory (default: rnaseq_analysis_results)"
    echo ""
    echo "Options:"
    echo "  --method METHOD Analysis method: 'r' (DESeq2/edgeR), 'python' (PyDESeq2), or 'both' (default)"
    echo "  --type TYPE      Count file type: 'star', 'rsem', or 'auto' (default: auto-detect)"
    echo ""
    echo "Examples:"
    echo "  $0 /path/to/star/counts"
    echo "  $0 /path/to/rsem/counts --type rsem --method python"
    echo "  $0 /path/to/counts my_results --method both"
    echo ""
    echo "The script will:"
    echo "  1. Build gene count matrix from STAR or RSEM output"
    echo "  2. Create sample metadata"
    echo "  3. Run differential expression analysis (R and/or Python)"
    echo "  4. Generate quality control plots"
    echo "  5. Create summary reports"
}

# Parse command line arguments
if [ $# -eq 0 ]; then
    show_usage
    exit 1
fi

COUNT_DIR="$1"
shift

# Parse optional arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --method)
            METHOD="$2"
            shift 2
            ;;
        --type)
            COUNT_TYPE="$2"
            shift 2
            ;;
        *)
            if [ -z "$OUTPUT_DIR" ] || [ "$OUTPUT_DIR" = "rnaseq_analysis_results" ]; then
                OUTPUT_DIR="$1"
            fi
            shift
            ;;
    esac
done

# Validate input directory
if [ ! -d "$COUNT_DIR" ]; then
    print_error "Directory does not exist: $COUNT_DIR"
    exit 1
fi

# Validate method
if [[ ! "$METHOD" =~ ^(r|python|both)$ ]]; then
    print_error "Invalid method: $METHOD. Must be 'r', 'python', or 'both'"
    exit 1
fi

# Check for count files
STAR_FILES=$(find "$COUNT_DIR" -name "*ReadsPerGene.out.tab" 2>/dev/null | wc -l)
RSEM_FILES=$(find "$COUNT_DIR" -name "*.genes.results" 2>/dev/null | wc -l)

if [ "$COUNT_TYPE" = "auto" ]; then
    if [ "$STAR_FILES" -gt 0 ]; then
        COUNT_TYPE="star"
        COUNT_FILES=$STAR_FILES
    elif [ "$RSEM_FILES" -gt 0 ]; then
        COUNT_TYPE="rsem"
        COUNT_FILES=$RSEM_FILES
    else
        print_error "No count files found in $COUNT_DIR"
        exit 1
    fi
elif [ "$COUNT_TYPE" = "star" ]; then
    COUNT_FILES=$STAR_FILES
    if [ "$COUNT_FILES" -eq 0 ]; then
        print_error "No ReadsPerGene.out.tab files found in $COUNT_DIR"
        exit 1
    fi
elif [ "$COUNT_TYPE" = "rsem" ]; then
    COUNT_FILES=$RSEM_FILES
    if [ "$COUNT_FILES" -eq 0 ]; then
        print_error "No .genes.results files found in $COUNT_DIR"
        exit 1
    fi
fi

print_status "Found $COUNT_FILES $COUNT_TYPE count files in $COUNT_DIR"

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

python3 "$SCRIPT_DIR/build_count_matrix.py" "$COUNT_DIR" -o "$COUNT_MATRIX_DIR" --type "$COUNT_TYPE"

if [ $? -ne 0 ]; then
    print_error "Failed to build count matrix"
    exit 1
fi

print_success "Count matrix created successfully"

# Step 2: Run statistical analysis
if [[ "$METHOD" =~ ^(r|both)$ ]]; then
    print_status "Step 2a: Running R-based analysis (DESeq2/edgeR)..."
    
    if ! command -v Rscript &> /dev/null; then
        print_warning "Rscript not found. Skipping R-based analysis."
        if [ "$METHOD" = "r" ]; then
            print_error "R method requested but Rscript not available."
            exit 1
        fi
    else
        # Check for required R packages
        print_status "Checking required R packages..."
        
        Rscript -e "
        required_packages <- c('DESeq2', 'edgeR', 'ggplot2', 'pheatmap', 'RColorBrewer', 'dplyr', 'EnhancedVolcano')
        missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
        
        if (length(missing_packages) > 0) {
            cat('Missing R packages:', paste(missing_packages, collapse = ', '), '\n')
            quit(status = 1)
        } else {
            cat('All required R packages are available\n')
            quit(status = 0)
        }
        " 2>/dev/null
        
        if [ $? -eq 0 ]; then
            if [ ! -f "$SCRIPT_DIR/rnaseq_analysis.R" ]; then
                print_error "rnaseq_analysis.R not found in $SCRIPT_DIR"
            else
                Rscript "$SCRIPT_DIR/rnaseq_analysis.R" \
                    "$COUNT_MATRIX_DIR/gene_count_matrix.tsv" \
                    "$COUNT_MATRIX_DIR/sample_metadata.tsv" \
                    "$OUTPUT_DIR/statistical_analysis_r"
                
                if [ $? -eq 0 ]; then
                    print_success "R-based analysis completed successfully"
                else
                    print_warning "R-based analysis failed, but continuing..."
                fi
            fi
        else
            print_warning "Some required R packages are missing. Skipping R-based analysis."
        fi
    fi
fi

if [[ "$METHOD" =~ ^(python|both)$ ]]; then
    print_status "Step 2b: Running Python-based analysis (PyDESeq2)..."
    
    if ! command -v python3 &> /dev/null; then
        print_error "python3 not found. Please install Python 3."
        if [ "$METHOD" = "python" ]; then
            exit 1
        fi
    else
        # Check if PyDESeq2 is installed
        python3 -c "import pydeseq2" 2>/dev/null
        if [ $? -ne 0 ]; then
            print_warning "PyDESeq2 not found. Install with: pip install pydeseq2"
            if [ "$METHOD" = "python" ]; then
                print_error "Python method requested but PyDESeq2 not available."
                exit 1
            fi
        else
            if [ ! -f "$SCRIPT_DIR/pydeseq2_analysis.py" ]; then
                print_error "pydeseq2_analysis.py not found in $SCRIPT_DIR"
            else
                python3 "$SCRIPT_DIR/pydeseq2_analysis.py" \
                    "$COUNT_MATRIX_DIR/gene_count_matrix.tsv" \
                    "$COUNT_MATRIX_DIR/sample_metadata.tsv" \
                    -o "$OUTPUT_DIR/statistical_analysis_python"
                
                if [ $? -eq 0 ]; then
                    print_success "Python-based analysis completed successfully"
                else
                    print_warning "Python-based analysis failed, but continuing..."
                fi
            fi
        fi
    fi
fi

# Step 3: Create final summary
print_status "Step 3: Creating final summary..."

cat > "$OUTPUT_DIR/README.md" << EOF
# RNA-seq Analysis Results

This directory contains the results of RNA-seq differential expression analysis.

## Files Generated

### Count Matrix and Metadata
- \`count_matrices/gene_count_matrix.tsv\` - Gene count matrix (genes × samples)
- \`count_matrices/sample_metadata.tsv\` - Sample information and grouping
- \`count_matrices/count_summary.txt\` - Summary statistics of count data

### Statistical Analysis Results

#### R-based Analysis (DESeq2/edgeR)
- \`statistical_analysis_r/deseq2_results.tsv\` - DESeq2 differential expression results
- \`statistical_analysis_r/edger_results.tsv\` - edgeR differential expression results
- \`statistical_analysis_r/analysis_summary.txt\` - Summary of R analysis results

#### Python-based Analysis (PyDESeq2)
- \`statistical_analysis_python/pydeseq2_results.tsv\` - PyDESeq2 differential expression results
- \`statistical_analysis_python/analysis_summary.txt\` - Summary of Python analysis results

### Quality Control Plots
- \`statistical_analysis_*/qc_total_counts.pdf\` - Total read counts per sample
- \`statistical_analysis_*/qc_gene_counts.pdf\` - Gene count distributions
- \`statistical_analysis_*/heatmap_top_variable_genes.pdf\` - Heatmap of top variable genes

### Analysis Plots
- \`statistical_analysis_*/deseq2_ma_plot.pdf\` or \`pydeseq2_ma_plot.pdf\` - MA plot
- \`statistical_analysis_*/deseq2_volcano_plot.pdf\` or \`pydeseq2_volcano_plot.pdf\` - Volcano plot

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

### DESeq2 (R)
- Uses negative binomial model with shrinkage estimation
- Filtering: Genes with < 10 total counts were filtered out
- Significance: FDR < 0.05 considered significant

### edgeR (R)
- Uses negative binomial model with empirical Bayes methods
- Filtering: Genes with < 10 total counts were filtered out
- Significance: FDR < 0.05 considered significant

### PyDESeq2 (Python)
- Python implementation of DESeq2 method
- Filtering: Genes with < 10 total counts were filtered out
- Significance: padj < 0.05 considered significant

## Contact

For questions about this analysis, please refer to the pipeline documentation.
EOF

print_success "Analysis pipeline completed successfully!"
print_status "Results are available in: $OUTPUT_DIR"
print_status "Check README.md for detailed information about the results"

# Show quick summary
if [ -f "$OUTPUT_DIR/statistical_analysis_r/analysis_summary.txt" ]; then
    echo ""
    print_status "R Analysis Summary:"
    echo "==================="
    tail -n 15 "$OUTPUT_DIR/statistical_analysis_r/analysis_summary.txt"
fi

if [ -f "$OUTPUT_DIR/statistical_analysis_python/analysis_summary.txt" ]; then
    echo ""
    print_status "Python Analysis Summary:"
    echo "========================="
    tail -n 15 "$OUTPUT_DIR/statistical_analysis_python/analysis_summary.txt"
fi

echo ""
print_success "Pipeline completed! Check $OUTPUT_DIR for all results."
