# PyDESeq2 Integration Guide

This guide explains how to use PyDESeq2 for differential expression analysis in the RNA-seq pipeline.

## Overview

PyDESeq2 is a Python implementation of the DESeq2 method for bulk RNA-seq differential expression analysis. It provides similar functionality to the R-based DESeq2 but allows you to stay within the Python ecosystem.

## Installation

### Option 1: Using Conda Environment

Update your conda environment:

```bash
conda env update -f environment.yml
```

This will install PyDESeq2 and required dependencies (matplotlib, seaborn, scipy, scikit-learn).

### Option 2: Using pip

```bash
pip install pydeseq2 matplotlib seaborn
```

## Usage

### Step 1: Build Count Matrix

First, build a count matrix from your STAR or RSEM output files:

```bash
# For STAR counts (default)
python build_count_matrix.py star_counts/ -o count_matrices

# For Trinity/RSEM counts
python build_count_matrix.py trinity_counts/ -o count_matrices --type trinity

# Auto-detect (prefers STAR if both found)
python build_count_matrix.py counts/ -o count_matrices --type auto
```

This creates:
- `count_matrices/gene_count_matrix.tsv` - Gene count matrix
- `count_matrices/sample_metadata.tsv` - Sample metadata with grouping information

### Step 2: Run PyDESeq2 Analysis

#### Option A: Standalone Script

```bash
python pydeseq2_analysis.py \
    count_matrices/gene_count_matrix.tsv \
    count_matrices/sample_metadata.tsv \
    -o pydeseq2_results
```

#### Option B: Integrated Pipeline Script

The `run_rnaseq_analysis.sh` script now supports PyDESeq2:

```bash
# Run only PyDESeq2
./run_rnaseq_analysis.sh counts/ results/ --method python

# Run both R and Python analyses
./run_rnaseq_analysis.sh counts/ results/ --method both

# Specify count type
./run_rnaseq_analysis.sh trinity_counts/ results/ --type trinity --method python
```

### Step 3: Review Results

The analysis generates:

1. **Results Table**: `pydeseq2_results/pydeseq2_results.tsv`
   - Contains log2FoldChange, pvalue, padj, and other statistics
   - Similar format to DESeq2 results

2. **Quality Control Plots**:
   - `qc_total_counts.pdf` - Total read counts per sample
   - `qc_gene_counts.pdf` - Gene count distributions

3. **Analysis Plots**:
   - `pydeseq2_ma_plot.pdf` - MA plot showing fold changes
   - `pydeseq2_volcano_plot.pdf` - Volcano plot for significance
   - `heatmap_top_variable_genes.pdf` - Heatmap of top 50 variable genes

4. **Summary Report**: `analysis_summary.txt`

## Design Formula

The script automatically detects the design formula from your metadata:

- **Preferred**: `treatment` column
- **Alternative**: `group` column
- **Fallback**: `condition` column
- **Custom**: Use `--design` flag to specify

Example with custom design:

```bash
python pydeseq2_analysis.py \
    count_matrices/gene_count_matrix.tsv \
    count_matrices/sample_metadata.tsv \
    -o results \
    --design "treatment+condition"
```

## Sample Metadata Format

Your `sample_metadata.tsv` should have at least:
- `sample` column with sample names matching count matrix columns
- A grouping column (`treatment`, `group`, or `condition`)

Example:

```tsv
sample  treatment    condition   replicate
DC1L1   DC1         L           1
DC1L2   DC1         L           2
DC1R1   DC1         R           1
DGL1    DG          L           1
```

## Comparison with R DESeq2

| Feature | R DESeq2 | PyDESeq2 |
|---------|----------|----------|
| Method | Original DESeq2 | Python re-implementation |
| Dependencies | R + Bioconductor | Python only |
| Results | Highly validated | Similar but may have minor differences |
| Features | Full feature set | Core features (expanding) |
| Performance | Optimized | Good performance |

**Note**: PyDESeq2 is a re-implementation, so results may have minor numerical differences from R DESeq2, but should be very similar.

## Troubleshooting

### Import Errors

If you get import errors, ensure PyDESeq2 is installed:

```bash
python -c "import pydeseq2; print(pydeseq2.__version__)"
```

### Memory Issues

For large datasets, you may need to increase memory. The script uses `n_cpus=1` by default. You can modify the script to use more CPUs if needed.

### Design Formula Errors

If you get errors about the design formula:
1. Check that your metadata has the expected columns
2. Use `--design` to explicitly specify the formula
3. Ensure factor levels are appropriate (at least 2 levels)

### Count Matrix Issues

- Ensure sample names in count matrix match metadata `sample` column
- Check that count matrix has integer values (not TPM/FPKM)
- Verify no missing values or infinite values

## Integration with STAR Pipeline

After running STAR alignment:

```bash
# 1. Build count matrix from STAR output
python build_count_matrix.py star_alignments/ -o count_matrices

# 2. Run PyDESeq2 analysis
python pydeseq2_analysis.py \
    count_matrices/gene_count_matrix.tsv \
    count_matrices/sample_metadata.tsv \
    -o pydeseq2_results
```

## Integration with Trinity/RSEM Pipeline

After running Trinity + RSEM:

```bash
# 1. Build count matrix from Trinity/RSEM output
python build_count_matrix.py trinity_results/ -o count_matrices --type trinity

# 2. Run PyDESeq2 analysis
python pydeseq2_analysis.py \
    count_matrices/gene_count_matrix.tsv \
    count_matrices/sample_metadata.tsv \
    -o pydeseq2_results
```

## Next Steps

After differential expression analysis:

1. **Filter Results**: Focus on genes with padj < 0.05 and meaningful fold changes
2. **Functional Analysis**: Use significant genes for GO/KEGG pathway analysis
3. **Visualization**: Create publication-quality figures
4. **Validation**: Validate key findings with qPCR or other methods

## References

- PyDESeq2 GitHub: https://github.com/owkin/PyDESeq2
- PyDESeq2 Documentation: https://pydeseq2.readthedocs.io/
- Original DESeq2 Paper: Love et al. (2014) Genome Biology

