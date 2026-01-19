# PyDESeq2 Integration Guide

This guide explains how to use PyDESeq2 for differential expression analysis and CYP gene heatmap generation (Figure 6A style).

---

## Understanding the Data Flow

### Key Files Explained

| File | Type | What It Contains | When Created |
|------|------|------------------|--------------|
| `gene_count_matrix.tsv` | **INPUT** | Raw read counts per gene per sample (integers) | Built from STAR/RSEM output |
| `sample_metadata.tsv` | **INPUT** | Sample information (condition, replicate, etc.) | Built alongside count matrix |
| `pydeseq2_results.tsv` | **OUTPUT** | DE statistics (log2FC, pvalue, padj) per gene | After PyDESeq2 analysis |
| `cyp_heatmap_matrix.tsv` | **OUTPUT** | Transformed values used for heatmap plotting | After CYP heatmap generation |

### The Complete Pipeline

```
┌─────────────────────────────────────────────────────────────────────────┐
│                        RNA-seq Analysis Pipeline                         │
└─────────────────────────────────────────────────────────────────────────┘

  FASTQ files                                                              
      │                                                                    
      ▼                                                                    
┌─────────────────┐                                                        
│  STAR Alignment │  (per sample)                                          
└────────┬────────┘                                                        
         │                                                                 
         ▼                                                                 
  *_ReadsPerGene.out.tab  ←── Raw counts per sample (multiple files)       
      │                                                                    
      ▼                                                                    
┌─────────────────────────┐                                                
│  build_count_matrix.py  │  ←── Combines all samples into one matrix      
└────────────┬────────────┘                                                
             │                                                             
             ▼                                                             
┌───────────────────────────────────────┐                                  
│  gene_count_matrix.tsv  (INPUT)       │  ←── Genes × Samples matrix      
│  sample_metadata.tsv    (INPUT)       │      (raw integer counts)        
└───────────────────┬───────────────────┘                                  
                    │                                                      
                    ▼                                                      
┌───────────────────────────────────────┐                                  
│       pydeseq2_analysis.py            │  ←── Statistical analysis        
└───────────────────┬───────────────────┘                                  
                    │                                                      
        ┌───────────┴───────────┐                                          
        ▼                       ▼                                          
┌───────────────────┐  ┌────────────────────┐                              
│ pydeseq2_results  │  │   CYP Heatmap      │ (if --cyp-family-map given)  
│     .tsv          │  │   .pdf/.png        │                              
│                   │  │                    │                              
│ (DE statistics:   │  │ (Figure 6A style:  │                              
│  log2FC, padj,    │  │  root-upregulated  │                              
│  baseMean, etc.)  │  │  CYP genes)        │                              
└───────────────────┘  └────────────────────┘                              
```

### What's in Each File?

#### 1. `gene_count_matrix.tsv` (INPUT - Raw Counts)

This is the **starting point** for PyDESeq2. Contains raw, unnormalized read counts:

```
gene_id       DC1L1   DC1L2   DC1R1   DC1R2   DGL1    DGR1
LOC108175001    523     612     1847    1923    489     2103
LOC108175002    156      89      203     178     145      187
LOC108175003      0       2        1       0       0        3
...
```

- **Rows**: Gene IDs (e.g., LOC108175001)
- **Columns**: Sample names (e.g., DC1L1, DC1R1)
- **Values**: Raw integer counts (how many reads mapped to each gene)
- **Created by**: `build_count_matrix.py`

#### 2. `pydeseq2_results.tsv` (OUTPUT - DE Statistics)

This is the **result** of differential expression analysis:

```
gene_id         baseMean    log2FoldChange   lfcSE    stat      pvalue      padj
LOC108175001    1249.32          1.847       0.234    7.89    3.02e-15    1.23e-12
LOC108175002     159.67          0.312       0.156    2.00    4.55e-02    1.45e-01
LOC108196229    2847.91          4.523       0.312   14.50    2.31e-47    8.92e-44
...
```

- **baseMean**: Average normalized expression across all samples
- **log2FoldChange**: Log2 ratio (contrast_A / contrast_B). Positive = higher in A (e.g., root)
- **pvalue**: Raw p-value from statistical test
- **padj**: Adjusted p-value (FDR corrected) - **use this for significance**
- **Created by**: `pydeseq2_analysis.py`

#### 3. `cyp_heatmap_matrix.tsv` (OUTPUT - For Plotting)

Transformed expression values used to create the heatmap:

```
gene_id         DC1L1   DC1L2   DGL1    DGL2    DC1R1   DC1R2   DGR1
LOC108196229    -1.23   -1.45   -0.89   -1.02    1.34    1.56    1.12
LOC108201234    -0.67   -0.78   -0.45   -0.56    0.89    0.92    0.78
...
```

- **Values**: log2(normalized+1), then row-centered or z-scored
- **Columns**: Ordered by condition (leaf first, then root)
- **Rows**: Only CYP DEGs, ordered by family
- **Created by**: `pydeseq2_analysis.py` (when `--cyp-family-map` provided)

---

## Installation

### Option 1: Using Conda Environment

```bash
conda env update -f environment.yml
```

### Option 2: Using pip

```bash
pip install pydeseq2 matplotlib seaborn scipy
```

---

## Usage

### Step 1: Build Count Matrix (Required First!)

Before running PyDESeq2, you MUST create the count matrix from alignment output:

```bash
# From STAR alignment output
python build_count_matrix.py /path/to/star_counts/ -o count_matrices

# From Trinity/RSEM output
python build_count_matrix.py /path/to/rsem_counts/ -o count_matrices --type rsem
```

**This creates:**
- `count_matrices/gene_count_matrix.tsv` - The input for PyDESeq2
- `count_matrices/sample_metadata.tsv` - Sample information

### Step 2: Run PyDESeq2 Analysis

#### Basic DE Analysis (no CYP heatmap):

```bash
python pydeseq2_analysis.py \
    count_matrices/gene_count_matrix.tsv \
    count_matrices/sample_metadata.tsv \
    -o pydeseq2_results
```

####  CYP Heatmap:

```bash
python pydeseq2_analysis.py \
    count_matrices/gene_count_matrix.tsv \
    count_matrices/sample_metadata.tsv \
    -o pydeseq2_results \
    --contrast-factor condition \
    --contrast-A root \
    --contrast-B leaf \
    --cyp-family-map cyp_families.tsv \
    --root-up-only \
    --lfc 2.0 \
    --padj 0.05 \
    --scale center
```

### Step 3: Using SLURM (sbatch)

```bash
# Basic analysis
sbatch scripts/run_pydeseq2_analysis.sbatch

# With CYP heatmap
CYP_FAMILY_MAP=/path/to/cyp_families.tsv \
ROOT_UP_ONLY=true \
sbatch scripts/run_pydeseq2_analysis.sbatch
```

---

## Command-Line Arguments

### Required Arguments

| Argument | Description |
|----------|-------------|
| `count_matrix` | Path to gene count matrix TSV (from build_count_matrix.py) |
| `metadata` | Path to sample metadata TSV |

### Contrast Settings

| Argument | Default | Description |
|----------|---------|-------------|
| `--contrast-factor` | `condition` | Metadata column for comparison |
| `--contrast-A` | `root` | Numerator (positive log2FC means higher here) |
| `--contrast-B` | `leaf` | Denominator |
| `--design` | auto | Design formula (e.g., `condition` or `treatment+batch`) |

### DEG Filtering

| Argument | Default | Description |
|----------|---------|-------------|
| `--padj` | `0.05` | Adjusted p-value cutoff |
| `--lfc` | `2.0` | Absolute log2FC cutoff (0 to disable) |
| `--root-up-only` | `False` | Keep only genes upregulated in contrast-A |

### CYP Heatmap Options

| Argument | Default | Description |
|----------|---------|-------------|
| `--cyp-family-map` | - | TSV with gene_id, cyp_family columns |
| `--cyp-genes` | - | Simple text file with CYP gene IDs |
| `--sample-order` | - | Text file with sample names in desired order |
| `--scale` | `center` | Scaling: `center` or `zscore` |
| `--row-cluster` | `True` | Cluster genes hierarchically |
| `--col-cluster` | `False` | Cluster samples (False for Figure 6A) |

---

## Output Files

### Core DE Results

| File | Description |
|------|-------------|
| `pydeseq2_results.tsv` | Full DE results for ALL genes |
| `deg_filtered.tsv` | DEGs passing padj/lfc filters |
| `analysis_summary.txt` | Summary statistics |

### QC and Analysis Plots

| File | Description |
|------|-------------|
| `qc_total_counts.pdf` | Total reads per sample |
| `qc_gene_counts.pdf` | Gene count distributions |
| `pydeseq2_ma_plot.pdf` | MA plot (log2FC vs mean expression) |
| `pydeseq2_volcano_plot.pdf` | Volcano plot (-log10 padj vs log2FC) |

### CYP Heatmap (when `--cyp-family-map` provided)

| File | Description |
|------|-------------|
| `cyp_deg_filtered.tsv` | CYP DEGs with family annotations |
| `cyp_heatmap_matrix.tsv` | Transformed values (for custom plotting) |
| `cyp_heatmap.pdf` | High-quality heatmap |
| `cyp_heatmap.png` | Quick-view heatmap |

---

## CYP Family Extraction (Automatic)

**You don't need to create the CYP family map manually!** The pipeline automatically extracts CYP genes from your GTF file.

### How It Works

1. Set `GTF_FILE` to your annotation file
2. The sbatch script runs `extract_cyp_families.py` automatically
3. Creates `cyp_families.tsv` in the output directory

```bash
# Just provide the GTF - everything else is automatic!
GTF_FILE=/projects/tholl_lab_1/daisy_analysis/04_reference/dc_genomic.gtf \
ROOT_UP_ONLY=true \
sbatch scripts/run_pydeseq2_analysis.sbatch
```

### Manual Extraction (Optional)

If you need to run extraction separately:

```bash
python extract_cyp_families.py /path/to/genomic.gtf -o cyp_families.tsv
```

### Output Format

The auto-generated `cyp_families.tsv`:

```tsv
gene_id	cyp_family
LOC108193079	CYP71
LOC108196352	CYP86
LOC108201942	CYP707
LOC108193840	CYP72
```

### What Gets Extracted

- ✅ Cytochrome P450 genes (CYP71, CYP72, CYP86, etc.)
- ❌ Cyclophilins (peptidyl-prolyl isomerases with CYP in name) - filtered out

---

## Sample Metadata Format

Your `sample_metadata.tsv` should have:

```tsv
sample	condition	treatment	replicate
DC1L1	leaf	DC1	1
DC1L2	leaf	DC1	2
DC1R1	root	DC1	1
DC1R2	root	DC1	2
DGL1	leaf	DG	1
DGR1	root	DG	1
```

**Required:**
- `sample` - Must match column names in count matrix
- A grouping column (e.g., `condition`, `treatment`, or `group`)

---

## Example Workflows

### Workflow 1: Basic DE Analysis

```bash
# Step 1: Build count matrix
python build_count_matrix.py star_output/ -o 06_analysis/count_matrices

# Step 2: Run DE analysis
python pydeseq2_analysis.py \
    06_analysis/count_matrices/gene_count_matrix.tsv \
    06_analysis/count_matrices/sample_metadata.tsv \
    -o 06_analysis/pydeseq2_results

# Step 3: Check results
head 06_analysis/pydeseq2_results/pydeseq2_results.tsv
```

### Workflow 2: Figure 6A CYP Heatmap (Automatic)

```bash
# Step 1: Build count matrix (if not already done)
python build_count_matrix.py star_output/ -o 06_analysis/count_matrices

# Step 2: Run with GTF - CYP families extracted automatically!
GTF_FILE=/path/to/genomic.gtf \
ROOT_UP_ONLY=true \
sbatch scripts/run_pydeseq2_analysis.sbatch

# That's it! The script automatically:
# - Extracts CYP genes from GTF
# - Filters to root-upregulated DEGs
# - Creates Figure 6A-style heatmap

# Step 3: View heatmap
open 06_analysis/pydeseq2_results/cyp_heatmap.pdf
```

### Workflow 3: Manual CYP Extraction (if needed)

```bash
# Extract CYP genes from GTF manually
python extract_cyp_families.py \
    /path/to/genomic.gtf \
    -o cyp_families.tsv

# Then run with pre-made file
CYP_FAMILY_MAP=cyp_families.tsv \
ROOT_UP_ONLY=true \
sbatch scripts/run_pydeseq2_analysis.sbatch
```

---

## Troubleshooting

### "Count matrix not found"

The count matrix must be built FIRST from STAR/RSEM output:

```bash
python build_count_matrix.py /path/to/star_counts/ -o count_matrices
```

### "No CYP genes found in DEG list"

Check that:
1. Gene IDs in `--cyp-family-map` match your count matrix exactly
2. Your filters (padj, lfc) aren't too strict
3. CYP genes are actually differentially expressed

### "No samples match condition"

Check that `--contrast-A` and `--contrast-B` values match your metadata:

```bash
# See what's in your metadata
head count_matrices/sample_metadata.tsv
```

### Memory Issues

For large datasets, request more memory in SLURM:

```bash
#SBATCH --mem=64G
```

---

## Comparison: gene_count_matrix.tsv vs pydeseq2_results.tsv

| Aspect | gene_count_matrix.tsv | pydeseq2_results.tsv |
|--------|----------------------|---------------------|
| **Role** | INPUT | OUTPUT |
| **Contents** | Raw read counts | Statistical results |
| **Rows** | All genes | All tested genes |
| **Columns** | Sample names | Statistics (log2FC, padj, etc.) |
| **Values** | Integers (0, 1, 2, ...) | Floats (statistics) |
| **Created by** | build_count_matrix.py | pydeseq2_analysis.py |
| **Used for** | Starting point for DE analysis | Identifying significant genes |

**Key Point**: You need `gene_count_matrix.tsv` to create `pydeseq2_results.tsv`. The count matrix is the raw data; the results file contains the statistical analysis of that data.

---

## References

- PyDESeq2 GitHub: https://github.com/owkin/PyDESeq2
- PyDESeq2 Documentation: https://pydeseq2.readthedocs.io/
- Original DESeq2 Paper: Love et al. (2014) Genome Biology
