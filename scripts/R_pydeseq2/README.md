# R_pydeseq2 — Student R DESeq2 Workflow Replicated in Python

A single Python script (`R_pydeq2.py`) that replicates all 14 steps of the
previous student's R-based RNA-seq analysis using PyDESeq2 and Python plotting
libraries.

---

## Quick Start

```bash
# On the HPC — default settings (DC species, Leaf vs Root — Root is baseline):
sbatch scripts/R_pydeseq2/run_R_pydeq2.sbatch DC

# Custom contrast (Ethylene vs Normal — Normal is baseline):
CONTRAST_A=E CONTRAST_B=N sbatch scripts/R_pydeseq2/run_R_pydeq2.sbatch DC

# Skip slow steps:
SKIP_PHYLO=1 SKIP_NCBI=1 sbatch scripts/R_pydeseq2/run_R_pydeq2.sbatch DC
```

---

## What Input Files Do You Need?

Check that these files exist **before** running. The table shows where each
file comes from and which upstream script creates it.

| File | How it was made | HPC path |
|------|----------------|----------|
| `gene_count_matrix.tsv` | `sbatch scripts/featurecounts.sbatch DC` | `03_count_tables/00_1_DC/` |
| `sample_metadata.tsv` | Created by `build_count_matrix.py` (inside featurecounts.sbatch) | `03_count_tables/00_1_DC/` |
| `P450_list_RefSeq.txt` | Previous student (Geneious curation) | `07_NRdatabase/sukman_database/` |
| `dc_genomic.gtf` | NCBI reference download | `04_reference/` |
| `GCF_001625215.2_DH1_v3.0_protein.faa` | NCBI reference download | `04_reference/` |
| Methyltransferase list (optional) | Your curation or HMMER scan | `07_NRdatabase/omt_database/` |

### Prerequisite order (if starting from scratch)

```
1. featurecounts.sbatch DC
   └─→ gene_count_matrix.tsv + sample_metadata.tsv

2. (Optional) Reference files already on HPC:
   └─→ dc_genomic.gtf, GCF_..._protein.faa

3. (Optional) Gene lists from previous student or your HMMER pipeline:
   └─→ P450_list_RefSeq.txt, cyp_master_list.csv, etc.
```

If `gene_count_matrix.tsv` and `sample_metadata.tsv` exist, you can run the
script immediately — P450/MT filtering and phylo steps are skipped gracefully
if their input lists don't exist.

---

## What Does Each Step Do?

| Step | R equivalent | Python output file |
|------|-------------|-------------------|
| 0-1 | `read.delim()` + `rowSums > 20` filter | (loaded in memory) |
| 2-4 | `DESeqDataSetFromMatrix()` + `DESeq()` + `counts(dds, normalized=TRUE)` | `normalized_counts_refseq.csv`, `res_all.csv` |
| 5a | `vst()` + `plotPCA()` | `pca_plot.png/pdf` |
| 5b | `plotDispEsts()` | `dispersion_plot.png/pdf` |
| 7 | `subset(res, padj < 0.05 & abs(log2FC) >= 2)` | `res_all_significant.csv`, `res_Upregulated.csv`, `res_Downregulated.csv` |
| 8 | `pheatmap(log2(sigCounts+1), scale='row')` | `heatmap_all_significant_genes.png/pdf` |
| 9a | `EnhancedVolcano()` | `volcano_plot.png/pdf` |
| 9b | `plotMA(res)` | `ma_plot.png/pdf` |
| 10 | `res.df[gene_id %in% P450_list]` → `pheatmap` | `heatmap_P450.png/pdf`, `P450_expression_refseq.csv` |
| 11 | Same for methyltransferases | `heatmap_MT.png/pdf`, `MT_expression_refseq.csv` |
| 12 | `read.tree()` + `ggtree()` + `gheatmap()` | `heatmap_phylogeny_P450.png/pdf` |
| 13 | `rentrez::entrez_summary()` | `P450_genomic_locations_description.csv` |
| 14 | Distance calculation + cluster plots | `Genes_with_distance_P450.csv`, `P450_distance_*.png` |

---

## All Configurable Parameters

Set as environment variables before `sbatch`:

| Variable | Default | Description |
|----------|---------|-------------|
| `CONTRAST_A` | `L` | Numerator condition (positive log2FC = higher in Leaf) |
| `CONTRAST_B` | `R` | Denominator/baseline condition (Root) |
| `CONTRAST_FACTOR` | `condition` | Metadata column name |
| `MIN_COUNTS` | `20` | Minimum total counts to keep a gene |
| `PADJ` | `0.05` | Adjusted p-value cutoff |
| `LFC` | `2.0` | Log2 fold change cutoff |
| `THRESHOLD` | `50000` | Genomic cluster distance (bp) |
| `P450_LIST` | `P450_list_RefSeq.txt` | Path to P450 gene list |
| `MT_LIST` | (empty) | Path to methyltransferase gene list |
| `SKIP_PHYLO` | `0` | Set to `1` to skip tree building |
| `SKIP_NCBI` | `0` | Set to `1` to skip NCBI queries |
| `EMAIL` | `daisycortesj@vt.edu` | Email for NCBI Entrez |

---

## Required Python Packages

Listed in `requirements.txt`. Install with:

```bash
pip install -r scripts/R_pydeseq2/requirements.txt
```

| Package | What it replaces from R |
|---------|------------------------|
| `pydeseq2` | `DESeq2` |
| `numpy` + `pandas` | Base R data frames |
| `matplotlib` + `seaborn` | `pheatmap`, `ggplot2`, `EnhancedVolcano` |
| `scipy` | Hierarchical clustering |
| `biopython` (optional) | `rentrez`, `ape::read.tree` |

External tools (optional, for phylo tree building):
- `mafft` — multiple sequence alignment
- `FastTree` — phylogenetic tree inference

On the HPC these are loaded by the sbatch wrapper:
```bash
module load MAFFT/7.505-GCC-11.3.0-with-extensions
module load FastTree/2.1.11-GCCcore-11.3.0
```

---

## Output Directory

Results go to `06_analysis/R_pydeseq2_{SPECIES}/`

Download to your laptop:
```bash
scp -r daisycortesj@tinkercliffs.arc.vt.edu:/projects/tholl_lab_1/daisy_analysis/06_analysis/R_pydeseq2_DC/ ~/Desktop/
```

---

## Phylo Trees + Genomic Clustering Are Built In

Steps 12-14 of `R_pydeq2.py` already include phylogenetic tree building,
NCBI gene info fetching, and genomic clustering — no separate script needed.

To control these steps:

```bash
# Run everything including phylo + NCBI + genomic (default):
sbatch scripts/R_pydeseq2/run_R_pydeq2.sbatch DC

# Skip phylo tree building (faster, skips MAFFT/FastTree):
SKIP_PHYLO=1 sbatch scripts/R_pydeseq2/run_R_pydeq2.sbatch DC

# Skip NCBI queries (no internet needed):
SKIP_NCBI=1 sbatch scripts/R_pydeseq2/run_R_pydeq2.sbatch DC

# Change genomic cluster distance (default 50 kb):
THRESHOLD=100000 sbatch scripts/R_pydeseq2/run_R_pydeq2.sbatch DC
```

---

## Running With Your Existing Pipeline (Starting From Count Matrix)

If you already have `gene_count_matrix.tsv` from `featurecounts.sbatch`,
you only need one command:

```bash
sbatch scripts/R_pydeseq2/run_R_pydeq2.sbatch DC
```

This runs all 14 steps in a single job — DESeq2, QC plots, heatmaps,
volcano, P450 filtering, phylo tree, NCBI fetch, and genomic clustering.

### Using phylo_heatmap.py and genomic_cluster_analysis.py Independently

These standalone scripts are used by the **CYP heatmap pipeline** (not
R_pydeq2). They take the output from `run_cyp_express_extract.sbatch`.
See `ReadME/CYP_HEATMAP_PIPELINE.md` for full run order, or use:

```bash
# Run both phylo + genomic as post-analysis for the CYP pipeline:
sbatch scripts/post_pydeq2.sbatch DC cyp

# Or individually:
sbatch scripts/run_phylo_heatmap.sbatch DC cyp
sbatch scripts/run_genomic_clustering.sbatch DC cyp
```

---

## Files in This Directory

| File | Purpose |
|------|---------|
| `R_pydeq2.py` | Main Python script — all 14 steps (including phylo + genomic) |
| `run_R_pydeq2.sbatch` | SLURM wrapper for HPC submission |
| `requirements.txt` | Python package dependencies |
| `README.md` | This file |
