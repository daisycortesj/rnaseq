# CYP450 Heatmap Pipeline

Focused CYP-only analysis: take the HMMER-confirmed CYP master list, keep only
genes that are expressed in your samples, BLAST them against swissprot, filter
for differentially expressed genes, and generate heatmaps/plots.

All commands run on the ARC HPC from `/projects/tholl_lab_1/daisy_analysis`.

---

## Prerequisites (already done, do not rerun)

| File | How it was made | Location |
|------|----------------|----------|
| `cyp_master_list.csv` | `sbatch scripts/run_cyp450_database.sbatch` | `07_NRdatabase/cyp450_database/` |
| `pydeseq2_results_UNFILTERED.tsv` | `sbatch scripts/run_pydeseq2_step1_analysis.sbatch DC` | `06_analysis/pydeseq2_DC_step1_unfiltered/` |
| `GCF_001625215.2_DH1_v3.0_protein.faa` | NCBI download | `04_reference/` |
| `gene_count_matrix.tsv` | `sbatch scripts/build_count_matrix.sbatch` | `03_count_tables/00_1_DC/` |
| `sample_metadata.tsv` | Created manually | `03_count_tables/00_1_DC/` |
| swissprot BLAST database | `sbatch scripts/download_blastdb.sbatch` | `07_NRdatabase/blastdb/` |

---

## Pipeline Overview

```
Step 1+2: run_cyp_express_extract.sbatch
          → cyp_expressed_list.tsv (CYP genes in your count matrix)
          → cyp_proteins.fasta (protein sequences for those genes)

Step 3:   blastp_discoveryfilter.sbatch DC swissprot cyp_proteins.fasta
          → BLAST results for CYP proteins

Step 4:   run_combine_filter.sbatch DC swissprot discovery standard cyp
          → cyp_blast_annotated.tsv (expression + BLAST combined)
          → cyp_blast_filtered_standard.tsv (only significant DE CYPs)

Step 5:   run_pydeseq2_step3_plots.sbatch
          → heatmap, volcano, MA plots
```

---

## Commands to Run

### Step 1+2: Intersect + extract proteins

```bash
cd /projects/tholl_lab_1/daisy_analysis
sbatch 05_rnaseq-code/scripts/run_cyp_express_extract.sbatch
```

Wait for it to finish. Check output:

```bash
ls -la 07_NRdatabase/cyp450_database/cyp_expressed_list.tsv
ls -la 07_NRdatabase/cyp450_database/cyp_proteins.fasta
grep -c "^>" 07_NRdatabase/cyp450_database/cyp_proteins.fasta
```

### Step 3: BLAST CYP proteins against swissprot

Uses the existing `blastp_discoveryfilter.sbatch` with a custom query FASTA
(3rd argument):

```bash
sbatch 05_rnaseq-code/scripts/blastp_discoveryfilter.sbatch DC swissprot \
    /projects/tholl_lab_1/daisy_analysis/07_NRdatabase/cyp450_database/cyp_proteins.fasta
```

This is fast (~5 min for ~396 proteins vs hours for all 33K).

Output goes to: `06_analysis/blastp_DC_cyp_proteins/`

Check when done:

```bash
ls -la 06_analysis/blastp_DC_cyp_proteins/
wc -l 06_analysis/blastp_DC_cyp_proteins/blastp_DC_swissprot_discovery.tsv
```

### Step 4: Combine BLAST + expression, then filter

Uses the existing `run_combine_filter.sbatch` with `cyp` as the 5th argument.
This tells it to use the CYP expressed list and CYP BLAST results, then
filters for significant DE genes:

```bash
sbatch 05_rnaseq-code/scripts/run_combine_filter.sbatch DC swissprot discovery standard cyp
```

To change the filter cutoff (e.g. lenient):

```bash
sbatch 05_rnaseq-code/scripts/run_combine_filter.sbatch DC swissprot discovery lenient cyp
```

Output (in `07_NRdatabase/cyp450_database/`):
- `cyp_blast_annotated.tsv` — all CYP genes with expression + BLAST
- `cyp_blast_filtered_standard.tsv` — only significant DE CYP genes

Check when done:

```bash
wc -l 07_NRdatabase/cyp450_database/cyp_blast_annotated.tsv
wc -l 07_NRdatabase/cyp450_database/cyp_blast_filtered_standard.tsv
head -5 07_NRdatabase/cyp450_database/cyp_blast_annotated.tsv | column -t -s $'\t'
```

### Step 5: Generate plots (heatmap, volcano, MA)

Uses the existing step 3 plots script. Pass the filtered file:

```bash
sbatch 05_rnaseq-code/scripts/run_pydeseq2_step3_plots.sbatch DC \
    /projects/tholl_lab_1/daisy_analysis/07_NRdatabase/cyp450_database/cyp_blast_filtered_standard.tsv
```

Output: `06_analysis/pydeseq2_DC_step3_plots_cyp_blast_filtered_*/`

---

## What Each Step Does

### Steps 1+2 — Intersect + extract (`run_cyp_express_extract.sbatch`)

**Step 1** takes the 396 HMMER-confirmed CYP gene_ids from `cyp_master_list.csv`
and looks for them in the PyDESeq2 unfiltered results. Keeps only CYP genes
that are in your count matrix (= expressed in your samples). Attaches expression
stats (baseMean, log2FoldChange, padj) to each gene.

**Step 2** reads the expressed list, looks up each gene's protein_id (XP_*),
and extracts those protein sequences from the full carrot protein FASTA.
Writes a small subset FASTA (~396 sequences) with gene_id (LOC*) as headers.

**Output:**
- `cyp_expressed_list.tsv` — CYP genes with expression stats + protein_id
- `cyp_proteins.fasta` — protein sequences for BLAST

### Step 3 — BLAST (`blastp_discoveryfilter.sbatch`)

Same BLAST script you already use for all proteins, but now with a 3rd argument
pointing to the CYP subset FASTA. Same settings (e-value 1e-4, discovery mode).
Only ~396 proteins instead of 33K, so it runs in minutes.

**Output:** `blastp_DC_swissprot_discovery.tsv` in `06_analysis/blastp_DC_cyp_proteins/`

### Step 4 — Combine + filter (`run_combine_filter.sbatch` with `cyp`)

Same sbatch you already use, but with `cyp` as the 5th argument. This tells
it to use the CYP expressed list as the DESeq input and the CYP BLAST results.
It does two things in one job:
1. Merges the best BLAST hit per gene onto the expressed list
2. Filters for significant DE genes (padj + log2FC cutoffs)

**Output:**
- `cyp_blast_annotated.tsv` — all CYP genes with expression + BLAST
- `cyp_blast_filtered_standard.tsv` — only significant DE CYP genes

### Step 5 — Plots (`run_pydeseq2_step3_plots.sbatch`)

Same plots script you already use. Generates heatmap, volcano plot, MA plot,
PCA, and sample correlation heatmap from the filtered CYP gene list.

**Output:** PDFs and PNGs

---

## Output Files Summary

| File | Location | Description |
|------|----------|-------------|
| `cyp_expressed_list.tsv` | `07_NRdatabase/cyp450_database/` | CYP genes in count matrix with expression stats |
| `cyp_proteins.fasta` | `07_NRdatabase/cyp450_database/` | Protein sequences for BLAST (gene_id headers) |
| `blastp_DC_swissprot_discovery.tsv` | `06_analysis/blastp_DC_cyp_proteins/` | Raw BLAST hits |
| `cyp_blast_annotated.tsv` | `07_NRdatabase/cyp450_database/` | Expression + BLAST combined |
| `cyp_blast_filtered_standard.tsv` | `07_NRdatabase/cyp450_database/` | Only significant DE CYP genes |
| heatmap, volcano, MA, PCA | `06_analysis/pydeseq2_DC_step3_plots_*/` | Final plots (PDF + PNG) |

---

## Sanity Checks

After steps 1+2:

```bash
wc -l 07_NRdatabase/cyp450_database/cyp_expressed_list.tsv
grep -c "^>" 07_NRdatabase/cyp450_database/cyp_proteins.fasta
```

After step 3:

```bash
cut -f1 06_analysis/blastp_DC_cyp_proteins/blastp_DC_swissprot_discovery.tsv | sort -u | wc -l
```

After step 4:

```bash
wc -l 07_NRdatabase/cyp450_database/cyp_blast_annotated.tsv
wc -l 07_NRdatabase/cyp450_database/cyp_blast_filtered_standard.tsv
```

**Good result:** Most ~396 CYP genes in expressed list, most get BLAST hits,
some subset (20-80) pass the DE filter.

**Bad result:**
- 0 genes in intersection → gene_id format mismatch
- 0 BLAST hits → check `$BLASTDB` path
- 0 genes pass filter → try lenient mode

---

## Scripts Used

| Script | New or existing | Role |
|--------|----------------|------|
| `scripts/run_cyp_express_extract.sbatch` | New | Steps 1+2: intersect + extract |
| `scripts/cyp_intersect_pydeseq2.py` | New | Called by step 1 |
| `scripts/cyp_extract_proteins.py` | New | Called by step 2 |
| `scripts/blastp_discoveryfilter.sbatch` | Existing (updated) | Step 3: BLAST with custom query |
| `scripts/run_combine_filter.sbatch` | Existing (updated) | Step 4: combine + filter with `cyp` source |
| `scripts/run_pydeseq2_step3_plots.sbatch` | Existing | Step 5: heatmap, volcano, MA |
