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
Step 1: cyp_master_list.csv + pydeseq2_results_UNFILTERED.tsv
        → cyp_expressed_list.tsv
        (only CYP genes that are in your count matrix)

Step 2: cyp_expressed_list.tsv + full protein FASTA
        → cyp_proteins.fasta
        (protein sequences for only those expressed CYP genes)

Step 3: cyp_proteins.fasta vs swissprot
        → cyp_blast_results.tsv
        (BLAST annotations for each CYP protein)

Step 4: cyp_expressed_list.tsv + cyp_blast_results.tsv
        → cyp_blast_annotated.tsv
        (expression stats + BLAST annotations in one file)

Step 5a: cyp_blast_annotated.tsv
         → cyp_blast_annotated_FILTERED.tsv
         (only significantly DE CYP genes)

Step 5b: cyp_blast_annotated_FILTERED.tsv + count matrix + metadata
         → heatmap, volcano, MA plots
```

---

## Commands to Run

### Steps 1-4: one sbatch job

```bash
cd /projects/tholl_lab_1/daisy_analysis
sbatch 05_rnaseq-code/scripts/run_cyp_heatmap.sbatch
```

Check job status:

```bash
squeue -u $USER
```

When it finishes, verify output:

```bash
ls -la 07_NRdatabase/cyp450_database/cyp_expressed_list.tsv
ls -la 07_NRdatabase/cyp450_database/cyp_proteins.fasta
ls -la 07_NRdatabase/cyp450_database/cyp_blast_results.tsv
ls -la 07_NRdatabase/cyp450_database/cyp_blast_annotated.tsv
```

### Step 5a: filter for DE genes

Uses the existing PyDESeq2 step 2 filter script. Pass the annotated file as a custom input path:

```bash
sbatch 05_rnaseq-code/scripts/run_pydeseq2_step2_filter.sbatch DC \
    /projects/tholl_lab_1/daisy_analysis/07_NRdatabase/cyp450_database/cyp_blast_annotated.tsv
```

Default filter: padj < 0.05 and |log2FC| > 2.0.

To use different cutoffs:

```bash
PADJ=0.05 LFC=1.0 sbatch 05_rnaseq-code/scripts/run_pydeseq2_step2_filter.sbatch DC \
    /projects/tholl_lab_1/daisy_analysis/07_NRdatabase/cyp450_database/cyp_blast_annotated.tsv
```

Output goes to: `06_analysis/pydeseq2_DC_step2_filtered/cyp_blast_annotated_FILTERED.tsv`

### Step 5b: generate plots (heatmap, volcano, MA)

Uses the existing PyDESeq2 step 3 plots script. Pass the filtered file:

```bash
sbatch 05_rnaseq-code/scripts/run_pydeseq2_step3_plots.sbatch DC \
    /projects/tholl_lab_1/daisy_analysis/06_analysis/pydeseq2_DC_step2_filtered/cyp_blast_annotated_FILTERED.tsv
```

Output goes to: `06_analysis/pydeseq2_DC_step3_plots_cyp_blast_annotated_*/`

---

## What Each Step Does

### Step 1 — Intersect CYP master list with PyDESeq2

**Script:** `scripts/cyp_intersect_pydeseq2.py`

**Input:**
- `cyp_master_list.csv` — 396 CYP genes (Both + HMMER_only evidence)
- `pydeseq2_results_UNFILTERED.tsv` — all ~30K genes with expression stats

**What it does:**
- Takes the 396 HMMER-confirmed CYP gene_ids (LOC*)
- Looks for them in the PyDESeq2 results
- Keeps only genes found in both (= expressed in your samples)
- Combines CYP info (evidence, description, protein_id) with expression
  stats (baseMean, log2FoldChange, padj)
- Labels each gene as root_up, leaf_up, ns, or no_data
- If any gene is missing its protein_id, fills it from the GTF

**Output:** `cyp_expressed_list.tsv`

**Columns:** gene_id, evidence, gene_name, chromosome, start, stop, strand,
description, matched_keyword, protein_id, hmmer_evalue, hmmer_score,
hmmer_bias, domain_name, domain_accession, baseMean, log2FoldChange,
pvalue, padj, direction

### Step 2 — Extract CYP protein sequences

**Script:** `scripts/cyp_extract_proteins.py`

**Input:**
- `cyp_expressed_list.tsv` — from step 1
- `GCF_001625215.2_DH1_v3.0_protein.faa` — full carrot protein FASTA (27MB, ~33K proteins)

**What it does:**
- Reads the expressed list to get each gene's protein_id (XP_*)
- Scans the full protein FASTA and pulls out only matching sequences
- Rewrites each FASTA header to use gene_id (LOC*) instead of
  protein_id, so BLAST results match directly to gene_ids downstream

**Output:** `cyp_proteins.fasta` (~396 sequences, small file)

**Example header:** `>LOC108194654 cytochrome P450 71A1-like [Daucus carota]`

### Step 3 — BLASTp (discovery mode)

**Tool:** `blastp` (runs directly in the sbatch, not a separate script)

**Input:**
- `cyp_proteins.fasta` — from step 2
- swissprot database — in `07_NRdatabase/blastdb/`

**Settings:** same as `blastp_discoveryfilter.sbatch`:
- e-value: 1e-4
- max_target_seqs: 25
- qcov_hsp_perc: 40
- outfmt 6 (tabular, 17 columns)

**Runtime:** ~5 minutes for ~396 proteins (vs hours for all 33K)

**Output:** `cyp_blast_results.tsv` (raw BLAST hits, multiple per gene)

### Step 4 — Combine BLAST + expression

**Script:** `scripts/combine_blast_deseq.py` (existing script, not modified)

**Input:**
- `cyp_expressed_list.tsv` — passed as `--deseq`
- `cyp_blast_results.tsv` — passed as `--blast`

**What it does:**
- Keeps the best BLAST hit per gene (lowest e-value)
- Parses blast_description and blast_species from the hit
- Merges BLAST annotations onto the expressed list by gene_id
- Saves one file with expression + BLAST columns together

**Output:** `cyp_blast_annotated.tsv`

### Step 5a — Filter for DE genes

**Script:** `pydeseq2_filter_results.py` via `run_pydeseq2_step2_filter.sbatch` (existing)

**Input:** `cyp_blast_annotated.tsv` (passed as custom file path)

**What it does:**
- Removes genes with NA padj
- Keeps genes with padj < 0.05 AND |log2FC| > 2.0 (adjustable)
- These are your significantly differentially expressed CYP genes

**Output:** `cyp_blast_annotated_FILTERED.tsv` (in `06_analysis/pydeseq2_DC_step2_filtered/`)

### Step 5b — Generate plots

**Script:** `run_pydeseq2_step3_plots.sbatch` (existing)

**Input:** filtered file from 5a, plus count matrix and metadata (auto-detected)

**What it generates:**
- MA plot (red/blue up/down genes)
- Volcano plot (4-color, labeled)
- PCA plot (top 500 variable genes)
- Sample correlation heatmap
- CYP heatmap (per-sample expression, root vs leaf)
- OMT heatmap (if any OMT genes are in the file)

**Output:** PDFs and PNGs in `06_analysis/pydeseq2_DC_step3_plots_*/`

---

## Output Files Summary

| File | Location | Description |
|------|----------|-------------|
| `cyp_expressed_list.tsv` | `07_NRdatabase/cyp450_database/` | CYP genes in your count matrix with expression stats |
| `cyp_proteins.fasta` | `07_NRdatabase/cyp450_database/` | Protein sequences for BLAST (gene_id headers) |
| `cyp_blast_results.tsv` | `07_NRdatabase/cyp450_database/` | Raw BLAST hits (outfmt 6) |
| `cyp_blast_annotated.tsv` | `07_NRdatabase/cyp450_database/` | Expression + BLAST in one file |
| `cyp_blast_annotated_FILTERED.tsv` | `06_analysis/pydeseq2_DC_step2_filtered/` | Only significantly DE CYP genes |
| heatmap, volcano, MA, PCA | `06_analysis/pydeseq2_DC_step3_plots_*/` | Final plots (PDF + PNG) |

---

## Sanity Checks

After steps 1-4 finish:

```bash
# How many CYP genes are in the expressed list?
# (subtract 1 for the header line)
wc -l 07_NRdatabase/cyp450_database/cyp_expressed_list.tsv

# How many proteins were extracted?
grep -c "^>" 07_NRdatabase/cyp450_database/cyp_proteins.fasta

# How many got BLAST hits?
cut -f1 07_NRdatabase/cyp450_database/cyp_blast_results.tsv | sort -u | wc -l

# Quick look at the combined annotated file
head -5 07_NRdatabase/cyp450_database/cyp_blast_annotated.tsv | column -t -s $'\t'
```

After step 5a:

```bash
# How many CYP genes are significantly DE?
wc -l 06_analysis/pydeseq2_DC_step2_filtered/cyp_blast_annotated_FILTERED.tsv
```

**Good result:** Most of the ~396 CYP genes appear in the expressed list, most
get swissprot BLAST hits with descriptions like "cytochrome P450 71A1", and
some subset (maybe 20-80) pass the DE filter.

**Bad result:**
- 0 genes in intersection → gene_id format mismatch (LOC* vs something else)
- 0 BLAST hits → swissprot DB not found (check `$BLASTDB` path)
- 0 genes pass filter → try lower cutoffs: `PADJ=0.05 LFC=1.0`

---

## Scripts Used

| Script | New or existing | Role |
|--------|----------------|------|
| `scripts/cyp_intersect_pydeseq2.py` | New | Step 1: intersect master list with PyDESeq2 |
| `scripts/cyp_extract_proteins.py` | New | Step 2: extract protein sequences |
| `blastp` | Existing tool | Step 3: BLAST against swissprot |
| `scripts/combine_blast_deseq.py` | Existing | Step 4: merge BLAST + expression |
| `pydeseq2_filter_results.py` | Existing | Step 5a: filter for DE genes |
| `scripts/run_pydeseq2_step3_plots.sbatch` | Existing | Step 5b: heatmap, volcano, MA |
| `scripts/run_cyp_heatmap.sbatch` | New | Runs steps 1-4 in one job |
