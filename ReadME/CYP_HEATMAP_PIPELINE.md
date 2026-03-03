# CYP450 Heatmap Pipeline

Three paths to go from a CYP gene list to heatmaps/plots. Pick the one that
fits your goal — all produce output compatible with the step 3 plotting script.

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
| Previous student's CYP list (optional) | Geneious list | `07_NRdatabase/sukman_database/P450_list_RefSeq.txt` |

---

## Pipeline Overview — Three Paths

```
═══ PATH A: DIRECT (gene list → PyDESeq2 → filter → plots) ════════════════

  run_filter_count_genelist.sbatch
    → Runs PyDESeq2 on full count matrix, filters to your gene list,
      applies DE cutoffs (padj < 0.05, |log2FC| > 2.0)
    → geneious_candidates.tsv (or cyp_candidates.tsv)

  run_pydeseq2_step3_plots.sbatch DC <output.tsv>
    → heatmap, volcano, MA plots

  SCRIPTS: filter_count_by_genelist.py
  BEST FOR: Geneious P450 list, quick exploration, any gene list

═══ PATH B: SHORT (HMMER list → intersect → plots) ════════════════════════

  run_cyp_express_extract.sbatch
    → cyp_expressed_list.tsv  (YOUR HMMER/GTF list ∩ PyDESeq2 step 1)
    → geneious_expressed_list.tsv  (Geneious list, if PREV_LIST set)

  run_pydeseq2_step3_plots.sbatch DC cyp_expressed
    → heatmap, volcano, MA plots
    (plots script filters internally: padj < 0.05, |log2FC| > 2.0)

  SCRIPTS: cyp_intersect_pydeseq2.py, cyp_extract_proteins.py
  BEST FOR: Your HMMER/GTF CYP master list, protein extraction for BLAST

═══ PATH C: FULL (HMMER list → BLAST → combine → plots) ═══════════════════

  run_cyp_express_extract.sbatch  (same as Path B step 1)
  blastp_discoveryfilter.sbatch DC swissprot cyp_proteins.fasta
  run_combine_filter.sbatch DC swissprot discovery standard cyp
  run_pydeseq2_step3_plots.sbatch DC cyp_discovery
    → heatmap, volcano, MA plots (with swissprot annotations on labels)

  BEST FOR: Publication-quality figures with protein names
```

**Which path to choose?**
- **Path A** — simplest, self-contained. Give it any gene list + count matrix,
  it runs PyDESeq2 fresh, filters, and outputs candidates ready for plots.
  Use this for the Geneious P450 list or any quick analysis.
- **Path B** — uses pre-computed step 1 results. Faster if step 1 already ran.
  Also extracts protein FASTA (needed if you want BLAST later).
- **Path C** — full pipeline with swissprot BLAST. Heatmap labels show protein
  names instead of just LOC IDs. Use for publication figures.

---

## Commands to Run

### PATH A: Direct (gene list → PyDESeq2 → filter → plots)

**One sbatch does everything — runs PyDESeq2 on the full count matrix, filters
to your gene list, applies DE cutoffs, outputs curated candidates.**

**Geneious P450 list (default):**

```bash
cd /projects/tholl_lab_1/daisy_analysis
sbatch 05_rnaseq-code/scripts/run_filter_count_genelist.sbatch
```

**Your CYP master list:**

```bash
GENE_LIST=07_NRdatabase/cyp450_database/cyp_master_list.csv \
OUT_NAME=cyp_candidates.tsv \
  sbatch 05_rnaseq-code/scripts/run_filter_count_genelist.sbatch
```

**Custom cutoffs:**

```bash
PADJ=0.01 LFC=3.0 sbatch 05_rnaseq-code/scripts/run_filter_count_genelist.sbatch
```

Check output:

```bash
wc -l 07_NRdatabase/sukman_database/geneious_candidates.tsv
head -5 07_NRdatabase/sukman_database/geneious_candidates.tsv | column -t -s $'\t'
```

Then generate plots:

```bash
sbatch 05_rnaseq-code/scripts/run_pydeseq2_step3_plots.sbatch DC \
  /projects/tholl_lab_1/daisy_analysis/07_NRdatabase/sukman_database/geneious_candidates.tsv
```

> **How it works:** Loads the full count matrix (~33K genes), runs PyDESeq2
> using the metadata (root vs leaf), computes baseMean/log2FoldChange/padj
> for every gene, then filters to only the genes in your list that pass
> padj < 0.05 AND |log2FC| > 2.0. Output has raw counts + DE stats.

---

### PATH B: Short (intersect with pre-computed step 1)

### Step 1+2: Intersect + extract proteins

**Option A — Your list only (default):**

```bash
cd /projects/tholl_lab_1/daisy_analysis
sbatch 05_rnaseq-code/scripts/run_cyp_express_extract.sbatch
```

**Option B — Also run the Geneious list separately:**

```bash
cd /projects/tholl_lab_1/daisy_analysis
PREV_LIST=07_NRdatabase/sukman_database/P450_list_RefSeq.txt \
  sbatch 05_rnaseq-code/scripts/run_cyp_express_extract.sbatch
```

Both lists are intersected with the **same** PyDESeq2 step 1 results, but produce
**separate output files** so you can compare them independently:

| Output | Source |
|---|---|
| `cyp_expressed_list.tsv` + `cyp_proteins.fasta` | Your HMMER/GTF list (~396 genes) |
| `geneious_expressed_list.tsv` + `geneious_proteins.fasta` | Geneious list (~50 genes) |

Wait for it to finish. Check output:

```bash
wc -l 07_NRdatabase/cyp450_database/cyp_expressed_list.tsv
wc -l 07_NRdatabase/cyp450_database/geneious_expressed_list.tsv
grep -c "^>" 07_NRdatabase/cyp450_database/cyp_proteins.fasta
grep -c "^>" 07_NRdatabase/cyp450_database/geneious_proteins.fasta
```

### SHORT PATH — Skip straight to plots (no BLAST, Path B only)

After Step 1+2 finishes, you can go directly to plots. The plots script
filters internally — by default padj < 0.05 and |log2FC| > 2.0 — so you do NOT
need step 2 or any separate filtering step.

**Your list:**

```bash
sbatch 05_rnaseq-code/scripts/run_pydeseq2_step3_plots.sbatch DC cyp_expressed
```

**Geneious list** (if you ran with `PREV_LIST`):

```bash
sbatch 05_rnaseq-code/scripts/run_pydeseq2_step3_plots.sbatch DC \
  /projects/tholl_lab_1/daisy_analysis/07_NRdatabase/cyp450_database/geneious_expressed_list.tsv
```

Custom cutoffs (either list):

```bash
PADJ_CUTOFF=0.01 LFC_CUTOFF=1.5 \
  sbatch 05_rnaseq-code/scripts/run_pydeseq2_step3_plots.sbatch DC cyp_expressed
```

Output: `06_analysis/pydeseq2_DC_step3_plots_cyp_expressed_*/` (yours) and
`06_analysis/pydeseq2_DC_step3_plots_geneious_expressed_list_*/` (Geneious)

> **That's it for the short path.** Each list gets its own heatmap, volcano,
> and MA plot so you can compare them side by side.
>
> If you want swissprot protein names on your plots, continue with the full path below.

---

### PATH C: Full — Steps 3-5 (with BLAST annotations)

### Step 3: BLAST CYP proteins against swissprot

Uses the existing BLAST sbatch scripts with a custom query FASTA (3rd argument).
Choose **discovery** (permissive, more hits) or **strict** (high-confidence only):

**Discovery mode** (recommended first pass):

```bash
sbatch 05_rnaseq-code/scripts/blastp_discoveryfilter.sbatch DC swissprot \
    /projects/tholl_lab_1/daisy_analysis/07_NRdatabase/cyp450_database/cyp_proteins.fasta
```

**Strict mode** (stringent e-value, high coverage):

```bash
sbatch 05_rnaseq-code/scripts/blastp_strictfilter.sbatch DC swissprot \
    /projects/tholl_lab_1/daisy_analysis/07_NRdatabase/cyp450_database/cyp_proteins.fasta
```

Both are fast (~5 min for ~396 proteins vs hours for all 33K).

Output goes to: `06_analysis/blastp_DC_cyp_proteins/`

Check when done:

```bash
ls -la 06_analysis/blastp_DC_cyp_proteins/
wc -l 06_analysis/blastp_DC_cyp_proteins/blastp_DC_swissprot_discovery.tsv
wc -l 06_analysis/blastp_DC_cyp_proteins/blastp_DC_swissprot_strict.tsv
```

### Step 4: Combine BLAST + expression, then filter

Uses the existing `run_combine_filter.sbatch` with `cyp` as the 5th argument.
The 3rd argument (`discovery` or `strict`) must match the BLAST mode you ran in step 3:

**If you ran discovery BLAST:**

```bash
sbatch 05_rnaseq-code/scripts/run_combine_filter.sbatch DC swissprot discovery standard cyp
```

**If you ran strict BLAST:**

```bash
sbatch 05_rnaseq-code/scripts/run_combine_filter.sbatch DC swissprot strict standard cyp
```

To change the DE filter cutoff (e.g. lenient):

```bash
sbatch 05_rnaseq-code/scripts/run_combine_filter.sbatch DC swissprot discovery lenient cyp
```

Output (in `07_NRdatabase/cyp450_database/`):
- `cyp_discovery_annotated.tsv` (or `cyp_strict_annotated.tsv`) — all CYP genes with expression + BLAST
- `cyp_discovery_filtered_standard.tsv` (or `cyp_strict_filtered_standard.tsv`) — only significant DE CYP genes

Check when done:

```bash
wc -l 07_NRdatabase/cyp450_database/cyp_discovery_annotated.tsv
wc -l 07_NRdatabase/cyp450_database/cyp_discovery_filtered_standard.tsv
head -5 07_NRdatabase/cyp450_database/cyp_discovery_annotated.tsv | column -t -s $'\t'
```

### Step 5: Generate plots (heatmap, volcano, MA)

```bash
# Discovery BLAST results:
sbatch 05_rnaseq-code/scripts/run_pydeseq2_step3_plots.sbatch DC cyp_discovery

# Strict BLAST results:
sbatch 05_rnaseq-code/scripts/run_pydeseq2_step3_plots.sbatch DC cyp_strict
```

Output: `06_analysis/pydeseq2_DC_step3_plots_cyp_discovery_*/` or `06_analysis/pydeseq2_DC_step3_plots_cyp_strict_*/`

### Optional: Re-filter with different cutoffs (pydeseq2_step2)

Step 4 already filters using `standard` mode (padj < 0.05, |log2FC| > 2.0).
If you want **different cutoffs** without re-running the combine step, use
`run_pydeseq2_step2_filter.sbatch`. This reads the `_annotated.tsv` (the
combined but unfiltered file from step 4) and applies new cutoffs.

**How it differs from Step 4 filtering:**
- Step 4 combines + filters in one job. You pick the filter mode (`standard`,
  `strict`, `lenient`) when you run it.
- `pydeseq2_step2_filter` re-filters the already-combined file. Useful when
  you want to try different cutoffs without re-running BLAST+combine.

```bash
# Default cutoffs (padj<0.05, |log2FC|>2.0):
sbatch 05_rnaseq-code/scripts/run_pydeseq2_step2_filter.sbatch DC cyp_discovery

# Custom cutoffs:
PADJ=0.01 LFC=3.0 sbatch 05_rnaseq-code/scripts/run_pydeseq2_step2_filter.sbatch DC cyp_discovery

# Strict BLAST version:
sbatch 05_rnaseq-code/scripts/run_pydeseq2_step2_filter.sbatch DC cyp_strict
```

Output: `06_analysis/pydeseq2_DC_step2_filtered/cyp_discovery_FILTERED.tsv`

---

## What Each Step Does

### Path A — Direct filter (`run_filter_count_genelist.sbatch`)

Runs PyDESeq2 from scratch on the **full** count matrix (~33K genes) using the
sample metadata (root vs leaf). This computes baseMean, log2FoldChange, pvalue,
and padj for every gene. Then it filters to only the genes in your gene list
and applies DE cutoffs (default padj < 0.05, |log2FC| > 2.0).

The output TSV has raw counts per sample *plus* DE stats for each candidate
gene. Feed it directly to `run_pydeseq2_step3_plots.sbatch` for heatmaps.

**Why run PyDESeq2 on the full matrix?** Size factor estimation and dispersion
estimates are more accurate with all ~33K genes. The script uses the full
matrix for statistics, then subsets the results to your gene list.

**Output:** `geneious_candidates.tsv` (or whatever `OUT_NAME` you set)

### Path B Steps 1+2 — Intersect + extract (`run_cyp_express_extract.sbatch`)

**Step 1** takes the HMMER-confirmed CYP gene_ids from `cyp_master_list.csv`
and looks for them in the PyDESeq2 unfiltered results. Keeps only CYP genes
that are in your count matrix (= expressed in your samples). Attaches expression
stats (baseMean, log2FoldChange, padj) to each gene.

If `PREV_LIST` is set, Step 1 runs a **second time** for the Geneious list
(`P450_list_RefSeq.txt`) against the same PyDESeq2 results. Each list gets
its own separate output file — they are NOT merged together.

**Step 2** reads each expressed list, looks up each gene's protein_id (XP_*),
and extracts those protein sequences from the full carrot protein FASTA.
Writes a small subset FASTA with gene_id (LOC*) as headers.

**Output (your list):**
- `cyp_expressed_list.tsv` — your CYP genes with expression stats + protein_id
- `cyp_proteins.fasta` — your protein sequences for BLAST

**Output (Geneious list, only if `PREV_LIST` set):**
- `geneious_expressed_list.tsv` — Geneious P450 genes with expression stats
- `geneious_proteins.fasta` — Geneious protein sequences

### Step 3 — BLAST (`blastp_discoveryfilter.sbatch` or `blastp_strictfilter.sbatch`)

Same BLAST scripts you already use for all proteins, but now with a 3rd argument
pointing to the CYP subset FASTA. Only ~396 proteins instead of 33K, so it runs in minutes.

- **Discovery mode:** e-value 1e-4, 40% coverage, 25 hits per query (permissive)
- **Strict mode:** e-value 1e-8, 70% coverage, 5 hits per query (high-confidence)

**Output:** `blastp_DC_swissprot_discovery.tsv` and/or `blastp_DC_swissprot_strict.tsv`
in `06_analysis/blastp_DC_cyp_proteins/`

### Step 4 — Combine + filter (`run_combine_filter.sbatch` with `cyp`)

Same sbatch you already use, but with `cyp` as the 5th argument. This tells
it to use the CYP expressed list as the DESeq input and the CYP BLAST results.
It does two things in one job:
1. Merges the best BLAST hit per gene onto the expressed list
2. Filters for significant DE genes (padj + log2FC cutoffs)

**Output:**
- `cyp_discovery_annotated.tsv` (or `cyp_strict_annotated.tsv`) — all CYP genes with expression + BLAST
- `cyp_discovery_filtered_standard.tsv` (or `cyp_strict_filtered_standard.tsv`) — only significant DE CYP genes

### Step 5 — Plots (`run_pydeseq2_step3_plots.sbatch`)

Same plots script you already use. Generates heatmap, volcano plot, MA plot, PCA,
and sample correlation heatmap. Accepts `cyp_discovery` or `cyp_strict` as a named
input source — reads the `_annotated.tsv` file produced by step 4.

**Output:** PDFs and PNGs

### Optional — Re-filter (`run_pydeseq2_step2_filter.sbatch`)

Only needed if you want different padj/LFC cutoffs than what step 4 used.
Reads the `_annotated.tsv` from step 4 and applies new thresholds. This avoids
re-running the combine + BLAST steps — just a quick filter on the existing file.

**When to use:**
- Step 4 gave you 80 DE genes with `standard` but you want fewer → use `PADJ=0.01 LFC=3.0`
- Step 4 gave you 5 DE genes with `standard` but you want more → use `PADJ=0.1 LFC=1.0`
- You already ran step 4 and don't want to wait for combine again

**Output:** `cyp_discovery_FILTERED.tsv` in `06_analysis/pydeseq2_DC_step2_filtered/`

---

## Output Files Summary

### Path A outputs

| File | Location | Description |
|------|----------|-------------|
| `geneious_candidates.tsv` | `07_NRdatabase/sukman_database/` | Geneious P450 genes — raw counts + DE stats, filtered |
| `cyp_candidates.tsv` | `07_NRdatabase/cyp450_database/` | CYP master list genes (if run with CYP list) |

### Path B outputs

| File | Location | Description |
|------|----------|-------------|
| `cyp_expressed_list.tsv` | `07_NRdatabase/cyp450_database/` | Your CYP genes in count matrix with expression stats |
| `geneious_expressed_list.tsv` | `07_NRdatabase/cyp450_database/` | Geneious P450 genes in count matrix (if PREV_LIST set) |
| `cyp_proteins.fasta` | `07_NRdatabase/cyp450_database/` | Protein sequences for BLAST (gene_id headers) |

### Path C outputs (additional)

| File | Location | Description |
|------|----------|-------------|
| `blastp_DC_swissprot_discovery.tsv` | `06_analysis/blastp_DC_cyp_proteins/` | Raw BLAST hits (discovery) |
| `blastp_DC_swissprot_strict.tsv` | `06_analysis/blastp_DC_cyp_proteins/` | Raw BLAST hits (strict) |
| `cyp_discovery_annotated.tsv` | `07_NRdatabase/cyp450_database/` | Expression + discovery BLAST combined |
| `cyp_strict_annotated.tsv` | `07_NRdatabase/cyp450_database/` | Expression + strict BLAST combined |
| `cyp_discovery_filtered_standard.tsv` | `07_NRdatabase/cyp450_database/` | Significant DE CYPs (discovery) |
| `cyp_strict_filtered_standard.tsv` | `07_NRdatabase/cyp450_database/` | Significant DE CYPs (strict) |

### All paths → plots

| File | Location | Description |
|------|----------|-------------|
| heatmap, volcano, MA, PCA | `06_analysis/pydeseq2_DC_step3_plots_*/` | Final plots (PDF + PNG) |

---

## Sanity Checks

After Path A:

```bash
wc -l 07_NRdatabase/sukman_database/geneious_candidates.tsv
head -5 07_NRdatabase/sukman_database/geneious_candidates.tsv | column -t -s $'\t'
```

After Path B steps 1+2:

```bash
wc -l 07_NRdatabase/cyp450_database/cyp_expressed_list.tsv
grep -c "^>" 07_NRdatabase/cyp450_database/cyp_proteins.fasta

# If you ran with PREV_LIST:
wc -l 07_NRdatabase/cyp450_database/geneious_expressed_list.tsv
grep -c "^>" 07_NRdatabase/cyp450_database/geneious_proteins.fasta
```

After step 3:

```bash
cut -f1 06_analysis/blastp_DC_cyp_proteins/blastp_DC_swissprot_discovery.tsv | sort -u | wc -l
```

After step 4 (discovery example):

```bash
wc -l 07_NRdatabase/cyp450_database/cyp_discovery_annotated.tsv
wc -l 07_NRdatabase/cyp450_database/cyp_discovery_filtered_standard.tsv
```

**Good result:** Most CYP genes in expressed list, most get BLAST hits,
some subset (20-80) pass the DE filter. The Geneious list is smaller (~50 genes)
so expect fewer DE genes from that set.

**Bad result:**
- 0 genes in intersection → gene_id format mismatch (check LOC prefix)
- 0 BLAST hits → check `$BLASTDB` path
- 0 genes pass filter → try lenient mode

---

## Scripts Used

| Script | Path | Role |
|--------|------|------|
| `run_filter_count_genelist.sbatch` | Path A | Gene list → PyDESeq2 → filter → candidates |
| `filter_count_by_genelist.py` | Path A | Python script called by above batch |
| `run_cyp_express_extract.sbatch` | Path B/C | Steps 1+2: intersect + extract proteins |
| `cyp_intersect_pydeseq2.py` | Path B/C | Called by step 1 — supports `--prev-list` + `--database` |
| `cyp_extract_proteins.py` | Path B/C | Called by step 2 — extracts protein FASTA |
| `blastp_discoveryfilter.sbatch` | Path C | Step 3: discovery BLAST with custom query |
| `blastp_strictfilter.sbatch` | Path C | Step 3: strict BLAST with custom query |
| `run_combine_filter.sbatch` | Path C | Step 4: combine BLAST+expression, filter DE genes |
| `run_pydeseq2_step3_plots.sbatch` | All paths | Step 5: plots (heatmap, volcano, MA, PCA) |
| `run_pydeseq2_step2_filter.sbatch` | Path B/C | Optional: re-filter with different cutoffs |
