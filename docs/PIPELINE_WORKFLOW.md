# RNA-seq Differential Expression and Gene Family Analysis Pipeline

Complete step-by-step guide for analyzing CYP (cytochrome P450) and OMT
(O-methyltransferase) gene families across *Daucus carota* (DC) and
*Daucus glaber* (DG) using the DC reference genome.

All commands are run from the HPC with `conda activate rnaseq`.

---

## Pipeline Overview

```
Alignment ──> Counts ──> DESeq2 Step 1 ──> BLAST ──> Combine ──> Step 2 (filter) ──> Step 3 (plots)
                             │                │         │                               ↑
                             │                │         └── also directly to Step 3 ────┘
                             │                │
                             │           HMMER (domain verification)
                             │
                       build_count_matrix
```

**Two paths through Steps 2 and 3:**

- **Path A (Filtered):** Step 8 combined annotated file --> Step 2 (filter) --> Step 3 (plots)
- **Path B (Unfiltered):** Step 8 combined annotated file --> Step 3 (plots directly)

Step 3 auto-detects CYP and OMT genes from BLAST annotations and generates
separate heatmaps for each family. No standalone heatmap scripts needed.

Run each step for **both DC and DG** before moving to the combined comparison.

---

## Progress Checklist

| Step | DC | DG | Script |
|------|----|----|--------|
| 0. Build count matrix | | | `04_counting/run_count_matrix.sbatch` |
| 1. PyDESeq2 step 1 (unfiltered) | | | `05_pydeseq2/run_step1_analysis.sbatch` |
| 2. PyDESeq2 step 2 (filter) | | | `05_pydeseq2/run_step2_filter.sbatch` |
| 3. PyDESeq2 step 3 (plots) | | | `05_pydeseq2/run_step3_plots.sbatch` |
| 4. Prepare BLAST input | | | `06_blast/run_blast_input.sbatch` |
| 5. Translate CDS to protein | | | `06_blast/run_translate_cds.sbatch` |
| 6. BLASTp (discovery) | | | `06_blast/run_blastp_discovery.sbatch` |
| 7. HMMER domain scan | | | `07_domains/run_hmmer.sbatch` |
| 8. Combine BLAST + DESeq2 | | | `06_blast/run_combine_blast_deseq.sbatch` |
| 9. Step 2 on annotated data | | | `05_pydeseq2/run_step2_filter.sbatch` |
| 10. Step 3 plots from annotated | | | `05_pydeseq2/run_step3_plots.sbatch` |
| 11. Compare DC vs DG | -- | -- | `09_comparison/run_compare_species.sbatch` |
| 12. Gene families combined DC+DG | -- | -- | `08_gene_families/run_gene_families_combined.sbatch` |

---

## Step 0: Build Count Matrix

**Prerequisite:** STAR alignment output (`ReadsPerGene.out.tab` files).

```bash
sbatch scripts/04_counting/run_count_matrix.sbatch DC
sbatch scripts/04_counting/run_count_matrix.sbatch DG
```

**Output:**
- `03_count_tables/00_1_DC/gene_count_matrix.tsv`
- `03_count_tables/00_1_DC/sample_metadata.tsv`

---

## Step 1: PyDESeq2 Unfiltered Analysis

Runs differential expression on ALL genes. Produces unfiltered results and
`all_gene_ids.txt` needed for BLAST input.

```bash
sbatch scripts/05_pydeseq2/run_step1_analysis.sbatch DC
sbatch scripts/05_pydeseq2/run_step1_analysis.sbatch DG
```

**Output:**
- `06_analysis/pydeseq2_DC_step1_unfiltered/pydeseq2_results_UNFILTERED.tsv`
- `06_analysis/pydeseq2_DC_step1_unfiltered/all_gene_ids.txt`

---

## Step 2: PyDESeq2 Filter (Flexible Input)

Filters any DESeq2-like TSV by padj and log2FC. Accepts multiple input sources.

```bash
# Filter raw unfiltered results (default):
sbatch scripts/05_pydeseq2/run_step2_filter.sbatch DC

# Filter combined BLAST+DESeq2 annotated file:
sbatch scripts/05_pydeseq2/run_step2_filter.sbatch DC combined_annotated

# Filter gene-family verified file:
sbatch scripts/05_pydeseq2/run_step2_filter.sbatch DC gene_families

# Custom thresholds:
PADJ=0.01 LFC=3.0 sbatch scripts/05_pydeseq2/run_step2_filter.sbatch DC combined_annotated
```

**Input sources:**
- `unfiltered` (default) -- `pydeseq2_{SP}_step1_unfiltered/pydeseq2_results_UNFILTERED.tsv`
- `combined_annotated` -- `combined_{SP}/{SP}_swissprot_discovery_annotated.tsv`
- `gene_families` -- `gene_families_{SP}/{SP}_CYP_OMT_combined.tsv`
- or any explicit file path

**Output:**
- `06_analysis/pydeseq2_DC_step2_filtered/pydeseq2_results_FILTERED.tsv` (for unfiltered source)
- `06_analysis/pydeseq2_DC_step2_filtered/combined_annotated_FILTERED.tsv` (for combined_annotated source)
- `06_analysis/pydeseq2_DC_step2_filtered/gene_families_FILTERED.tsv` (for gene_families source)

---

## Step 3: PyDESeq2 Plots (MA, Volcano, CYP/OMT Heatmaps)

Generates all plots in one step. Auto-detects CYP and OMT genes from BLAST
annotations (or `gene_family` column) and generates separate heatmaps for each.
HMMER domain labels are added automatically if the HMMER file exists.

```bash
# From Step 2 filtered results (default):
sbatch scripts/05_pydeseq2/run_step3_plots.sbatch DC

# Directly from combined annotated file (skip Step 2):
sbatch scripts/05_pydeseq2/run_step3_plots.sbatch DC combined_annotated

# From gene-family verified file:
sbatch scripts/05_pydeseq2/run_step3_plots.sbatch DC gene_families

# From unfiltered results (all genes):
sbatch scripts/05_pydeseq2/run_step3_plots.sbatch DC unfiltered
```

**Input sources:** same as Step 2 plus `filtered` (default).

**Output (varies by input):**
- `06_analysis/pydeseq2_DC_step3_plots/ma_plot.pdf`
- `06_analysis/pydeseq2_DC_step3_plots/volcano_plot.pdf`
- `06_analysis/pydeseq2_DC_step3_plots/cyp_heatmap.pdf` (if CYP genes detected)
- `06_analysis/pydeseq2_DC_step3_plots/omt_heatmap.pdf` (if OMT genes detected)
- `06_analysis/pydeseq2_DC_step3_plots/cyp_gene_list.tsv`
- `06_analysis/pydeseq2_DC_step3_plots/omt_gene_list.tsv`

When using `combined_annotated` as input, the output goes to
`pydeseq2_DC_step3_plots_annotated/` to keep it separate from the default.

**Download plots:**
```bash
scp daisycortesj@tinkercliffs.arc.vt.edu:/projects/tholl_lab_1/daisy_analysis/06_analysis/pydeseq2_DC_step3_plots_annotated/*.pdf ~/Desktop/
```

---

## Steps 4-5: BLAST Input Preparation

```bash
sbatch scripts/06_blast/run_blast_input.sbatch DC
sbatch scripts/06_blast/run_blast_input.sbatch DG
sbatch scripts/06_blast/run_translate_cds.sbatch DC
sbatch scripts/06_blast/run_translate_cds.sbatch DG
```

**Output:** `06_analysis/blast_input_DC/all_genes_protein.fasta`

---

## Step 6: BLASTp Search

```bash
sbatch scripts/06_blast/run_blastp_discovery.sbatch DC swissprot
sbatch scripts/06_blast/run_blastp_discovery.sbatch DG swissprot
```

If the job times out, resubmit the same command -- it resumes automatically.

**Output:** `06_analysis/blastp_DC/blastp_DC_swissprot_discovery.tsv`

---

## Step 7: HMMER Domain Scan

```bash
sbatch scripts/07_domains/run_hmmer.sbatch DC all_genes_protein.fasta
sbatch scripts/07_domains/run_hmmer.sbatch DG all_genes_protein.fasta
```

**Output:** `06_analysis/hmmer_DC/all_genes_protein_pfam_domains.txt`

HMMER is optional but improves confidence scoring and adds domain labels to
Step 3 heatmaps. Steps 6 and 7 can run in parallel.

---

## Step 8: Combine BLAST + DESeq2

Merges BLAST annotations with DESeq2 expression statistics.

```bash
sbatch scripts/06_blast/run_combine_blast_deseq.sbatch DC swissprot discovery
sbatch scripts/06_blast/run_combine_blast_deseq.sbatch DG swissprot discovery
```

**Output:**
- `06_analysis/combined_DC/DC_swissprot_discovery_annotated.tsv` (all genes, annotated)

This annotated file is the main input for Steps 9-10 below.

---

## Steps 9-10: Filter and Plot the Annotated Data

This is where everything comes together. Two paths:

### Path A: Filter first, then plot

```bash
# Filter the annotated file (keeps padj < 0.05, |log2FC| > 2):
sbatch scripts/05_pydeseq2/run_step2_filter.sbatch DC combined_annotated

# Generate plots from filtered annotated data:
sbatch scripts/05_pydeseq2/run_step3_plots.sbatch DC combined_annotated
```

### Path B: Plot directly (no filter)

```bash
# Generate plots from all annotated genes (unfiltered):
sbatch scripts/05_pydeseq2/run_step3_plots.sbatch DC combined_annotated
```

Both paths auto-detect CYP/OMT genes and generate separate heatmaps.
HMMER domain labels are added if the HMMER file exists from Step 7.

---

## Step 11: Compare DC vs DG

Cross-species comparison of annotated results.

```bash
sbatch scripts/09_comparison/run_compare_species.sbatch DC DG swissprot discovery
```

**Output:**
- `06_analysis/comparison_DC_vs_DG/DC_vs_DG_swissprot_discovery_comparison.tsv`
- Subsets: `*_shared_same_direction.tsv`, `*_DC_only.tsv`, `*_DG_only.tsv`

---

## Step 12: Combined Gene Family Extraction (DC + DG)

Extracts verified CYP and OMT genes from both species with BLAST+HMMER
confidence scoring. Produces merged tables and per-species heatmaps.

```bash
# BLAST only:
sbatch scripts/08_gene_families/run_gene_families_combined.sbatch DC DG swissprot discovery

# BLAST + HMMER:
sbatch scripts/08_gene_families/run_gene_families_combined.sbatch DC DG swissprot discovery full
```

**Output:**
- `06_analysis/gene_families_DC_DG/DC_DG_CYP_verified.tsv`
- `06_analysis/gene_families_DC_DG/DC_DG_OMT_verified.tsv`
- `06_analysis/gene_families_DC_DG/DC_DG_CYP_OMT_combined.tsv`
- Category subsets: `*_shared_same_direction.tsv`, `*_DC_only.tsv`, `*_DG_only.tsv`
- Heatmaps: `cyp_heatmap_DC.pdf`, `cyp_heatmap_DG.pdf`, `omt_heatmap_DC.pdf`, `omt_heatmap_DG.pdf`

---

## Quick Reference: Full Run Order

```bash
# ── STEP 0: Count matrices ──
sbatch scripts/04_counting/run_count_matrix.sbatch DC
sbatch scripts/04_counting/run_count_matrix.sbatch DG

# ── STEP 1: DESeq2 unfiltered ──
sbatch scripts/05_pydeseq2/run_step1_analysis.sbatch DC
sbatch scripts/05_pydeseq2/run_step1_analysis.sbatch DG

# ── STEPS 2-3 (early): Quick stats-only filter + plots ──
sbatch scripts/05_pydeseq2/run_step2_filter.sbatch DC
sbatch scripts/05_pydeseq2/run_step2_filter.sbatch DG
sbatch scripts/05_pydeseq2/run_step3_plots.sbatch DC
sbatch scripts/05_pydeseq2/run_step3_plots.sbatch DG

# ── STEPS 4-5: BLAST input ──
sbatch scripts/06_blast/run_blast_input.sbatch DC
sbatch scripts/06_blast/run_blast_input.sbatch DG
sbatch scripts/06_blast/run_translate_cds.sbatch DC
sbatch scripts/06_blast/run_translate_cds.sbatch DG

# ── STEPS 6-7: BLAST + HMMER (can run in parallel) ──
sbatch scripts/06_blast/run_blastp_discovery.sbatch DC swissprot
sbatch scripts/06_blast/run_blastp_discovery.sbatch DG swissprot
sbatch scripts/07_domains/run_hmmer.sbatch DC all_genes_protein.fasta
sbatch scripts/07_domains/run_hmmer.sbatch DG all_genes_protein.fasta

# ── STEP 8: Combine BLAST + DESeq2 ──
sbatch scripts/06_blast/run_combine_blast_deseq.sbatch DC swissprot discovery
sbatch scripts/06_blast/run_combine_blast_deseq.sbatch DG swissprot discovery

# ── STEPS 9-10: Filter + plot the annotated data ──
# Path A (filtered):
sbatch scripts/05_pydeseq2/run_step2_filter.sbatch DC combined_annotated
sbatch scripts/05_pydeseq2/run_step3_plots.sbatch DC combined_annotated
# Path B (unfiltered, skip step 2):
sbatch scripts/05_pydeseq2/run_step3_plots.sbatch DC combined_annotated

# ── Same for DG ──
sbatch scripts/05_pydeseq2/run_step2_filter.sbatch DG combined_annotated
sbatch scripts/05_pydeseq2/run_step3_plots.sbatch DG combined_annotated

# ── STEP 11: Cross-species comparison ──
sbatch scripts/09_comparison/run_compare_species.sbatch DC DG swissprot discovery

# ── STEP 12: Combined gene family extraction ──
sbatch scripts/08_gene_families/run_gene_families_combined.sbatch DC DG swissprot discovery full
```

---

## Directory Structure

```
/projects/tholl_lab_1/daisy_analysis/
├── 03_count_tables/
│   ├── 00_1_DC/          gene_count_matrix.tsv, sample_metadata.tsv
│   └── 00_2_DG/          gene_count_matrix.tsv, sample_metadata.tsv
├── 05_rnaseq-code/       this repository (scripts/)
├── 06_analysis/
│   ├── pydeseq2_DC_step1_unfiltered/   UNFILTERED results + all_gene_ids.txt
│   ├── pydeseq2_DC_step2_filtered/     FILTERED results (stats-only or annotated)
│   ├── pydeseq2_DC_step3_plots/        Plots from filtered results
│   ├── pydeseq2_DC_step3_plots_annotated/  Plots from annotated data
│   ├── blast_input_DC/                 CDS + protein FASTA files
│   ├── blastp_DC/                      BLAST results
│   ├── hmmer_DC/                       HMMER Pfam domain results
│   ├── combined_DC/                    BLAST+DESeq2 merged tables
│   ├── gene_families_DC/               CYP/OMT verified lists (optional)
│   ├── comparison_DC_vs_DG/            Cross-species DE comparison
│   └── gene_families_DC_DG/            Combined CYP/OMT across species
└── 07_NRdatabase/        BLAST databases
```

---

## Troubleshooting

**"File not found" errors:** Each step depends on the previous one. Check the
progress checklist and make sure all prerequisites are complete.

**BLAST times out:** Resubmit the same `sbatch` command. Built-in resume support.

**No CYP/OMT heatmaps generated:** The input file needs either a
`blast_description` column (from Step 8 combine) or a `gene_family` column
(from gene family extraction). Raw DESeq2 results from Step 1/2 won't have
BLAST annotations, so only MA and volcano plots are generated.

**Low confidence scores in gene families:** Re-run with HMMER results available
(`full` mode). Genes confirmed by both BLAST and HMMER get `high` confidence.
