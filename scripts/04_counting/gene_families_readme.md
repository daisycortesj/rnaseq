# Gene Families Pipeline (CYP + OMT)

Three paths to go from a gene list to heatmaps/plots. Works for **both**
CYP (cytochrome P450) and OMT (O-methyltransferase) gene families — just
change the gene family argument.

All commands run on the ARC HPC from `/projects/tholl_lab_1/daisy_analysis`.

---

## Prerequisites (already done, do not rerun)

| File | How it was made | Location |
|------|----------------|----------|
| `cyp_master_list.csv` | `sbatch scripts/08_gene_families/run_gene_database.sbatch cyp` | `07_NRdatabase/cyp450_database/` |
| `omt_master_list.csv` | `sbatch scripts/08_gene_families/run_gene_database.sbatch omt` | `07_NRdatabase/omt_database/` |
| `pydeseq2_results_UNFILTERED.tsv` | `sbatch scripts/05_pydeseq2/run_step1_analysis.sbatch DC` | `06_analysis/pydeseq2_DC_step1_unfiltered/` |
| `GCF_001625215.2_DH1_v3.0_protein.faa` | NCBI download | `04_reference/` |
| `gene_count_matrix.tsv` | `sbatch scripts/build_count_matrix.sbatch` | `03_count_tables/00_1_DC/` |
| `sample_metadata.tsv` | Created manually | `03_count_tables/00_1_DC/` |
| swissprot BLAST database | `sbatch scripts/download_blastdb.sbatch` | `07_NRdatabase/blastdb/` |
| `P450_list_RefSeq.txt` | Previous student (Geneious) — CYP gene IDs | `07_NRdatabase/sukman_database/` |
| `Methyltransferase_list.txt` | Previous student — OMT gene IDs | `07_NRdatabase/sukman_database/` |
| `P450_expression_refseq.txt` | Previous student — all expressed P450s | `07_NRdatabase/sukman_database/` |
| `MT_expression_refseq.txt` | Previous student — all expressed MTs | `07_NRdatabase/sukman_database/` |
| `P450_expression_refseq_logfold2.txt` | Previous student — DE filtered | `07_NRdatabase/sukman_database/` |
| `P450_list_RefSeq_log2fold.txt` | Previous student — DE gene IDs | `07_NRdatabase/sukman_database/` |
| `P450_upregulated.txt` | Previous student — upregulated gene IDs | `07_NRdatabase/sukman_database/` |

---

## Pipeline Overview — Three Paths

```
═══ PATH A: DIRECT (gene list → PyDESeq2 → filter → plots) ════════════════

  run_filter_genelist.sbatch DC [CYP|OMT]
    → Runs PyDESeq2 on full count matrix, filters to your gene list,
      applies DE cutoffs (padj < 0.05, |log2FC| > 2.0)
    → CYP: geneious_candidates_DC.tsv
    → OMT: omt_candidates_DC.tsv

  run_step3_plots.sbatch DC <output.tsv>
    → heatmap, volcano, MA plots

  SCRIPTS: filter_count_by_genelist.py
  BEST FOR: Geneious gene list, quick exploration, any gene list

═══ PATH B: SHORT (HMMER list → intersect → plots) ════════════════════════

  run_gene_extract.sbatch [cyp|omt]
    → CYP: cyp_expressed_list.tsv + cyp_proteins.fasta
    → OMT: omt_expressed_list.tsv + omt_proteins.fasta

  run_step3_plots.sbatch DC [cyp_expressed|omt_expressed]
    → heatmap, volcano, MA plots
    (plots script filters internally: padj < 0.05, |log2FC| > 2.0)

  SCRIPTS: gene_intersect_pydeseq2.py, gene_extract_proteins.py
  BEST FOR: Your HMMER/GTF master list, protein extraction for BLAST

═══ PATH C: FULL (HMMER list → BLAST → combine → plots) ═══════════════════

  run_gene_extract.sbatch [cyp|omt]           (same as Path B step 1)
  run_blastp_discovery.sbatch DC swissprot <proteins.fasta>
  run_combine_blast_deseq.sbatch DC swissprot discovery standard [cyp|omt]
  run_step3_plots.sbatch DC [cyp_discovery|omt_discovery]
    → heatmap, volcano, MA plots (with swissprot annotations on labels)

  BEST FOR: Publication-quality figures with protein names
```

**Which path to choose?**
- **Path A** — simplest, self-contained. Give it any gene list + count matrix,
  it runs PyDESeq2 fresh, filters, and outputs candidates ready for plots.
  Use this for the Geneious P450/OMT list or any quick analysis.
- **Path B** — uses pre-computed step 1 results. Faster if step 1 already ran.
  Also extracts protein FASTA (needed if you want BLAST later).
- **Path C** — full pipeline with swissprot BLAST. Heatmap labels show protein
  names instead of just LOC IDs. Use for publication figures.

---

## Commands to Run

### PATH A: Direct (gene list → PyDESeq2 → filter → plots)

**One sbatch does everything — runs PyDESeq2 on the full count matrix, filters
to your gene list, applies DE cutoffs, outputs curated candidates.**

#### CYP (default)

```bash
cd /projects/tholl_lab_1/daisy_analysis

# Step 1: PyDESeq2 unfiltered (same for CYP and OMT)
sbatch scripts/05_pydeseq2/run_step1_analysis.sbatch DC

# Step 2: Filter by CYP gene list (P450_list_RefSeq.txt)
sbatch scripts/08_gene_families/run_filter_genelist.sbatch DC
#   → 07_NRdatabase/sukman_database/geneious_candidates_DC.tsv

# Step 3: Plots
sbatch scripts/05_pydeseq2/run_step3_plots.sbatch DC \
  /projects/tholl_lab_1/daisy_analysis/07_NRdatabase/sukman_database/geneious_candidates_DC.tsv
```

#### OMT

```bash
cd /projects/tholl_lab_1/daisy_analysis

# Step 1: PyDESeq2 unfiltered (same for CYP and OMT — skip if already done)
sbatch scripts/05_pydeseq2/run_step1_analysis.sbatch DC

# Step 2: Filter by OMT gene list (Methyltransferase_list.txt)
sbatch scripts/08_gene_families/run_filter_genelist.sbatch DC OMT
#   → 07_NRdatabase/sukman_database/omt_candidates_DC.tsv

# Step 3: Plots
sbatch scripts/05_pydeseq2/run_step3_plots.sbatch DC \
  /projects/tholl_lab_1/daisy_analysis/07_NRdatabase/sukman_database/omt_candidates_DC.tsv
```

Check output:

```bash
# CYP
wc -l 07_NRdatabase/sukman_database/geneious_candidates_DC.tsv
head -5 07_NRdatabase/sukman_database/geneious_candidates_DC.tsv | column -t -s $'\t'

# OMT
wc -l 07_NRdatabase/sukman_database/omt_candidates_DC.tsv
head -5 07_NRdatabase/sukman_database/omt_candidates_DC.tsv | column -t -s $'\t'
```

---

### PATH B: Short (HMMER list → intersect → plots)

#### CYP

```bash
cd /projects/tholl_lab_1/daisy_analysis

# Step 1+2: Intersect CYP master list with PyDESeq2 + extract proteins
sbatch scripts/08_gene_families/run_gene_extract.sbatch cyp
#   → 07_NRdatabase/cyp450_database/cyp_expressed_list.tsv
#   → 07_NRdatabase/cyp450_database/cyp_proteins.fasta

# Step 3: Plots (heatmap, volcano, MA — filters internally: padj < 0.05, |log2FC| > 2.0)
sbatch scripts/05_pydeseq2/run_step3_plots.sbatch DC cyp_expressed
#   → 06_analysis/pydeseq2_DC_step3_plots_cyp_expressed/
```

#### OMT

```bash
cd /projects/tholl_lab_1/daisy_analysis

# Step 1+2: Intersect OMT master list with PyDESeq2 + extract proteins
sbatch scripts/08_gene_families/run_gene_extract.sbatch omt
#   → 07_NRdatabase/omt_database/omt_expressed_list.tsv
#   → 07_NRdatabase/omt_database/omt_proteins.fasta

# Step 3: Plots (heatmap, volcano, MA — filters internally: padj < 0.05, |log2FC| > 2.0)
sbatch scripts/05_pydeseq2/run_step3_plots.sbatch DC omt_expressed
#   → 06_analysis/pydeseq2_DC_step3_plots_omt_expressed/
```

Check output:

```bash
# CYP
wc -l 07_NRdatabase/cyp450_database/cyp_expressed_list.tsv
grep -c "^>" 07_NRdatabase/cyp450_database/cyp_proteins.fasta

# OMT
wc -l 07_NRdatabase/omt_database/omt_expressed_list.tsv
grep -c "^>" 07_NRdatabase/omt_database/omt_proteins.fasta
```

Custom cutoffs (either family):

```bash
PADJ_CUTOFF=0.01 LFC_CUTOFF=1.5 \
  sbatch scripts/05_pydeseq2/run_step3_plots.sbatch DC cyp_expressed

PADJ_CUTOFF=0.01 LFC_CUTOFF=1.5 \
  sbatch scripts/05_pydeseq2/run_step3_plots.sbatch DC omt_expressed
```

> **That's it for the short path.** Each family gets its own heatmap, volcano,
> and MA plot. If you want swissprot protein names on your plots, continue
> with Path C below.

---

### PATH C: Full (HMMER list → BLAST → combine → plots)

#### CYP

```bash
cd /projects/tholl_lab_1/daisy_analysis

# Step 1+2: Intersect + extract (same as Path B)
sbatch scripts/08_gene_families/run_gene_extract.sbatch cyp

# Step 3: BLASTp against swissprot
sbatch scripts/06_blast/run_blastp_discovery.sbatch DC swissprot \
    /projects/tholl_lab_1/daisy_analysis/07_NRdatabase/cyp450_database/cyp_proteins.fasta

# Step 4: Combine BLAST + DESeq2 + filter
sbatch scripts/06_blast/run_combine_blast_deseq.sbatch DC swissprot discovery standard cyp

# Step 5: Plots (with swissprot protein names on labels)
sbatch scripts/05_pydeseq2/run_step3_plots.sbatch DC cyp_discovery
#   → 06_analysis/pydeseq2_DC_step3_plots_cyp_discovery/
```

#### OMT

```bash
cd /projects/tholl_lab_1/daisy_analysis

# Step 1+2: Intersect + extract (same as Path B)
sbatch scripts/08_gene_families/run_gene_extract.sbatch omt

# Step 3: BLASTp against swissprot
sbatch scripts/06_blast/run_blastp_discovery.sbatch DC swissprot \
    /projects/tholl_lab_1/daisy_analysis/07_NRdatabase/omt_database/omt_proteins.fasta

# Step 4: Combine BLAST + DESeq2 + filter
sbatch scripts/06_blast/run_combine_blast_deseq.sbatch DC swissprot discovery standard omt

# Step 5: Plots (with swissprot protein names on labels)
sbatch scripts/05_pydeseq2/run_step3_plots.sbatch DC omt_discovery
#   → 06_analysis/pydeseq2_DC_step3_plots_omt_discovery/
```

Check output:

```bash
# CYP
wc -l 07_NRdatabase/cyp450_database/cyp_discovery_annotated.tsv
head -5 07_NRdatabase/cyp450_database/cyp_discovery_annotated.tsv | column -t -s $'\t'

# OMT
wc -l 07_NRdatabase/omt_database/omt_discovery_annotated.tsv
head -5 07_NRdatabase/omt_database/omt_discovery_annotated.tsv | column -t -s $'\t'
```

**Strict mode** (instead of discovery, use `strict` in steps 3-5):

```bash
# CYP strict
sbatch scripts/06_blast/run_blastp_discovery.sbatch DC swissprot \
    /projects/tholl_lab_1/daisy_analysis/07_NRdatabase/cyp450_database/cyp_proteins.fasta
sbatch scripts/06_blast/run_combine_blast_deseq.sbatch DC swissprot strict standard cyp
sbatch scripts/05_pydeseq2/run_step3_plots.sbatch DC cyp_strict

# OMT strict
sbatch scripts/06_blast/run_blastp_discovery.sbatch DC swissprot \
    /projects/tholl_lab_1/daisy_analysis/07_NRdatabase/omt_database/omt_proteins.fasta
sbatch scripts/06_blast/run_combine_blast_deseq.sbatch DC swissprot strict standard omt
sbatch scripts/05_pydeseq2/run_step3_plots.sbatch DC omt_strict
```

---

## Verification (any path)

After any path, verify your results against the previous student's data.
The verify script runs 6 checks and produces a comparison table + summary.

### CYP verification

```bash
cd /projects/tholl_lab_1/daisy_analysis

# Against Sukman database (default)
RESULTS=06_analysis/pydeseq2_DC_step3_plots_geneious_candidates_DC_padj005_lfc20/cyp_gene_list.tsv \
RESULTS_UNFILTERED=06_analysis/pydeseq2_DC_step1_unfiltered/pydeseq2_results_UNFILTERED.tsv \
OUTPUT=07_NRdatabase/sukman_database/verify_DC_cyp_filtered_sukman.tsv \
OUTPUT_UNFILTERED=07_NRdatabase/sukman_database/verify_DC_cyp_unfiltered_sukman.tsv \
  sbatch scripts/11_verify/run_verify_genelist.sbatch

# Against Ahmed database
DATABASE=ahmed \
RESULTS=06_analysis/pydeseq2_DC_step3_plots_geneious_candidates_DC_padj005_lfc20/cyp_gene_list.tsv \
RESULTS_UNFILTERED=06_analysis/pydeseq2_DC_step1_unfiltered/pydeseq2_results_UNFILTERED.tsv \
OUTPUT=07_NRdatabase/ahmed_database/verify_DC_cyp_filtered_ahmed.tsv \
OUTPUT_UNFILTERED=07_NRdatabase/ahmed_database/verify_DC_cyp_unfiltered_ahmed.tsv \
  sbatch scripts/11_verify/run_verify_genelist.sbatch
```

### OMT verification

```bash
cd /projects/tholl_lab_1/daisy_analysis

# Against Sukman database
GENE_FAMILY=OMT \
RESULTS=06_analysis/pydeseq2_DC_step3_plots_omt_candidates_DC_padj005_lfc20/omt_gene_list.tsv \
RESULTS_UNFILTERED=06_analysis/pydeseq2_DC_step1_unfiltered/pydeseq2_results_UNFILTERED.tsv \
OUTPUT=07_NRdatabase/sukman_database/verify_DC_omt_filtered_sukman.tsv \
OUTPUT_UNFILTERED=07_NRdatabase/sukman_database/verify_DC_omt_unfiltered_sukman.tsv \
  sbatch scripts/11_verify/run_verify_genelist.sbatch
```

> **Note:** The OMT verification currently only has `Methyltransferase_list.txt`
> (Check 1) and `MT_expression_refseq.txt` (Check 3). Checks 2, 4, 5, 6 are
> skipped because the Sukman database doesn't have OMT-specific log2fold,
> upregulated, or downregulated files yet. When those files are added, the
> checks will run automatically.

### Download results

```bash
# CYP
scp daisycortesj@tinkercliffs.arc.vt.edu:/projects/tholl_lab_1/daisy_analysis/07_NRdatabase/sukman_database/verify_DC_cyp_*.tsv ~/Desktop/

# OMT
scp daisycortesj@tinkercliffs.arc.vt.edu:/projects/tholl_lab_1/daisy_analysis/07_NRdatabase/sukman_database/verify_DC_omt_*.tsv ~/Desktop/
```

---

## What Each Step Does

### Path A — Direct filter (`run_filter_genelist.sbatch`)

Runs PyDESeq2 from scratch on the **full** count matrix (~33K genes) using the
sample metadata (root vs leaf). Computes baseMean, log2FoldChange, pvalue, and
padj for every gene. Then filters to only genes in your gene list and applies
DE cutoffs (default padj < 0.05, |log2FC| > 2.0).

The 2nd argument picks the gene family:
- `DC` or `DC CYP` → uses `P450_list_RefSeq.txt` → `geneious_candidates_DC.tsv`
- `DC OMT` → uses `Methyltransferase_list.txt` → `omt_candidates_DC.tsv`

**Why run PyDESeq2 on the full matrix?** Size factor estimation and dispersion
estimates are more accurate with all ~33K genes. The script uses the full
matrix for statistics, then subsets the results to your gene list.

### Path B Steps 1+2 — Intersect + extract (`run_gene_extract.sbatch`)

**Step 1** takes the HMMER-confirmed gene_ids from the master list
(`cyp_master_list.csv` or `omt_master_list.csv`) and looks for them in the
PyDESeq2 unfiltered results. Keeps only genes that are in your count matrix
(= expressed in your samples). Attaches expression stats (baseMean,
log2FoldChange, padj) to each gene.

**Step 2** reads the expressed list, looks up each gene's protein_id (XP_*),
and extracts those protein sequences from the full carrot protein FASTA.
Writes a small subset FASTA with gene_id (LOC*) as headers.

**CYP output:** `cyp_expressed_list.tsv` + `cyp_proteins.fasta` in `07_NRdatabase/cyp450_database/`
**OMT output:** `omt_expressed_list.tsv` + `omt_proteins.fasta` in `07_NRdatabase/omt_database/`

### Step 3 — BLAST (`run_blastp_discovery.sbatch`)

Same BLAST scripts you already use for all proteins, but with a 3rd argument
pointing to the CYP/OMT subset FASTA. Only a few hundred proteins instead of
33K, so it runs in minutes.

- **Discovery mode:** e-value 1e-4, 40% coverage, 25 hits per query (permissive)
- **Strict mode:** e-value 1e-8, 70% coverage, 5 hits per query (high-confidence)

### Step 4 — Combine + filter (`run_combine_blast_deseq.sbatch`)

The 5th argument (`cyp` or `omt`) tells it which expressed list and BLAST
results to use. It merges the best BLAST hit per gene onto the expressed list,
then filters for significant DE genes (padj + log2FC cutoffs).

### Step 5 — Plots (`run_step3_plots.sbatch`)

Same plots script for both families. Generates heatmap, volcano plot, MA plot,
PCA, and sample correlation heatmap. Accepts these named input sources:

| Input source | Gene family | Description |
|---|---|---|
| `cyp_expressed` | CYP | No BLAST, direct from intersect |
| `cyp_discovery` | CYP | With swissprot annotations (discovery BLAST) |
| `cyp_strict` | CYP | With swissprot annotations (strict BLAST) |
| `omt_expressed` | OMT | No BLAST, direct from intersect |
| `omt_discovery` | OMT | With swissprot annotations (discovery BLAST) |
| `omt_strict` | OMT | With swissprot annotations (strict BLAST) |

---

## Output Files Summary

### Path A outputs

| File | Location | Description |
|------|----------|-------------|
| `geneious_candidates_DC.tsv` | `07_NRdatabase/sukman_database/` | CYP genes — raw counts + DE stats, filtered |
| `omt_candidates_DC.tsv` | `07_NRdatabase/sukman_database/` | OMT genes — raw counts + DE stats, filtered |

### Path B outputs

| File | Location | Description |
|------|----------|-------------|
| `cyp_expressed_list.tsv` | `07_NRdatabase/cyp450_database/` | CYP genes in count matrix with expression stats |
| `cyp_proteins.fasta` | `07_NRdatabase/cyp450_database/` | CYP protein sequences for BLAST |
| `omt_expressed_list.tsv` | `07_NRdatabase/omt_database/` | OMT genes in count matrix with expression stats |
| `omt_proteins.fasta` | `07_NRdatabase/omt_database/` | OMT protein sequences for BLAST |

### Path C outputs (additional)

| File | Location | Description |
|------|----------|-------------|
| `cyp_discovery_annotated.tsv` | `07_NRdatabase/cyp450_database/` | CYP expression + discovery BLAST combined |
| `cyp_discovery_filtered_standard.tsv` | `07_NRdatabase/cyp450_database/` | Significant DE CYPs (discovery) |
| `omt_discovery_annotated.tsv` | `07_NRdatabase/omt_database/` | OMT expression + discovery BLAST combined |
| `omt_discovery_filtered_standard.tsv` | `07_NRdatabase/omt_database/` | Significant DE OMTs (discovery) |

### Verification outputs

| File | Location | Description |
|------|----------|-------------|
| `verify_DC_cyp_filtered_sukman.tsv` | `07_NRdatabase/sukman_database/` | CYP filtered results vs Sukman |
| `verify_DC_cyp_filtered_SUMMARY_sukman.txt` | `07_NRdatabase/sukman_database/` | CYP human-readable summary |
| `verify_DC_omt_filtered_sukman.tsv` | `07_NRdatabase/sukman_database/` | OMT filtered results vs Sukman |
| `verify_DC_omt_filtered_SUMMARY_sukman.txt` | `07_NRdatabase/sukman_database/` | OMT human-readable summary |

---

## Previous Student Reference Files

All in `07_NRdatabase/sukman_database/`:

### CYP (P450) files

| File | Contents |
|------|----------|
| `P450_list_RefSeq.txt` | Full P450 gene list (gene IDs, one per line) |
| `P450_list_RefSeq_log2fold.txt` | Subset: DE genes only (gene IDs) |
| `P450_expression_refseq.txt` | All expressed P450 genes with baseMean, log2FC, padj |
| `P450_expression_refseq_logfold2.txt` | DE-filtered (|log2FC| > 2) with full DESeq2 stats |
| `P450_upregulated.txt` | Upregulated genes (IDs without LOC prefix) |
| `P450_downregulated_list.txt` | Downregulated genes |

### OMT (Methyltransferase) files

| File | Contents |
|------|----------|
| `Methyltransferase_list.txt` | Full OMT gene list (gene IDs, one per line) |
| `MT_expression_refseq.txt` | All expressed OMT genes with expression data |

> **Note:** The OMT database currently only has the gene list and expression
> file. When additional files (log2fold, upregulated, downregulated) are added,
> the verification script will pick them up automatically.

---

## Contrast Direction: R vs L (DC) or L vs R (DG)

| Contrast | Positive log2FC means | Negative log2FC means | Use for |
|----------|----------------------|----------------------|---------|
| `CONTRAST_A=R CONTRAST_B=L` (default) | Up in **root** | Up in **leaf** | DC (root focus) |
| `CONTRAST_A=L CONTRAST_B=R` | Up in **leaf** | Up in **root** | DG (leaf focus) |

The default across the pipeline is `R vs L`. To switch to `L vs R`, prepend
`CONTRAST_A=L CONTRAST_B=R` before any `sbatch` command. You must use the
**same contrast** for all steps in a run (step 1, filter, step 3, verify).

---

## Scripts Used

| Script | Path | Role |
|--------|------|------|
| `run_filter_genelist.sbatch` | Path A | Gene list → PyDESeq2 → filter → candidates |
| `filter_count_by_genelist.py` | Path A | Python script called by above batch |
| `run_gene_extract.sbatch` | Path B/C | Steps 1+2: intersect + extract proteins |
| `gene_intersect_pydeseq2.py` | Path B/C | Called by step 1 — intersects master list with PyDESeq2 |
| `gene_extract_proteins.py` | Path B/C | Called by step 2 — extracts protein FASTA |
| `run_blastp_discovery.sbatch` | Path C | Step 3: discovery BLAST with custom query |
| `run_combine_blast_deseq.sbatch` | Path C | Step 4: combine BLAST+expression, filter DE genes |
| `run_step3_plots.sbatch` | All paths | Step 5: plots (heatmap, volcano, MA, PCA) |
| `run_verify_genelist.sbatch` | All paths | 6-check verification vs previous student |
| `verify_genelist.py` | All paths | Python verification script |

---

## Changelog

### 2026-04-01: Added OMT support across pipeline

**What changed:** The pipeline previously only supported CYP (P450) genes.
Now supports both CYP and OMT (O-methyltransferase) via a gene family
argument/env var.

**Files changed:**

| File | What was changed |
|------|-----------------|
| `scripts/08_gene_families/run_filter_genelist.sbatch` | Added 2nd positional arg (`CYP` or `OMT`). OMT auto-selects `Methyltransferase_list.txt` as gene list and `omt_candidates_DC.tsv` as output. |
| `scripts/11_verify/run_verify_genelist.sbatch` | Added `GENE_FAMILY` env var. OMT uses `Methyltransferase_list.txt` + `MT_expression_refseq.txt`. Output files named with `_omt_` tag. Ahmed database gracefully falls back to sukman for OMT. |
| `scripts/11_verify/verify_genelist.py` | Added `--gene-family` arg. Summary header now says "CYP VERIFICATION SUMMARY" or "OMT VERIFICATION SUMMARY" (was hardcoded "P450"). |
| `scripts/11_verify/verify_genelist_beginner.py` | Same `--gene-family` arg as above. |
| `scripts/04_counting/CYP_heatmap.md` | Replaced by this file (`gene_families_readme.md`). |

Look for `# CHANGED: 2026-04-01` comments in the modified files for exact change locations.

### 2026-03-28: Contrast direction bugfix

**Problem:** Several scripts had hardcoded `R vs L` contrast assumptions.

**Files changed:**

| File | What was fixed |
|------|---------------|
| `scripts/08_gene_families/filter_count_by_genelist.py` | Print labels were hardcoded "Upregulated (root)". Now reads `--contrast-a`/`--contrast-b`. |
| `scripts/11_verify/verify_genelist.py` | Direction checks assumed positive log2FC = up in root. Added `--contrast-a`/`--contrast-b` args. |
| `scripts/11_verify/run_verify_genelist.sbatch` | Now passes `--contrast-a`/`--contrast-b` to the Python script. |
