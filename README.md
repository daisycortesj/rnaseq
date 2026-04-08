# RNA-seq Differential Expression & Gene Family Pipeline

Analysis pipeline for CYP (cytochrome P450) and OMT (O-methyltransferase) gene
families in *Daucus carota* (DC) and *Daucus glaber* (DG). 

---

## Quick Start

```bash
# 1. Clone and activate environment
git clone <repo-url> && cd rnaseq
conda env create -f environment.yml
conda activate rnaseq

# 2. Install the rna_pipeline package (one-time)
pip install -e .

# 3. Verify all tools are installed
rna-pipeline --mode check-tools

# 4. Run the pipeline step by step (from scripts/)
sbatch scripts/01_qc/run_fastqc.sbatch
sbatch scripts/01_qc/run_fastp.sbatch
# ... continue through numbered stages
```

See [docs/PIPELINE_WORKFLOW.md](docs/PIPELINE_WORKFLOW.md) for the full
step-by-step runbook with commands and expected outputs.

---

## Repository Layout

```
rnaseq/
├── README.md                 ← You are here
├── environment.yml           ← Conda environment (all tools + Python deps)
├── pyproject.toml            ← Makes rna_pipeline installable via pip
│
├── scripts/                  ← USER-FACING: sbatch jobs + helper Python
│   ├── 01_qc/               │  FastQC, fastp
│   ├── 02_alignment/         │  STAR index, align, samtools sort
│   ├── 03_assembly/          │  Trinity de novo assembly
│   ├── 04_counting/          │  featureCounts, build count matrix
│   ├── 05_pydeseq2/          │  PyDESeq2 3-step workflow + R_pydeseq2 all-in-one
│   ├── 06_blast/             │  CDS extraction, BLASTp/BLASTx
│   ├── 07_domains/           │  HMMER Pfam, PROSITE motifs
│   ├── 08_gene_families/     │  Gene family extraction (CYP, OMT)
│   ├── 09_comparison/        │  Cross-species DC vs DG
│   ├── 10_post_analysis/     │  Phylogenetic trees, genomic clustering
│   ├── 11_verify/            │  Verification against prior results
│   └── setup/                │  One-time database downloads
│
├── rna_pipeline/             ← CORE LIBRARY: pipeline engine + tool wrappers
│   ├── cli.py                │  Command-line interface (rna-pipeline command)
│   ├── main.py               │  Workflow orchestration + check-tools mode
│   ├── tools/                │  Python wrappers for all pipeline tools:
│   │                         │    qc, star, trinity, featurecounts,
│   │                         │    blast, hmmer, prosite, pydeseq2
│   ├── runners/              │  Subprocess execution
│   └── utils/                │  I/O and system helpers
│
├── docs/                     ← Documentation
│   ├── PIPELINE_WORKFLOW.md  │  Master runbook (Steps 0–12)
│   ├── HMMER_PROSITE_GUIDE.md   Domain analysis setup & usage
│   ├── DIRECTORY_STRUCTURE.md   Detailed file map
│   ├── CONCEPTS.md           │  RNA-seq background for beginners
│   └── SETUP.md              │  Environment & HPC setup
│
└── archive/                  ← Deprecated scripts & old docs (kept for reference)
```

---

## Pipeline at a Glance

```
 01 QC ──→ 02 Alignment ──→ 03 Assembly ──→ 04 Counting
                                                  │
                                                  ▼
                                          05 PyDESeq2 Step 1 (statistics)
                                                  │
                    ┌──────────────────────┬──────────────────────────┐
                    │                      │                          │
               PATH A (quick)     PATH B (short family)    PATH C (full family)
                    │                      │                          │
                    ▼                      ▼                          ▼
           05 Step 2 (filter)    08 run_gene_extract        08 run_gene_extract
                    │                 (cyp or omt)               (cyp or omt)
                    ▼                      │                          │
           05 Step 3 (plots)               │                          ▼
           MA, volcano, PCA,               │                 06 BLASTp (swissprot)
           correlation only                │                          │
                                           │                          ▼
                                           │                 06 Combine + Filter
                                           │                          │
                                           ▼                          ▼
                                   05 Step 3 (plots)         05 Step 3 (plots)
                                   CYP/OMT heatmaps          CYP/OMT heatmaps
                                   (no BLAST labels)         (with protein names)
                                                                      │
                                                              07 Domains (HMMER)
                                                                      │
                                                              09 Cross-Species
                                                                      │
                                                              10 Post-Analysis
                                                                      │
                                                              11 Verification
```

**Why three paths?**
- **Path A** — quick plots (MA, volcano, PCA, correlation). No BLAST needed.
- **Path B** — short gene-family path. Uses your HMMER/GTF master list
  intersected with PyDESeq2, goes straight to heatmaps (no BLAST labels).
- **Path C** — full gene-family path. Same as Path B but adds BLASTp
  annotation so heatmap labels show protein names (publication-quality).

All paths start from Step 1. Choose your path based on whether you need
gene family heatmaps and whether you want BLAST annotations on the labels.

Each numbered directory in `scripts/` matches a pipeline stage. Each
directory contains `.sbatch` files you submit with `sbatch` and `.py` helper
scripts that the batch jobs call automatically.

---

## Verify Your Setup

After activating the environment, run this to confirm every tool is installed:

```bash
conda activate rnaseq
rna-pipeline --mode check-tools
```

This checks all 9 tools across the pipeline (STAR, Trinity, featureCounts,
PyDESeq2, BLAST, HMMER, PROSITE, FastQC, samtools) and reports which are
available and which are missing. Each sbatch script also verifies its own
tool before running, so a missing tool will fail fast with a clear message
instead of crashing halfway through a job.

---

## Key Commands

### Stages 01–04: QC through Counting (run once per species)

All scripts accept an optional species code (e.g., `MF`, `DC`) to process
only that species. Omit the code to process all species.

| Stage | Command | Notes |
|-------|---------|-------|
| QC (raw) | `sbatch scripts/01_qc/run_fastqc.sbatch MF` | Check raw read quality |
| QC (clean) | `sbatch scripts/01_qc/run_fastqc.sbatch MF clean` | Check quality after fastp |
| Trim reads | `sbatch scripts/01_qc/run_fastp.sbatch MF` | Adapter removal + quality filtering |
| STAR index | `sbatch scripts/02_alignment/run_genome_index.sbatch carrot` | Build once per genome |
| Align | `sbatch scripts/02_alignment/run_star_align_all.sbatch MF` | Auto-detects clean vs raw reads |
| Assembly | `sbatch scripts/03_assembly/run_trinity_all.sbatch MF` | Auto-detects clean vs raw reads |
| Count reads + build matrix | `sbatch scripts/04_counting/run_featurecounts.sbatch DC` | |
| Rebuild metadata only | `sbatch scripts/04_counting/run_count_matrix.sbatch DC` | |

#### Trinity assembly scripts (`scripts/03_assembly/`)

| Script | What it does | When to use |
|--------|-------------|-------------|
| `run_trinity_all.sbatch` | Assembles all samples for one species (or all species). Processes samples one after another in a single job. | **Main script.** Use for most cases. |
| `run_trinity.sbatch` | Assembles one specific sample (e.g., just MFF1). | Test a single sample first, or re-run one that failed. |
| `run_trinity_array.sbatch` | Runs as a SLURM job array — each sample gets its own independent job on a separate node. | Run samples in parallel for faster wall time (uses more cluster resources). |
| `run_trinity_rsem.sbatch` | Quantifies expression (RSEM) using completed Trinity assemblies. Not assembly — this is a later step. | **After** all assemblies are done. Maps reads back to transcripts for gene counts. |

All assembly scripts auto-detect whether fastp-cleaned reads exist in
`01_processed/`. If they do, those are used. Otherwise raw reads from
`00_rawdata/` are used automatically.

### Stage 05: PyDESeq2 Step 1 (always run this first)

| Stage | Command |
|-------|---------|
| DESeq2 step 1 (DC, root focus) | `sbatch scripts/05_pydeseq2/run_step1_analysis.sbatch DC` |
| DESeq2 step 1 (DC, single-factor override) | `DESIGN="condition" sbatch scripts/05_pydeseq2/run_step1_analysis.sbatch DC` |
| DESeq2 step 1 (DG, leaf focus) | `CONTRAST_A=L CONTRAST_B=R sbatch scripts/05_pydeseq2/run_step1_analysis.sbatch DG` |

> **Contrast direction:** Default is `R vs L` (positive log2FC = up in root).
> For DG leaf-focus analysis, prepend `CONTRAST_A=L CONTRAST_B=R` to flip the
> contrast so positive log2FC = up in leaf. Use the **same contrast** for all
> steps in a run. See [CYP_heatmap.md](scripts/04_counting/CYP_heatmap.md)
> for full DC and DG command examples.

> **Multi-factor design (auto-detected for DC):** DC automatically uses
> `DESIGN="variety+condition"` because it has two varieties (DC1 + DC2).
> This tells PyDESeq2 to control for variety baseline differences when
> testing Root vs Leaf. All other species default to `DESIGN="condition"`.
> Override with `DESIGN="condition"` to force single-factor, or
> `DESIGN="variety+condition"` to force multi-factor. Requires a `variety`
> column in `sample_metadata.tsv` (auto-generated by `build_count_matrix.py`).

### Path A: Quick plots (MA, volcano, PCA — no gene family heatmaps)

| Stage | Command (DC) | Command (DG, leaf focus) |
|-------|-------------|------------------------|
| DESeq2 step 2 (filter) | `sbatch scripts/05_pydeseq2/run_step2_filter.sbatch DC` | `CONTRAST_A=L CONTRAST_B=R sbatch scripts/05_pydeseq2/run_step2_filter.sbatch DG` |
| DESeq2 step 3 (plots) | `sbatch scripts/05_pydeseq2/run_step3_plots.sbatch DC` | `CONTRAST_A=L CONTRAST_B=R sbatch scripts/05_pydeseq2/run_step3_plots.sbatch DG` |

### Stage 08 prerequisite: Build a gene-family database (run once per family)

Before running Path B or C you need a master list. This only needs to run
**once** per family — it creates the `<family>_master_list.csv` file.

| Family | Command |
|--------|---------|
| CYP | `sbatch scripts/08_gene_families/run_gene_database.sbatch cyp` |
| OMT | `sbatch scripts/08_gene_families/run_gene_database.sbatch omt` |

### Path B: Short gene-family heatmaps (HMMER list → intersect → plots)

Replace `cyp` / `omt` to switch families. Both use the same scripts.

**For CYP:**

| Step | Stage | Command |
|------|-------|---------|
| 1 | Intersect + extract | `sbatch scripts/08_gene_families/run_gene_extract.sbatch cyp` |
| 2 | DESeq2 step 3 (plots) | `sbatch scripts/05_pydeseq2/run_step3_plots.sbatch DC cyp_expressed` |

**For OMT:**

| Step | Stage | Command |
|------|-------|---------|
| 1 | Intersect + extract | `sbatch scripts/08_gene_families/run_gene_extract.sbatch omt` |
| 2 | DESeq2 step 3 (plots) | `sbatch scripts/05_pydeseq2/run_step3_plots.sbatch DC omt_expressed` |

### Path C: Full gene-family pipeline (HMMER list → BLAST → combine → plots)

Same idea — replace `cyp` / `omt` and update the FASTA path accordingly.

**For CYP:**

| Step | Stage | Command |
|------|-------|---------|
| 1 | Intersect + extract | `sbatch scripts/08_gene_families/run_gene_extract.sbatch cyp` |
| 2 | BLASTp (swissprot) | `sbatch scripts/06_blast/run_blastp_discovery.sbatch DC swissprot 07_NRdatabase/cyp450_database/cyp_proteins.fasta` |
| 3 | Combine + filter | `sbatch scripts/06_blast/run_combine_blast_deseq.sbatch DC swissprot discovery standard cyp` |
| 4 | DESeq2 step 3 (plots) | `sbatch scripts/05_pydeseq2/run_step3_plots.sbatch DC cyp_discovery` |

**For OMT:**

| Step | Stage | Command |
|------|-------|---------|
| 1 | Intersect + extract | `sbatch scripts/08_gene_families/run_gene_extract.sbatch omt` |
| 2 | BLASTp (swissprot) | `sbatch scripts/06_blast/run_blastp_discovery.sbatch DC swissprot 07_NRdatabase/omt_database/omt_proteins.fasta` |
| 3 | Combine + filter | `sbatch scripts/06_blast/run_combine_blast_deseq.sbatch DC swissprot discovery standard omt` |
| 4 | DESeq2 step 3 (plots) | `sbatch scripts/05_pydeseq2/run_step3_plots.sbatch DC omt_discovery` |

### Path B-full: All-gene BLAST plots (no family filter)

| Stage | Command |
|-------|---------|
| BLAST prep | `sbatch scripts/06_blast/run_blast_input.sbatch DC` |
| BLASTp | `sbatch scripts/06_blast/run_blastp_discovery.sbatch DC swissprot` |
| Combine BLAST+DESeq2 | `sbatch scripts/06_blast/run_combine_blast_deseq.sbatch DC swissprot discovery` |
| DESeq2 step 3 (plots) | `sbatch scripts/05_pydeseq2/run_step3_plots.sbatch DC combined_annotated` |

### Path 1 (DG): Gene list filter → plots → verify

```bash
# Filter gene list + run PyDESeq2 (leaf focus)
CONTRAST_A=L CONTRAST_B=R sbatch scripts/08_gene_families/run_filter_genelist.sbatch DG

# Generate plots from filtered candidates
CONTRAST_A=L CONTRAST_B=R sbatch scripts/05_pydeseq2/run_step3_plots.sbatch DG \
  /projects/tholl_lab_1/daisy_analysis/07_NRdatabase/sukman_database/geneious_candidates_DG.tsv

# Verify against previous student
CONTRAST_A=L CONTRAST_B=R \
RESULTS=06_analysis/pydeseq2_DG_step3_plots/cyp_gene_list.tsv \
RESULTS_UNFILTERED=06_analysis/pydeseq2_DG_step1_unfiltered/pydeseq2_results_UNFILTERED.tsv \
OUTPUT=07_NRdatabase/sukman_database/verification_DG_filtered.tsv \
OUTPUT_UNFILTERED=07_NRdatabase/sukman_database/verification_DG_unfiltered.tsv \
sbatch 05_rnaseq-code/scripts/11_verify/run_verify_genelist.sbatch
```

### Stages 07–10: Downstream analysis

| Stage | Command |
|-------|---------|
| HMMER | `sbatch scripts/07_domains/run_hmmer.sbatch DC` |
| Gene families | `sbatch scripts/08_gene_families/run_gene_families_combined.sbatch DC DG swissprot discovery full` |
| Compare species | `sbatch scripts/09_comparison/run_compare_species.sbatch DC DG swissprot discovery` |

### Stage 11: Verification against previous students

The verify script compares your PyDESeq2 results against a previous student's
P450 data. Set `DATABASE=` to choose which student to compare against:

| Database | What it contains | File(s) |
|----------|-----------------|---------|
| `sukman` (default) | 5 separate .txt files (gene list, DE, expression, log2fold, up/downregulated) | `07_NRdatabase/sukman_database/` |
| `ahmed` | Single CSV with upregulated P450 genes (full DESeq2 output) | `07_NRdatabase/ahmed_database/Upregulated_P450_values_ahmed.csv` |

**DC (root focus, default contrast R vs L):**

```bash
cd /projects/tholl_lab_1/daisy_analysis

# DC vs Sukman (default — just run with no extra variables)
sbatch 05_rnaseq-code/scripts/11_verify/run_verify_genelist.sbatch

# DC vs Ahmed
DATABASE=ahmed \
sbatch scripts/11_verify/run_verify_genelist.sbatch
```

**DG (leaf focus, flipped contrast L vs R):**

DG uses `CONTRAST_A=L CONTRAST_B=R` so positive log2FC = up in leaf.
You also need to point RESULTS to the DG output files and change the
output names to say DG instead of DC.

```bash
cd /projects/tholl_lab_1/daisy_analysis

# DG vs Sukman
CONTRAST_A=L CONTRAST_B=R \
RESULTS=06_analysis/pydeseq2_DG_step3_plots/cyp_gene_list.tsv \
RESULTS_UNFILTERED=06_analysis/pydeseq2_DG_step1_unfiltered/pydeseq2_results_UNFILTERED.tsv \
OUTPUT=07_NRdatabase/sukman_database/verify_DG_filtered_sukman.tsv \
OUTPUT_UNFILTERED=07_NRdatabase/sukman_database/verify_DG_unfiltered_sukman.tsv \
sbatch scripts/11_verify/run_verify_genelist.sbatch

# DG vs Ahmed
DATABASE=ahmed CONTRAST_A=L CONTRAST_B=R \
RESULTS=06_analysis/pydeseq2_DG_step3_plots/cyp_gene_list.tsv \
RESULTS_UNFILTERED=06_analysis/pydeseq2_DG_step1_unfiltered/pydeseq2_results_UNFILTERED.tsv \
OUTPUT=07_NRdatabase/ahmed_database/verify_DG_filtered_ahmed.tsv \
OUTPUT_UNFILTERED=07_NRdatabase/ahmed_database/verify_DG_unfiltered_ahmed.tsv \
sbatch scripts/11_verify/run_verify_genelist.sbatch
```

**SK (root focus, default contrast R vs L):**

```bash
cd /projects/tholl_lab_1/daisy_analysis

# SK vs Ahmed
DATABASE=ahmed \
RESULTS=06_analysis/pydeseq2_SK_step3_.../cyp_gene_list.tsv \
RESULTS_UNFILTERED=06_analysis/pydeseq2_SK_step1/pydeseq2_results_UNFILTERED.tsv \
OUTPUT=07_NRdatabase/ahmed_database/verify_SK_filtered_ahmed.tsv \
OUTPUT_UNFILTERED=07_NRdatabase/ahmed_database/verify_SK_unfiltered_ahmed.tsv \
sbatch scripts/11_verify/run_verify_genelist.sbatch
```

Output files (names end with which database was compared):

| Sample | Database | Filtered table | Summary |
|--------|----------|---------------|---------|
| DC | sukman | `verify_DC_filtered_sukman.tsv` | `verify_DC_filtered_SUMMARY_sukman.txt` |
| DC | ahmed | `verify_DC_filtered_ahmed.tsv` | `verify_DC_filtered_SUMMARY_ahmed.txt` |
| DG | sukman | `verify_DG_filtered_sukman.tsv` | `verify_DG_filtered_SUMMARY_sukman.txt` |
| DG | ahmed | `verify_DG_filtered_ahmed.tsv` | `verify_DG_filtered_SUMMARY_ahmed.txt` |

Six checks are run automatically (Check 6 is skipped for Ahmed since he
only has upregulated genes). The script automatically accounts for the
flipped contrast direction when comparing DG (L vs R) against previous
student data (R vs L). See comments in `run_verify_genelist.sbatch` for
details on each check.

---

## Changelog

### 2026-03-31: OMT support + script renames (cyp → gene)

**OMT support.** The gene-family pipeline now supports both CYP (cytochrome
P450) and OMT (O-methyltransferase) families. Every script now uses a single
codebase — pass `cyp` or `omt` as an argument to choose which family to run.
Both Path B (short) and Path C (full) work for either family.

**Script renames.** CYP-specific script names generalized to `gene_*`:

| Old name | New name |
|----------|----------|
| `run_cyp450_database.sbatch` | `run_gene_database.sbatch` |
| `run_cyp_extract.sbatch` | `run_gene_extract.sbatch` |
| `cyp_gff_keyword_search.py` | `gene_gff_keyword_search.py` |
| `cyp_hmmer_scan.py` | `gene_hmmer_scan.py` |
| `cyp_combine_results.py` | `gene_combine_results.py` |
| `cyp_intersect_pydeseq2.py` | `gene_intersect_pydeseq2.py` |
| `cyp_extract_proteins.py` | `gene_extract_proteins.py` |

**BLAST OMT parsing.** `combine_blast_deseq.py` now detects OMT genes from
BLAST descriptions (methyltransferase, COMT, CCoAOMT, etc.) and saves a
separate `*_OMT_only.tsv` alongside the existing `*_CYP_only.tsv`.

**New input sources for Step 3 plots:** `omt_expressed`, `omt_discovery`,
`omt_strict` — parallel to the existing `cyp_*` sources.

| File | Change |
|------|--------|
| `scripts/08_gene_families/run_gene_database.sbatch` | Renamed from `run_cyp450_database.sbatch`. Now accepts `FAMILY` arg (`cyp` or `omt`). Selects Pfam domain + keywords + output dir automatically. Replaces the separate `run_omt_database.sbatch`. |
| `scripts/08_gene_families/gene_gff_keyword_search.py` | Renamed from `cyp_gff_keyword_search.py`. Uses `--family` arg to pick CYP or OMT keywords. |
| `scripts/08_gene_families/gene_hmmer_scan.py` | Renamed from `cyp_hmmer_scan.py`. Uses `--pfam` and `--family-name` args for any Pfam domain. |
| `scripts/08_gene_families/gene_combine_results.py` | Renamed from `cyp_combine_results.py`. Uses `--family-name` arg for display labels. |
| `scripts/08_gene_families/run_gene_extract.sbatch` | Renamed from `run_cyp_extract.sbatch`. Now accepts `FAMILY` arg (`cyp` or `omt`). Routes to `cyp450_database/` or `omt_database/`. |
| `scripts/08_gene_families/gene_intersect_pydeseq2.py` | Renamed from `cyp_intersect_pydeseq2.py`. Updated docstring examples. |
| `scripts/08_gene_families/gene_extract_proteins.py` | Renamed from `cyp_extract_proteins.py`. Generalized banner, auto-detects family from path. |
| `scripts/06_blast/combine_blast_deseq.py` | `parse_blast_stitle()` now detects OMT subfamilies. Summary reports both CYP + OMT counts. Saves `*_OMT_only.tsv`. |
| `scripts/06_blast/run_combine_blast_deseq.sbatch` | `SOURCE` arg accepts `omt` alongside `cyp`. Routes to `omt_database/`. |
| `scripts/05_pydeseq2/run_step3_plots.sbatch` | Added `omt_expressed`, `omt_discovery`, `omt_strict` input sources. |

Look for `# CHANGED: 2026-03-31` comments in those files for exact locations.

### 2026-03-31: Naming convention, variety covariate, and MF fix

**Naming convention.** Documented a standard FASTQ naming format in
`config.sh` Section 4 and README "Adding New Samples":

```
{VARIETY}_{CONDITION}_{REPLICATE}_{READ}.fq.gz
   DC1   _    L     _     1     _ 1  .fq.gz
```

New samples should follow this format. Legacy names (DC1L1, DGL1, MFF1,
T_L_R1) are still auto-detected — no need to rename existing files.

**Multi-factor design.** Step 1 now auto-detects `DESIGN="variety+condition"`
for DC (which has DC1 + DC2). This tells PyDESeq2 to control for variety
baseline differences when testing Root vs Leaf. Other species default to
`DESIGN="condition"`. Override with `DESIGN="condition" sbatch ...` to force
single-factor.

**MF fix.** The condition regex now handles `F` (fruit) and `N` in addition
to `L` (leaf) and `R` (root). MF also gets `SAMPLE_CONDITIONS="MFF=F MFL=L"`
in config.sh as a fallback.

| File | Change |
|------|--------|
| `scripts/config.sh` | Added naming convention docs to Section 4. Added `SAMPLE_CONDITIONS="MFF=F MFL=L"` for MF. |
| `scripts/04_counting/build_count_matrix.py` | `extract_sample_info()` now supports two patterns: new underscore format (`DC1_L_1`) tried first, then legacy format (`DC1L1`). Regex expanded from `[LR]` to `[LRFN]`. Added `variety` column to both `extract_sample_info()` and `extract_sample_info_from_map()`. |
| `scripts/05_pydeseq2/run_step1_analysis.sbatch` | Added `DESIGN` variable (default: `auto`). Auto-detects `variety+condition` for DC, `condition` for all others. Passes `--design` to the Python script. |
| `README.md` | Replaced "Adding New Samples" with naming convention + config.sh steps. Added multi-factor design docs and command to Stage 05 table. |

### 2026-03-30: Simplify Step 1 PyDESeq2 + auto-fix verify output names

**Step 1 now uses a single `dds.deseq2()` call** instead of a 3-path
if/elif/else fallback. This matches how `R_pydeq2.py` and
`filter_count_by_genelist.py` already worked. `refit_cooks=True` is still
enabled. The old code is preserved in comments in case you need to restore it.

The verify script now **auto-detects** when you pass an unfiltered RESULTS
file and renames the output from `filtered` to `unfiltered` automatically.

Run Step 1 for DC:

```bash
cd /projects/tholl_lab_1/daisy_analysis
sbatch 05_rnaseq-code/scripts/05_pydeseq2/run_step1_analysis.sbatch DC
```

A new reference document `QUICK_PYDESEQ2_COMMANDS.md` explains how log2FC
is calculated across all five scripts and which settings differ between them.

| File | Change |
|------|--------|
| `scripts/05_pydeseq2/pydeseq2_run_analysis.py` | Replaced 3-path if/elif/else (deseq2/fit/step-by-step) with single `dds.deseq2()` call. Old code preserved in comments dated 2026-03-30. |
| `scripts/11_verify/run_verify_genelist.sbatch` | Added auto-fix: if RESULTS path contains "unfiltered" but OUTPUT says "filtered", the output name is corrected automatically. |
| `QUICK_PYDESEQ2_COMMANDS.md` | New file documenting how log2FC is calculated in Step 1, Step 2, Step 3, filter_count_by_genelist.py, and R_pydeq2.py. |

### 2026-03-30: Ahmed database support + unified verify script

The verification script now supports comparing against multiple databases
(Sukman or Ahmed) from a single sbatch file using `DATABASE=sukman` or
`DATABASE=ahmed`.

| File | Change |
|------|--------|
| `scripts/11_verify/run_verify_genelist.sbatch` | Added `DATABASE` variable with if/elif to switch between sukman and ahmed file paths. Output filenames now end with the database name (e.g. `verify_DC_filtered_ahmed.tsv`). Passes `--database-label` to the Python script. |
| `scripts/11_verify/verify_genelist.py` | Checks 5 and 6 now use `load_ids_from_list()` instead of manual line reading, so CSV files (like Ahmed's) work correctly. Added `--database-label` argument so the summary file ends with `_SUMMARY_ahmed.txt`. |

Look for `# CHANGED:` comments in those files for exact locations.

To add a new database in the future, copy the `ahmed` elif block in the
sbatch and change the paths. See the comments in the file for instructions.

### 2026-03-28: Contrast direction bugfix (L vs R support)

Three scripts had hardcoded `R vs L` assumptions that caused wrong labels and
verification results when running DG with `CONTRAST_A=L CONTRAST_B=R`:

| File | Fix |
|------|-----|
| `scripts/08_gene_families/filter_count_by_genelist.py` | Print labels now use the actual contrast instead of hardcoded "root"/"leaf" |
| `scripts/11_verify/verify_genelist.py` | Added `--contrast-a`/`--contrast-b` args; checks 4, 5, 6 and the comparison table now flip direction when contrast is L vs R |
| `scripts/11_verify/run_verify_genelist.sbatch` | Now passes `--contrast-a`/`--contrast-b` to the Python script (was reading env vars but not forwarding them) |

Look for `# CHANGED:` comments in those files for exact locations. The core
PyDESeq2 calculation, heatmaps, MA plots, and volcano plots were already correct.

---

## Adding New Samples

Two things to do: (1) name your FASTQ files correctly, (2) register the
sample group in `scripts/config.sh`.

### Step 1: Name your FASTQ files

Follow this convention for all new paired-end FASTQ files:

```
{VARIETY}_{CONDITION}_{REPLICATE}_{READ}.fq.gz
```

| Part | What it is | Examples |
|------|-----------|---------|
| VARIETY | Species code + variety number | DC1, DC2, DG, MF |
| CONDITION | Single letter for tissue | L (leaf), R (root), F (fruit) |
| REPLICATE | Replicate number | 1, 2, 3 |
| READ | Paired-end read | 1 or 2 |

**One variety** (most species):

```
DG_L_1_1.fq.gz   DG_L_1_2.fq.gz     ← DG, Leaf, replicate 1, reads 1 & 2
DG_L_2_1.fq.gz   DG_L_2_2.fq.gz     ← DG, Leaf, replicate 2
DG_R_1_1.fq.gz   DG_R_1_2.fq.gz     ← DG, Root, replicate 1
```

**Multiple varieties** (like DC with DC1 and DC2):

```
DC1_L_1_1.fq.gz  DC1_L_1_2.fq.gz    ← DC variety 1, Leaf, replicate 1
DC2_L_1_1.fq.gz  DC2_L_1_2.fq.gz    ← DC variety 2, Leaf, replicate 1
DC1_R_1_1.fq.gz  DC1_R_1_2.fq.gz    ← DC variety 1, Root, replicate 1
```

If you follow this format, the pipeline auto-detects variety, condition,
and replicate — no extra configuration needed.

> **Legacy names:** Older samples (DC1L1, DGL1, MFF1, T_L_R1) are still
> supported via regex auto-detection or `SAMPLE_CONDITIONS` in config.sh.
> You don't need to rename existing files.

### Step 2: Register the sample group in config.sh

Open `scripts/config.sh` and add a new block inside the `get_sample_info()`
function:

```bash
XX)
    SPECIES_DIR="00_5_NewName"
    SPECIES_NAME="Species name"
    GENOME_TYPE="carrot"
    ;;
```

Replace `XX` with a short code (2-4 letters), `SPECIES_DIR` with your data
folder name under `00_rawdata/`, and `GENOME_TYPE` with `carrot` or `nutmeg`.

If your files DON'T follow the naming convention (legacy/external data), add
`SAMPLE_CONDITIONS` to tell the pipeline how to parse them:

```bash
XX)
    SPECIES_DIR="00_5_NewName"
    SPECIES_NAME="Species name"
    GENOME_TYPE="carrot"
    SAMPLE_CONDITIONS="_L_=L _R_=R"
    ;;
```

The format is `substring=condition` — if a sample name contains the substring,
it gets that condition.

### Step 3: Rebuild the metadata (if you have multiple varieties)

If your sample group contains multiple varieties (e.g. DC1 and DC2),
`build_count_matrix.py` now auto-generates a `variety` column from the
sample names. Re-run the count matrix step to regenerate the metadata:

```bash
# Delete old metadata so it gets rebuilt
rm 03_count_tables/00_1_DC/sample_metadata.tsv

# Re-run (will regenerate gene_count_matrix.tsv + sample_metadata.tsv)
sbatch scripts/04_counting/run_count_matrix.sbatch DC
```

The resulting `sample_metadata.tsv` will look like:

```
sample  group  group_number  condition  replicate  variety  treatment  full_condition
DC1L1   DC     1             L          1          DC1      DC1        DC1_L
DC1L2   DC     1             L          2          DC1      DC1        DC1_L
DC1L3   DC     1             L          3          DC1      DC1        DC1_L
DC1R1   DC     1             R          1          DC1      DC1        DC1_R
DC1R2   DC     1             R          2          DC1      DC1        DC1_R
DC1R3   DC     1             R          3          DC1      DC1        DC1_R
DC2L1   DC     2             L          1          DC2      DC2        DC2_L
DC2L2   DC     2             L          2          DC2      DC2        DC2_L
DC2L3   DC     2             L          3          DC2      DC2        DC2_L
DC2R1   DC     2             R          1          DC2      DC2        DC2_R
DC2R2   DC     2             R          2          DC2      DC2        DC2_R
DC2R3   DC     2             R          3          DC2      DC2        DC2_R
```

Step 1 auto-detects `DESIGN="variety+condition"` for DC (see Stage 05
commands above). For species with only one variety, it uses `condition`.

---

## Documentation Guide

| Document | What it covers |
|----------|---------------|
| [PIPELINE_WORKFLOW.md](docs/PIPELINE_WORKFLOW.md) | Step-by-step commands, inputs, outputs, and a progress checklist |
| [QUICK_PYDESEQ2_COMMANDS.md](QUICK_PYDESEQ2_COMMANDS.md) | How log2FC is calculated across all 5 scripts and what settings differ |
| [HMMER_PROSITE_GUIDE.md](docs/HMMER_PROSITE_GUIDE.md) | Setting up and running domain/motif analysis |
| [CONCEPTS.md](docs/CONCEPTS.md) | RNA-seq background: FASTQ, QC, alignment, DE, BLAST |
| [SETUP.md](docs/SETUP.md) | Conda environment, HPC layout, first-time configuration |
| [DIRECTORY_STRUCTURE.md](docs/DIRECTORY_STRUCTURE.md) | Complete file tree with descriptions |

---

## HPC Data Layout

```
/projects/tholl_lab_1/daisy_analysis/
├── 00_rawdata/        FASTQ files (DC, DG, MF sample groups)
├── 01_processed/      QC reports, fastp output, Trinity assemblies
├── 02_mapped/         STAR BAM files + ReadsPerGene.out.tab
├── 03_count_tables/   gene_count_matrix.tsv + sample_metadata.tsv per species
├── 04_reference/      Genome FASTA, GTF, STAR indices, protein FAA
├── 05_rnaseq-code/    This repository (git clone)
├── 06_analysis/       All analysis output (DESeq2, BLAST, HMMER, plots)
└── 07_NRdatabase/     BLAST databases, gene-family databases
    ├── cyp450_database/   CYP master list, expressed list, proteins
    ├── omt_database/      OMT master list, expressed list, proteins
    ├── sukman_database/   Sukman's P450 reference files (.txt)
    └── ahmed_database/    Ahmed's upregulated P450 CSV
```
