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
│   ├── 04_busco/             │  BUSCO transcriptome completeness QC
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
 01 QC ──→ 01 RiboDetector ──→ 02 Alignment ──→ 03 Assembly ──→ 04 Counting
           (MF only — removes              │
            residual rRNA before           ▼
            alignment)                04 BUSCO
                                      (transcriptome
                                       completeness QC —
                                       run after Trinity)
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
| Remove rRNA | `sbatch scripts/01_qc/run_ribodetector.sbatch MF` | **MF only** — removes residual rRNA detected after STAR. Run after fastp, before STAR. Output goes to `01_processed/00_3_MF/ribofree/`. Install once: `pip install ribodetector` |
| STAR index | `sbatch scripts/02_alignment/run_genome_index.sbatch carrot` | Build once per genome |
| Align | `sbatch scripts/02_alignment/run_star_align_all.sbatch MF` | Auto-detects ribofree > clean > raw reads |
| Assembly | See **Trinity assembly** section below | Recommended: array (parallel) |
| Count reads + build matrix | `sbatch scripts/04_counting/run_featurecounts.sbatch DC` | |
| Rebuild metadata only | `sbatch scripts/04_counting/run_count_matrix.sbatch DC` | |

#### Trinity assembly (`scripts/03_assembly/`)

**Recommended: use the array script** — each sample gets its own node and
72-hour time limit, so they run in parallel and finish in ~2 days instead
of ~7 days sequentially.

```bash
# Step 1: Generate sample list (auto-detects clean vs raw reads)
bash scripts/03_assembly/run_trinity_array.sbatch --generate MF

# Step 2: Submit all samples as parallel jobs (replace 6 with your count)
sbatch --array=1-6 scripts/03_assembly/run_trinity_array.sbatch
```

If a job times out, just resubmit — Trinity resumes from checkpoint
automatically. Completed samples are skipped.

| Script | What it does | When to use |
|--------|-------------|-------------|
| **`run_trinity_array.sbatch`** | Each sample runs as an independent parallel job on its own node. **72h per sample.** | **Recommended.** Fastest option (~2 days total). |
| `run_trinity_all.sbatch` | All samples run sequentially in one job. Accepts optional species code (e.g., `MF`). | When you can't use multiple nodes, or only have 1-2 samples. |
| `run_trinity.sbatch` | Runs one specific sample (e.g., `sbatch run_trinity.sbatch 00_3_MF MFF1`). | Test a single sample first, or re-run one that failed. |
| `scripts/04_rsem/run_rsem.sbatch` | Quantifies expression (RSEM + Bowtie2) against pooled/CD-HIT Trinity. See [docs/README_trinity.md](docs/README_trinity.md). | **After** Trinity + CD-HIT. Produces `RSEM.gene.counts.matrix`. |
| `scripts/04_counting/run_count_matrix.sbatch MF --type rsem` | Converts RSEM output to `gene_count_matrix.tsv` + `sample_metadata.tsv` in `03_count_tables/00_5_MF_trinity/`. | **After** RSEM. Required before PyDESeq2 for Trinity projects. |

All assembly scripts auto-detect the best available reads in this priority order:

1. **RiboDetector-filtered reads** (`ribofree/`) — used if present (rRNA removed)
2. **Fastp-cleaned reads** (`*_clean.fq.gz`) — used if no ribofree reads
3. **Raw reads** (`00_rawdata/`) — fallback if nothing else exists

For MF samples, ribofree reads will be picked up automatically after running
`run_ribodetector.sbatch`. DC and DG use fastp-cleaned reads as usual.

### Stage 04: BUSCO Transcriptome Completeness QC (`scripts/04_busco/`)

After Trinity finishes, run **BUSCO** to check how complete your assembly is.
BUSCO compares your transcriptome against a curated set of "universal" genes
that should be present in every member of a taxonomic group and reports the
percentage found (Complete / Duplicated / Fragmented / Missing).

The script auto-picks the right lineage based on the species code:

| Species code | Genome type | BUSCO lineage | # genes | Why |
|--------------|-------------|---------------|---------|-----|
| `MF` | nutmeg | `embryophyta_odb10` | 1,614 | Nutmeg is a basal angiosperm — use land plants |
| `DC`, `DG`, `SK`, `DCDG` | carrot | `eudicots_odb10` | ~2,326 | Carrot is a eudicot — use the tighter set |

**One-time setup** — create a dedicated `busco_env` conda env (BUSCO doesn't
support Python 3.14 yet, so it gets its own env):

```bash
# Run ONCE on a login node — needs internet to download packages
bash scripts/04_busco/setup_busco_env.sh
```

This helper is idempotent (safe to re-run), auto-detects mamba for faster
solves, and verifies that BUSCO + HMMER + AUGUSTUS + BLAST+ + SEPP all work.

**Run BUSCO on your pooled Trinity assembly:**

| Species | Command |
|---------|---------|
| Nutmeg (MF) | `sbatch scripts/04_busco/run_busco.sbatch MF` |
| Carrot (DC) | `sbatch scripts/04_busco/run_busco.sbatch DC` |
| Carrot (DG) | `sbatch scripts/04_busco/run_busco.sbatch DG` |
| Carrot (SK) | `sbatch scripts/04_busco/run_busco.sbatch SK` |
| Combined (DCDG) | `sbatch scripts/04_busco/run_busco.sbatch DCDG` |

> **Important:** Before your FIRST BUSCO submission, create the output folder
> on the HPC (the `--chdir` line in the .sbatch header needs it to exist):
> ```bash
> mkdir -p /projects/tholl_lab_1/daisy_analysis/01_processed/00_5_BUSCO
> ```

**Where outputs go.** Each species writes to its own subfolder under
`01_processed/00_5_BUSCO/`, so runs don't overwrite each other:

```
01_processed/00_5_BUSCO/
├── busco_<jobid>.out                      ← SLURM stdout log
├── busco_<jobid>.err                      ← SLURM stderr log
├── busco_MF/busco_MF_results/
│   ├── short_summary*.txt                 ← read this first (C/D/F/M %)
│   ├── full_table.tsv                     ← every BUSCO gene + status
│   ├── missing_busco_list.tsv             ← genes BUSCO couldn't find
│   └── logs/busco.log                     ← full log (for debugging)
├── busco_DC/busco_DC_results/
└── busco_DG/busco_DG_results/
```

**Reading the result.** Open the `short_summary*.txt` — example for a good
plant transcriptome:

```
C:84.5%[S:21.3%,D:63.2%],F:6.1%,M:9.4%,n:1614
```

- **C > 70%** (ideally > 85%) — good completeness
- **High D% is NORMAL** for transcriptomes — Trinity makes multiple isoforms
  per gene, so the same BUSCO gene shows up several times
- **M < 20%** — acceptable missing-gene rate

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
| HMMER (carrot reference) | `sbatch scripts/07_domains/run_hmmer.sbatch DC` |
| HMMER (Trinity MF nutmeg — CYP/OMT) | `sbatch scripts/07_domains/run_hmmer_genefinder.sbatch MF` |
| BLAST tblastn (Trinity MF nutmeg — CYP/OMT) | `sbatch scripts/07_domains/run_blast_genefinder.sbatch MF` |
| Gene families | `sbatch scripts/08_gene_families/run_gene_families_combined.sbatch DC DG swissprot discovery full` |
| Compare species | `sbatch scripts/09_comparison/run_compare_species.sbatch DC DG swissprot discovery` |

> **Trinity nutmeg (MF):** For CYP/OMT discovery, use the targeted Trinity scripts
> instead of the full-genome carrot pipeline: BLAST tblastn (`run_blast_genefinder.sbatch`,
> after CD-HIT) and HMMER (`run_hmmer_genefinder.sbatch`, after TransDecoder). Combine
> both ID lists (union). See [docs/README_trinity.md](docs/README_trinity.md) Steps 3b, 9, and 10.

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
