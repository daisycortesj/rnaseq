# RNA-seq Differential Expression & Gene Family Pipeline

Analysis pipeline for CYP (cytochrome P450) and OMT (O-methyltransferase) gene
families in *Daucus carota* (DC) and *Daucus glaber* (DG). Designed for
Virginia Tech ARC HPC with SLURM.

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
│   ├── 08_gene_families/     │  CYP450 database, gene family extraction
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
                              ┌────────────────────┤
                              │                    │
                         PATH A (quick)       PATH B (full — includes CYP/OMT heatmaps)
                              │                    │
                              ▼                    ▼
                     05 Step 2 (filter)    06 BLAST (annotate genes)
                              │                    │
                              ▼                    ▼
                     05 Step 3 (plots)     06 Combine BLAST + DESeq2
                     MA, volcano, PCA,             │
                     correlation only              ▼
                                           05 Step 3 (plots)
                                           MA, volcano, PCA, correlation,
                                           + CYP heatmap, OMT heatmap
                                                   │
                                                   ▼
                                           07 Domains (HMMER / PROSITE)
                                                   │
                                           08 Gene Families (CYP, OMT)
                                                   │
                                           09 Cross-Species Comparison
                                                   │
                                           10 Post-Analysis
                                                   │
                                           11 Verification
```

**Why two paths?** PyDESeq2 Step 1 must run before BLAST (BLAST needs the
gene list that Step 1 produces). But Step 3 plots need BLAST annotations to
make CYP/OMT heatmaps. So the pipeline loops back: Step 1 → BLAST → Step 3.

If you only need MA plot, volcano, PCA, and sample correlation, Path A is
enough. If you need gene family heatmaps, you must take Path B through BLAST
first.

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

| Stage | Command |
|-------|---------|
| QC | `sbatch scripts/01_qc/run_fastqc.sbatch` |
| Trim reads | `sbatch scripts/01_qc/run_fastp.sbatch` |
| STAR index | `sbatch scripts/02_alignment/run_genome_index.sbatch carrot` |
| Align | `sbatch scripts/02_alignment/run_star_align_all.sbatch` |
| Count matrix | `sbatch scripts/04_counting/run_count_matrix.sbatch DC` |

### Stage 05: PyDESeq2 Step 1 (always run this first)

| Stage | Command |
|-------|---------|
| DESeq2 step 1 (DC, root focus) | `sbatch scripts/05_pydeseq2/run_step1_analysis.sbatch DC` |
| DESeq2 step 1 (DG, leaf focus) | `CONTRAST_A=L CONTRAST_B=R sbatch scripts/05_pydeseq2/run_step1_analysis.sbatch DG` |

> **Contrast direction:** Default is `R vs L` (positive log2FC = up in root).
> For DG leaf-focus analysis, prepend `CONTRAST_A=L CONTRAST_B=R` to flip the
> contrast so positive log2FC = up in leaf. Use the **same contrast** for all
> steps in a run. See [CYP_heatmap.md](scripts/04_counting/CYP_heatmap.md)
> for full DC and DG command examples.

### Path A: Quick plots (MA, volcano, PCA — no gene family heatmaps)

| Stage | Command (DC) | Command (DG, leaf focus) |
|-------|-------------|------------------------|
| DESeq2 step 2 (filter) | `sbatch scripts/05_pydeseq2/run_step2_filter.sbatch DC` | `CONTRAST_A=L CONTRAST_B=R sbatch scripts/05_pydeseq2/run_step2_filter.sbatch DG` |
| DESeq2 step 3 (plots) | `sbatch scripts/05_pydeseq2/run_step3_plots.sbatch DC` | `CONTRAST_A=L CONTRAST_B=R sbatch scripts/05_pydeseq2/run_step3_plots.sbatch DG` |

### Path B: Full plots with CYP/OMT heatmaps (requires BLAST)

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

When you have new samples to analyze, you need to register them in
`scripts/config.sh` (Section 4). This is the **only file** you need to edit.

### Step 1: Add a sample group block

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

### Step 2: Add SAMPLE_CONDITIONS (if your sample names are non-standard)

The pipeline auto-detects conditions from sample names like `DC1L1` (it finds
the `L` = Leaf). If your sample names use a different format (e.g., `T_L_R1`,
`CtrlA_1`), add one extra line telling the pipeline how to read them:

```bash
XX)
    SPECIES_DIR="00_5_NewName"
    SPECIES_NAME="Species name"
    GENOME_TYPE="carrot"
    SAMPLE_CONDITIONS="_L_=L _R_=R"
    ;;
```

The format is `substring=condition` — if a sample name contains the substring,
it gets that condition. Examples:

| Sample names | SAMPLE_CONDITIONS line |
|---|---|
| `T_L_R1`, `T_R_R2` | `SAMPLE_CONDITIONS="_L_=L _R_=R"` |
| `CtrlA1`, `TreatA1` | `SAMPLE_CONDITIONS="Ctrl=C Treat=T"` |
| `NormalS1`, `EthyleneS1` | `SAMPLE_CONDITIONS="Normal=N Ethylene=E"` |
| `DC1L1`, `DC1R1` | *(not needed — auto-detected)* |

If you skip this line and the auto-detection fails, the pipeline will produce
a metadata file with empty condition fields, and PyDESeq2 will crash. See
`config.sh` Section 4 comments for more details.

---

## Documentation Guide

| Document | What it covers |
|----------|---------------|
| [PIPELINE_WORKFLOW.md](docs/PIPELINE_WORKFLOW.md) | Step-by-step commands, inputs, outputs, and a progress checklist |
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
└── 07_NRdatabase/     BLAST databases, CYP450 database
    ├── sukman_database/   Sukman's P450 reference files (.txt)
    └── ahmed_database/    Ahmed's upregulated P450 CSV
```
