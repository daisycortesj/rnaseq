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
 01 QC ──→ 02 Alignment ──→ 03 Assembly ──→ 04 Counting ──→ 05 PyDESeq2
                                                                │
 06 BLAST ←─────────────────────────────────────────────────────┘
     │
 07 Domains (HMMER / PROSITE)
     │
 08 Gene Families (CYP, OMT extraction)
     │
 09 Cross-Species Comparison (DC vs DG)
     │
 10 Post-Analysis (phylogenetic trees, genomic clustering)
     │
 11 Verification
```

Each numbered directory in `scripts/` matches a pipeline stage. Run them in
order. Each directory contains `.sbatch` files you submit with `sbatch` and
`.py` helper scripts that the batch jobs call automatically.

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

| Stage | Command |
|-------|---------|
| QC | `sbatch scripts/01_qc/run_fastqc.sbatch` |
| Trim reads | `sbatch scripts/01_qc/run_fastp.sbatch` |
| STAR index | `sbatch scripts/02_alignment/run_genome_index.sbatch carrot` |
| Align | `sbatch scripts/02_alignment/run_star_align_all.sbatch` |
| Count matrix | `sbatch scripts/04_counting/run_count_matrix.sbatch DC` |
| DESeq2 step 1 | `sbatch scripts/05_pydeseq2/run_step1_analysis.sbatch DC` |
| DESeq2 step 2 | `sbatch scripts/05_pydeseq2/run_step2_filter.sbatch DC` |
| DESeq2 step 3 | `sbatch scripts/05_pydeseq2/run_step3_plots.sbatch DC` |
| BLAST prep | `sbatch scripts/06_blast/run_blast_input.sbatch DC` |
| BLASTp | `sbatch scripts/06_blast/run_blastp_discovery.sbatch DC swissprot` |
| HMMER | `sbatch scripts/07_domains/run_hmmer.sbatch DC` |
| Gene families | `sbatch scripts/08_gene_families/run_gene_families_combined.sbatch DC DG swissprot discovery full` |
| Compare species | `sbatch scripts/09_comparison/run_compare_species.sbatch DC DG swissprot discovery` |

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
```
