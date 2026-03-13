# Environment & HPC Setup

One-time setup for running the RNA-seq pipeline on Virginia Tech ARC.

---

## 1. Conda Environment

The `environment.yml` at the repository root defines every tool and library.

```bash
# Create the environment (first time)
conda env create -f environment.yml

# Activate before any work
conda activate rnaseq

# Update after pulling new changes
conda env update -f environment.yml --prune
```

### What's installed

| Category | Packages |
|----------|----------|
| Core tools | STAR, Trinity, RSEM, Bowtie2, samtools |
| QC | FastQC, MultiQC, fastp |
| Python | pandas, numpy, matplotlib, seaborn, scipy, scikit-learn |
| Sequence analysis | BLAST, HMMER, EMBOSS (PROSITE), Entrez Direct |

PyDESeq2 is installed via pip inside the environment:

```bash
pip install pydeseq2
```

### Install the pipeline package

```bash
cd /path/to/rnaseq
pip install -e .
```

This makes `rna_pipeline` importable and provides the `rna-pipeline` CLI
command (`rna-pipeline --mode qc`, etc.).

---

## 2. HPC Data Directory

All data lives under `/projects/tholl_lab_1/daisy_analysis/`. The repository
itself is cloned into `05_rnaseq-code/`.

```
/projects/tholl_lab_1/daisy_analysis/
├── 00_rawdata/        Raw FASTQ files organized by species/sample group
│   ├── DC/            Daucus carota
│   ├── DG/            Daucus glaber
│   └── MF/            M. fistulosa (if applicable)
├── 01_processed/      fastp trimmed reads, FastQC/MultiQC reports
├── 02_mapped/         STAR alignment output (BAM, ReadsPerGene.out.tab)
├── 03_count_tables/   Count matrices built from STAR or featureCounts
├── 04_reference/      Reference genomes, GTF annotations, STAR indices
├── 05_rnaseq-code/    ← This repo
├── 06_analysis/       All downstream output (DESeq2, BLAST, heatmaps, etc.)
└── 07_NRdatabase/     BLAST databases (swissprot, NR), CYP450 reference DB
```

Most sbatch scripts define paths relative to a `BASE_DIR` variable set to
`/projects/tholl_lab_1/daisy_analysis` and a `CODE_DIR` variable pointing to
`05_rnaseq-code`.

---

## 3. SLURM Basics

All computation runs through SLURM batch jobs:

```bash
# Submit a job
sbatch scripts/01_qc/run_fastqc.sbatch

# Check your jobs
squeue -u $USER

# Cancel a job
scancel <JOB_ID>

# View job output (logs go to the directory specified by #SBATCH -o)
less slurm-<JOB_ID>.out
```

### Common SLURM parameters in our scripts

```bash
#SBATCH --account=tholl_lab_1         # Allocation account
#SBATCH --partition=normal_q          # Queue (normal_q, dev_q, etc.)
#SBATCH --cpus-per-task=8             # CPU cores
#SBATCH --mem=32G                     # RAM
#SBATCH --time=04:00:00               # Wall time limit
```

### Job dependencies

Chain pipeline steps so step N waits for step N-1:

```bash
JOB1=$(sbatch --parsable scripts/04_counting/run_featurecounts.sbatch)
sbatch --dependency=afterok:$JOB1 scripts/04_counting/run_count_matrix.sbatch DC
```

---

## 4. One-Time Database Setup

Before running BLAST or CYP450 analysis, download the required databases:

```bash
# SwissProt + NR databases
sbatch scripts/setup/download_blastdb.sbatch

# Protein reference databases
sbatch scripts/setup/download_protein_databases.sbatch

# Daucus carota protein DB extraction
sbatch scripts/setup/extract_daucus_carota_db.sbatch
```

### HMMER Pfam database

```bash
conda activate rnaseq
cd /projects/tholl_lab_1/daisy_analysis/07_NRdatabase
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
```

See [HMMER_PROSITE_GUIDE.md](HMMER_PROSITE_GUIDE.md) for full domain analysis
setup including PROSITE.

---

## 5. Verifying Your Setup

```bash
conda activate rnaseq

# Check key tools
which STAR && STAR --version
which fastqc && fastqc --version
which hmmscan
python -c "import pydeseq2; print('PyDESeq2 OK')"
python -c "from rna_pipeline.cli import main; print('rna_pipeline OK')"
```

If any command fails, re-run `conda env update -f environment.yml --prune`
and `pip install -e .`.
