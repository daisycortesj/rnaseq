# Trinity Assembly + Quality Check Workflow

**Scripts:** `scripts/03_assembly/trinity_pooled.sbatch`, `scripts/03_assembly/run_trinity.sbatch`  
**Related scripts:** `scripts/04_busco/run_busco.sbatch`, `scripts/04_busco/run_cdhit.sbatch`, `scripts/04_busco/run_busco_longest_isoform.sbatch`, `scripts/04_rsem/run_rsem.sbatch`  
**Where it fits:** fastp → (optional Kraken2) → **Trinity (here)** → BUSCO → CD-HIT → BUSCO again → **RSEM**

All paths below use `BASE_DIR=/projects/tholl_lab_1/daisy_analysis` from
`scripts/config.sh`. Submit jobs from the repo on ARC:

```bash
cd /projects/tholl_lab_1/daisy_analysis/05_rnaseq-code
```

---

## Why pooled Trinity?

Differential expression (DESeq2) needs **the same transcript IDs across all
samples**. If you run Trinity separately on each sample, you get different
assemblies with incompatible IDs — you cannot compare counts.

**Pooled assembly** = combine all reads from one species into one Trinity run.
That gives you one reference transcriptome. You then quantify each sample
against that shared assembly (e.g. with RSEM).

| Approach | Script | When to use |
|----------|--------|-------------|
| **Pooled** (recommended for DE) | `trinity_pooled.sbatch` | One assembly per species for DESeq2 |
| Per-sample | `run_trinity.sbatch` | Explore one sample, or re-run a failed job |

---

## Full workflow (what we have been running)

```
QC reads (fastp / Kraken2 / RiboDetector)
        ↓
Step 1  Pooled Trinity assembly
        ↓
Step 2  BUSCO on pooled Trinity        ← check completeness (high D% is normal)
        ↓
Step 3  CD-HIT deduplication           ← collapse near-duplicate isoforms
        ↓
Step 4  BUSCO on CD-HIT assembly       ← compare D% before vs after
        ↓
Optional (parallel with Step 5):
  4b  Longest isoform + BUSCO          ← build FASTA + third D% benchmark
        ↓
Step 5  RSEM + Bowtie2 per sample      ← gene counts for DESeq2
        ↓
Later   PyDESeq2 → annotation
```

---

## Before you start

### Dependencies

| Tool | How to load / install | Where used |
|------|----------------------|------------|
| Trinity | `module load trinity` (check with `module spider trinity`) | Step 1 |
| BUSCO | `bash scripts/04_busco/setup_busco_env.sh` (one time on login node) | Steps 2 & 4 |
| CD-HIT | `module load CD-HIT/4.8.1-GCC-12.3.0` (check with `module spider cd-hit`) | Step 3 |
| RSEM | `module load RSEM` (check with `module spider RSEM`) | Step 5 |
| Bowtie2 | `module load Bowtie2` (used internally by RSEM) | Step 5 |
| Trinity | `module load trinity` (only for `abundance_estimates_to_matrix.pl`) | Step 5 |

### Input reads (example: MF nutmeg)

Trinity needs paired-end FASTQ files. For MF, `trinity_pooled.sbatch` currently
uses:

- Fruit (MFF): Kraken2-filtered reads in `01_processed/00_00_kraken/`
- Leaf (MFL): fastp-cleaned reads in `01_processed/00_3_MF/`

See [README_kraken.md](../scripts/01_qc/README_kraken.md) if you need to run Kraken2
first.

### Create output folders (first time only)

```bash
mkdir -p /projects/tholl_lab_1/daisy_analysis/01_processed/00_5_BUSCO
mkdir -p /projects/tholl_lab_1/daisy_analysis/01_processed/00_6_cdhit
mkdir -p /projects/tholl_lab_1/daisy_analysis/01_processed/00_7_RSEM
```

---

## Step 1 — Pooled Trinity assembly

**Script:** `scripts/03_assembly/trinity_pooled.sbatch`  
**Runtime:** Can take several days for large datasets (MF uses 7-day wall time).

```bash
# Check no other Trinity job is already writing to the same folder
squeue -u $USER

# Submit
sbatch scripts/03_assembly/trinity_pooled.sbatch
```

**Monitor:**

```bash
squeue -u $USER
tail -f /projects/tholl_lab_1/daisy_analysis/01_processed/00_4_MF_trinity/trinity_pooled_<JOBID>.out
```

**Output (MF example):**

```
01_processed/00_4_MF_trinity/
├── MF_trinity_pooled.Trinity.fasta          ← final assembly (use this downstream)
├── MF_trinity_pooled.Trinity.fasta.gene_trans_map
├── MF_trinity_pooled/                       ← Trinity working files (large)
├── trinity_pooled_<JOBID>.out
└── trinity_pooled_<JOBID>.err
```

**Pooled assembly paths for all species:**

| Species | Code | Output FASTA |
|---------|------|--------------|
| Nutmeg | MF | `01_processed/00_4_MF_trinity/MF_trinity_pooled.Trinity.fasta` |
| Carrot | DC | `01_processed/00_1_DC_trinity/DC_trinity_pooled.Trinity.fasta` |
| Carrot | DG | `01_processed/00_2_DG_trinity/DG_trinity_pooled.Trinity.fasta` |
| Sukman | SK | `01_processed/00_4_Sukman_trinity/SK_trinity_pooled.Trinity.fasta` |
| DC + DG combined | DCDG | `01_processed/00_4_DC_DG_trinity/DCDG_trinity_pooled.Trinity.fasta` |

> **Do not submit two Trinity jobs to the same output folder.** They will
> corrupt each other. Check `squeue` first.

**Per-sample alternative** (one sample at a time):

```bash
sbatch scripts/03_assembly/run_trinity.sbatch 00_3_MF MFF1
```

---

## Step 2 — BUSCO on pooled Trinity

**Script:** `scripts/04_busco/run_busco.sbatch`  
**Purpose:** Measure assembly completeness (Complete / Duplicated / Fragmented / Missing).

**One-time BUSCO setup** (login node, needs internet):

```bash
bash scripts/04_busco/setup_busco_env.sh
```

**Run BUSCO** (pooled Trinity is the default input):

```bash
sbatch scripts/04_busco/run_busco.sbatch MF
sbatch scripts/04_busco/run_busco.sbatch DC
sbatch scripts/04_busco/run_busco.sbatch DG
```

**Lineage picked automatically:**

| Code | Lineage | Why |
|------|---------|-----|
| MF | `embryophyta_odb10` | Nutmeg is a basal angiosperm (land plants) |
| DC, DG, SK, DCDG | `eudicots_odb10` | Carrot is a eudicot (stricter set) |

**Output:**

```
01_processed/00_5_BUSCO/
├── busco_MF/
│   ├── short_summary*.txt       ← read this first
│   ├── run_<lineage>/
│   │   ├── full_table.tsv
│   │   └── missing_busco_list.tsv
│   └── logs/busco.log
├── busco_DC/
└── busco_<JOBID>.out / .err
```

**Reading `short_summary*.txt`:**

```
C:84.5%[S:21.3%,D:63.2%],F:6.1%,M:9.4%,n:1614
```

| Metric | Meaning | What to expect on pooled Trinity |
|--------|---------|----------------------------------|
| **C** (Complete) | Full-length BUSCO gene found | > 70% (ideally > 85%) = good |
| **D** (Duplicated) | Same gene found more than once | **High D% is normal** — Trinity makes many isoforms |
| **F** (Fragmented) | Partial match | Lower is better |
| **M** (Missing) | Gene not found | < 20% is acceptable |

If **C** looks good but **D** is very high, that is expected on raw Trinity.
That is why we run CD-HIT next.

---

## Step 3 — CD-HIT deduplication

**Script:** `scripts/04_busco/run_cdhit.sbatch`  
**Purpose:** Collapse near-duplicate transcripts at 95% identity so isoform
redundancy is reduced before re-running BUSCO.

```bash
sbatch scripts/04_busco/run_cdhit.sbatch MF
```

**What CD-HIT does:**

- Input: pooled Trinity FASTA (hundreds of thousands of transcripts)
- Clusters sequences at ≥ 95% identity
- Keeps one representative per cluster
- Output: smaller FASTA with fewer near-duplicates

**Output:**

```
01_processed/00_6_cdhit/
├── MF_trinity_cdhit95.fasta           ← deduplicated assembly
├── MF_trinity_cdhit95.fasta.clstr     ← which transcripts were merged
├── cdhit_<JOBID>.out
└── cdhit_<JOBID>.err
```

**Check reduction when the job finishes:**

```bash
grep -c '>' /projects/tholl_lab_1/daisy_analysis/01_processed/00_4_MF_trinity/MF_trinity_pooled.Trinity.fasta
grep -c '>' /projects/tholl_lab_1/daisy_analysis/01_processed/00_6_cdhit/MF_trinity_cdhit95.fasta
```

You should see a large drop in transcript count (e.g. 500k → 200k).

---

## Step 4 — BUSCO on CD-HIT assembly

**Script:** `scripts/04_busco/run_busco.sbatch` with second argument `cdhit`

```bash
sbatch scripts/04_busco/run_busco.sbatch MF cdhit
```

This uses the CD-HIT FASTA as input and writes to a **separate** output folder
so you can compare with Step 2:

| Run | Command | Input | Output folder |
|-----|---------|-------|---------------|
| Pooled Trinity | `sbatch ... run_busco.sbatch MF` | `MF_trinity_pooled.Trinity.fasta` | `busco_MF/` |
| CD-HIT | `sbatch ... run_busco.sbatch MF cdhit` | `MF_trinity_cdhit95.fasta` | `busco_MF_cdhit/` |

**What to compare between the two BUSCO runs:**

- **C%** should stay similar (completeness did not get worse)
- **D%** should drop on the CD-HIT run (fewer duplicate BUSCO hits)
- **M%** should stay similar

Alternative (same as `cdhit` argument):

```bash
USE_CDHIT=1 sbatch scripts/04_busco/run_busco.sbatch MF
```

---

## Optional — Longest isoform + BUSCO (three-way comparison)

This step is **optional QC for publication** — it does not replace CD-HIT for
RSEM. One script does both parts: builds the longest-isoform FASTA, then runs
BUSCO on it.

**Script:** `scripts/04_busco/run_busco_longest_isoform.sbatch`

```bash
sbatch scripts/04_busco/run_busco_longest_isoform.sbatch MF
```

**Step 1 (inside the job):** reads pooled Trinity FASTA + `gene_trans_map`,
keeps the longest transcript per gene →
`01_processed/00_6_cdhit/MF_longest_isoform.fasta` (~332k transcripts for MF)

**Step 2 (inside the job):** runs BUSCO on that file →
`busco_MF_longest_isoform/`

To regenerate the FASTA without re-running BUSCO, set `FORCE_REBUILD=1`:

```bash
FORCE_REBUILD=1 sbatch scripts/04_busco/run_busco_longest_isoform.sbatch MF
```

**Three BUSCO runs to compare:**

| Version | Command | Input | Output folder | MF example D% |
|---------|---------|-------|---------------|---------------|
| 1. Pooled Trinity | `run_busco.sbatch MF` | `MF_trinity_pooled.Trinity.fasta` | `busco_MF/` | 78.1% |
| 2. CD-HIT 95% | `run_busco.sbatch MF cdhit` | `MF_trinity_cdhit95.fasta` | `busco_MF_cdhit/` | 65.5% |
| 3. Longest isoform | `run_busco_longest_isoform.sbatch MF` | `MF_longest_isoform.fasta` | `busco_MF_longest_isoform/` | see job log |

**How to interpret:**

- **D drops a lot** (pooled → CD-HIT → longest): isoforms were the main driver
- **D stays elevated** on longest isoform: gene family expansion is real biology
- **C drops below ~88%** on longest isoform: check extraction log — don't use for DE

This BUSCO job can run **in parallel with RSEM** — they use different files.

---

## Step 5 — RSEM + Bowtie2 expression quantification

**Script:** `scripts/04_rsem/run_rsem.sbatch`  
**Purpose:** Align each sample's reads back to the pooled (CD-HIT filtered)
assembly and estimate gene-level read counts for differential expression.

RSEM uses **Bowtie2** internally for alignment, then applies the
expectation-maximization (EM) algorithm so multi-mapping reads are split
probabilistically across transcripts. The gene-level output
(`*.genes.results`) is what you feed into PyDESeq2.

**Prerequisites:** Steps 1–3 must be done (Trinity + CD-HIT). Step 4 (BUSCO on
CD-HIT) is optional but recommended for QC before spending time on RSEM.

```bash
# Default: quantify against CD-HIT assembly (recommended)
sbatch scripts/04_rsem/run_rsem.sbatch MF

# Explicit CD-HIT (same as default)
sbatch scripts/04_rsem/run_rsem.sbatch MF cdhit

# Use raw pooled Trinity instead (not recommended unless testing)
sbatch scripts/04_rsem/run_rsem.sbatch MF pooled
```

**What the script does:**

1. Builds an RSEM + Bowtie2 reference index from your assembly (one time)
2. Runs `rsem-calculate-expression` on each sample (skips samples already done)
3. Prints alignment rate per sample (target **> 70%**)
4. Combines all `.genes.results` files into one count matrix

**Input files (MF):**

| Item | Path |
|------|------|
| Assembly (default) | `01_processed/00_6_cdhit/MF_trinity_cdhit95.fasta` |
| Gene map | `01_processed/00_4_MF_trinity/MF_trinity_pooled.Trinity.fasta.gene_trans_map` |
| Fruit reads (MFF) | `01_processed/00_00_kraken/MFF*_filtered.fq.gz` |
| Leaf reads (MFL) | `01_processed/00_3_MF/MFL*_clean.fq.gz` |

**Output:**

```
01_processed/00_7_RSEM/
├── MF_rsem_ref/                         ← Bowtie2 index (built once)
├── MFF1/
│   ├── MFF1.genes.results               ← USE FOR DESEQ2
│   ├── MFF1.isoforms.results
│   └── MFF1.transcript.bam
├── MFF2/  MFF3/  MFL1/  MFL2/  MFL3/    ← same structure per sample
├── RSEM.gene.counts.matrix              ← combined matrix for PyDESeq2
├── RSEM.isoform.counts.matrix
└── rsem_<JOBID>.out / .err
```

**First few lines of `RSEM.gene.counts.matrix`:**

```
gene_id    MFF1    MFF2    MFF3    MFL1    MFL2    MFL3
TRINITY_DN0_c0_g1_i1    42.50    38.20    45.10    12.30    11.80    13.50
TRINITY_DN0_c0_g2_i1     0.00     0.00     0.00     5.20     4.90     5.10
```

**Strandedness:** MF libraries are **non-directional** (`--strandedness none`).
If you change library prep for a new project, update that flag in
`run_rsem.sbatch`.

**After RSEM finishes:**

```bash
# Check the count matrix exists
ls -lh /projects/tholl_lab_1/daisy_analysis/01_processed/00_7_RSEM/RSEM.gene.counts.matrix

# Next: differential expression
# See scripts/05_pydeseq2/
```

---

## Quick command cheat sheet (MF)

Run these in order on ARC:

```bash
cd /projects/tholl_lab_1/daisy_analysis/05_rnaseq-code

# 1. Pooled Trinity (days)
sbatch scripts/03_assembly/trinity_pooled.sbatch

# 2. BUSCO on pooled Trinity (~hours)
sbatch scripts/04_busco/run_busco.sbatch MF

# 3. CD-HIT (~hours)
sbatch scripts/04_busco/run_cdhit.sbatch MF

# 4. BUSCO on CD-HIT assembly (~hours)
sbatch scripts/04_busco/run_busco.sbatch MF cdhit

# Optional 4b: longest isoform + BUSCO (can run parallel with step 5)
sbatch scripts/04_busco/run_busco_longest_isoform.sbatch MF

# 5. RSEM + Bowtie2 per sample (~hours per run, all 6 samples in one job)
sbatch scripts/04_rsem/run_rsem.sbatch MF
```

**All species — CD-HIT + BUSCO cdhit:**

```bash
for CODE in MF DC DG SK DCDG; do
    sbatch scripts/04_busco/run_cdhit.sbatch "${CODE}"
done

# After CD-HIT jobs finish:
for CODE in MF DC DG SK DCDG; do
    sbatch scripts/04_busco/run_busco.sbatch "${CODE}" cdhit
done
```

---

## Troubleshooting

### Trinity job did not finish / timed out

- Check the `.err` log in the Trinity output folder.
- Do **not** delete the `MF_trinity_pooled/` working directory — Trinity can
  resume from checkpoint if you resubmit (for per-sample runs;
  `run_trinity.sbatch` handles this automatically).
- For pooled runs, you may need more time (`--time`) or more memory (`--mem`).

### BUSCO says assembly not found

**Pooled run (default):**

```bash
ls -lh /projects/tholl_lab_1/daisy_analysis/01_processed/00_4_MF_trinity/MF_trinity_pooled.Trinity.fasta
```

Trinity must finish first.

**CD-HIT run (`cdhit` argument):**

```bash
ls -lh /projects/tholl_lab_1/daisy_analysis/01_processed/00_6_cdhit/MF_trinity_cdhit95.fasta
```

Run `run_cdhit.sbatch` first.

### `module load CD-HIT/4.8.1-GCC-12.3.0` fails

```bash
module spider cd-hit
```

Update the module line in `run_cdhit.sbatch` to match your cluster.

### `busco` not found

```bash
bash scripts/04_busco/setup_busco_env.sh
```

Run on a login node with internet, then resubmit.

### CD-HIT runs out of memory

Edit `run_cdhit.sbatch`:

- Raise `#SBATCH --mem=64G` → `128G`
- Raise `CDHIT_MEM_MB=60000` → `120000`

### RSEM says assembly or gene map not found

**CD-HIT assembly (default):**

```bash
ls -lh /projects/tholl_lab_1/daisy_analysis/01_processed/00_6_cdhit/MF_trinity_cdhit95.fasta
```

Run `run_cdhit.sbatch` first.

**Gene map:**

```bash
ls -lh /projects/tholl_lab_1/daisy_analysis/01_processed/00_4_MF_trinity/MF_trinity_pooled.Trinity.fasta.gene_trans_map
```

Trinity pooled assembly must finish first. The gene map always comes from the
**original pooled Trinity output**, even when you quantify against the CD-HIT
FASTA.

### RSEM alignment rate is low (< 60%)

- Confirm you used the same read files as Trinity (MFF = Kraken-filtered,
  MFL = fastp-cleaned for MF).
- Check `--strandedness` matches your library prep.
- Re-run BUSCO on the CD-HIT assembly — if C% dropped badly, the assembly
  may be the problem, not RSEM.

### Longest isoform BUSCO fails at Step 1

Pooled Trinity must exist first:

```bash
ls -lh /projects/tholl_lab_1/daisy_analysis/01_processed/00_4_MF_trinity/MF_trinity_pooled.Trinity.fasta
```

### `abundance_estimates_to_matrix.pl` not found

```bash
module load trinity
```

Then re-run the matrix step manually (the script prints the exact command), or
resubmit the job after loading Trinity in the sbatch module block.

---

## What comes after this workflow

Once you are happy with assembly quality (BUSCO C% and acceptable D% after
CD-HIT) and RSEM alignment rates look good:

1. **PyDESeq2** — differential expression with `RSEM.gene.counts.matrix`
2. **TransDecoder** — predict ORFs from transcripts
3. **BLAST / annotation** — assign gene functions

For Kraken2 decontamination details, see
[README_kraken.md](../scripts/01_qc/README_kraken.md).

For full pipeline context, see the main [README.md](../README.md).
