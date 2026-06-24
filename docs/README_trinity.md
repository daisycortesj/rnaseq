# Trinity Assembly + Quality Check Workflow

**Scripts:** `scripts/03_assembly/trinity_pooled.sbatch`, `scripts/03_assembly/run_trinity.sbatch`  
**Related scripts:** `scripts/03_assembly/run_transdecoder.sbatch`, `scripts/07_domains/run_hmmer_genefinder.sbatch`, `scripts/04_busco/run_busco.sbatch`, `scripts/04_busco/run_cdhit.sbatch`, `scripts/04_busco/run_busco_longest_isoform.sbatch`, `scripts/04_rsem/run_rsem.sbatch`, `scripts/04_counting/run_count_matrix.sbatch`  
**Where it fits:** fastp → (optional Kraken2) → **Trinity (here)** → BUSCO → CD-HIT → BUSCO again → **RSEM** → **count matrix** → TransDecoder → **HMMER CYP/OMT** → PyDESeq2 → BLAST

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
Step 5  RSEM + Bowtie2 per sample      ← gene counts (decimals)
        ↓
Step 6  Build count matrix + metadata  ← integer counts + sample groups for PyDESeq2
        ↓
Step 7  TransDecoder ORF prediction    ← nucleotide → protein (slow; run once)
        ↓
Step 8  PyDESeq2 differential expression (Fruit vs Leaf)
        ↓
Step 9  HMMER CYP/OMT domain scan     ← find family candidates (fast; minutes)
        ↓
Step 10 BLAST + combine (Script 03–04) ← confirm and rank candidates
```

---

## Before you start

### Dependencies

| Tool | How to load / install | Where used |
|------|----------------------|------------|
| Trinity | `module load trinity` (check with `module spider trinity`) | Step 1 |
| BUSCO | `bash scripts/04_busco/setup_busco_env.sh` (one time on login node) | Steps 2 & 4 |
| CD-HIT | `module load CD-HIT/4.8.1-GCC-12.3.0` (check with `module spider cd-hit`) | Step 3 |
| RSEM + Bowtie2 | `conda activate rnaseq` (from `environment.yml` — no ARC module) | Step 5 |
| TransDecoder | `conda install -y -c bioconda transdecoder` in `rnaseq` env | Step 7 |
| HMMER | same `rnaseq` env (`conda install -c bioconda hmmer`) | Step 9 |
| Trinity utility | same `rnaseq` env (`abundance_estimates_to_matrix.pl`) | Step 5 |
| Python + pandas | same `rnaseq` env | Step 6 |

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

# Next: Step 6 — build PyDESeq2-ready files (see below)
```

---

## Step 6 — Build count matrix + sample metadata

**Script:** `scripts/04_counting/run_count_matrix.sbatch`  
**Purpose:** Convert RSEM output into the two files PyDESeq2 needs:

1. **`gene_count_matrix.tsv`** — genes × samples, **integer** counts  
2. **`sample_metadata.tsv`** — which sample belongs to which group (Fruit vs Leaf for MF)

RSEM writes **decimal** expected counts (e.g. `42.50`). This step rounds them to
whole numbers and builds the metadata table from your sample names.

**Prerequisite:** Step 5 (RSEM) must be done. You need
`01_processed/00_7_RSEM/RSEM.gene.counts.matrix`.

```bash
# Trinity pipeline — always use --type rsem
sbatch scripts/04_counting/run_count_matrix.sbatch MF --type rsem
```

**What the script does:**

1. Reads `RSEM.gene.counts.matrix` from `01_processed/00_7_RSEM/`
2. Rounds expected counts to integers
3. Builds `sample_metadata.tsv` using sample names + `config.sh`  
   (for MF: `MFF=F` means Fruit, `MFL=L` means Leaf)
4. Saves everything to `03_count_tables/00_5_MF_trinity/` (where PyDESeq2 looks for MF Trinity)

**Output:**

```
03_count_tables/00_5_MF_trinity/
├── gene_count_matrix.tsv    ← PyDESeq2 count input (integers)
├── sample_metadata.tsv      ← PyDESeq2 design input
└── count_summary.txt        ← quick stats to sanity-check
```

**First few lines of `gene_count_matrix.tsv`:**

```
gene_id    MFF1    MFF2    MFF3    MFL1    MFL2    MFL3
TRINITY_DN0_c0_g1_i1    43    38    45    12    12    14
TRINITY_DN0_c0_g2_i1     0     0     0     5     5     5
```

**First few lines of `sample_metadata.tsv`:**

```
sample  group  condition  replicate
MFF1    MFF1   F          1
MFF2    MFF2   F          2
MFL1    MFL1   L          1
```

The important column for PyDESeq2 is **`condition`** (`F` = Fruit, `L` = Leaf).
The script reads `MFF=F MFL=L` from `config.sh` automatically.

**After this step finishes:**

```bash
# Check files exist
ls -lh /projects/tholl_lab_1/daisy_analysis/03_count_tables/00_5_MF_trinity/

# Next: Step 7 — TransDecoder (after RSEM is done)
sbatch scripts/03_assembly/run_transdecoder.sbatch MF

# Step 8 — PyDESeq2 (can run while TransDecoder is still going)
CONTRAST_A=F CONTRAST_B=L sbatch scripts/05_pydeseq2/run_step1_analysis.sbatch MF
```

---

## Step 7 — TransDecoder ORF prediction (CYP/OMT pipeline)

**Script:** `scripts/03_assembly/run_transdecoder.sbatch`  
**Purpose:** Trinity outputs **nucleotide** transcripts (DNA letters: A, T, G, C).
BLAST and HMMER need **protein** sequences (amino acids). TransDecoder finds
the coding region (ORF) in each transcript and writes a translated FASTA.

**When to run:** **After Step 5 (RSEM) and Step 6 (count matrix).** The sbatch
script checks that `RSEM.gene.counts.matrix` exists before it starts.

> **Important:** TransDecoder does **not** read RSEM output files. It reads the
> CD-HIT FASTA from Step 3. We run it after RSEM so you finish counting reads
> before starting this slow annotation step. Step 8 (PyDESeq2) can run **while
> TransDecoder is still going** — they do not block each other.

This is the **slowest** step in the CYP450 + OMT annotation pipeline — run it
**once** and reuse the `.pep` file. The script skips automatically if output
already exists (pass `--force` to regenerate).

**Prerequisites:**

| Step | What must be done |
|------|-------------------|
| 1 + 3 | Trinity + CD-HIT (provides the input FASTA) |
| 5 | RSEM finished (`RSEM.gene.counts.matrix` exists) |
| 6 | Count matrix built (recommended before submitting) |

Use the **CD-HIT assembly** (default) so protein IDs match RSEM / PyDESeq2 counts.

**One-time install** (optional — the sbatch script auto-installs if missing):

```bash
conda activate rnaseq
conda install -y -c bioconda -c conda-forge transdecoder
which TransDecoder.LongOrfs TransDecoder.Predict
```

If you already installed TransDecoder, the job skips this step automatically.

**Submit on ARC (MF nutmeg):**

```bash
cd /projects/tholl_lab_1/daisy_analysis/05_rnaseq-code

# Default: CD-HIT assembly (same IDs as RSEM)
sbatch scripts/03_assembly/run_transdecoder.sbatch MF

# Re-run from scratch (deletes old .pep first)
sbatch scripts/03_assembly/run_transdecoder.sbatch MF cdhit --force
```

**Monitor:**

```bash
squeue -u $USER
tail -f /projects/tholl_lab_1/daisy_analysis/01_processed/00_6_cdhit/transdecoder_<JOBID>.out
```

**Input / output (MF example):**

| Item | Path |
|------|------|
| Input (default) | `01_processed/00_6_cdhit/MF_trinity_cdhit95.fasta` |
| **Main output** | `01_processed/00_6_cdhit/MF_trinity_cdhit95.fasta.transdecoder.pep` |
| Also created | `.transdecoder.cds`, `.transdecoder.gff3`, `.transdecoder.bed` |

**Protein header format** (used by later CYP/OMT scripts):

```
>TRINITY_DN123_c0_g1_i1.p1 type:complete len:412 ...
```

The transcript ID is everything before the trailing `.p1`. Later scripts remove
`.p1`, `.p2`, etc. to recover the Trinity ID for joining to DESeq2 counts.

**First few lines of `.transdecoder.pep`:**

```
>TRINITY_DN0_c0_g1_i1.p1 type:complete len:312 ...
MKVLFLLLIA...
>TRINITY_DN0_c0_g2_i1.p1 type:5prime_partial len:89 ...
MSTKELV...
```

---

## Step 8 — PyDESeq2 differential expression

After the count matrix (Step 6) is ready, run PyDESeq2. You can submit this
**before or while** TransDecoder (Step 7) is still running.

```bash
CONTRAST_A=F CONTRAST_B=L sbatch scripts/05_pydeseq2/run_step1_analysis.sbatch MF
```

See the main [README.md](../README.md) for the full BLAST/combine steps that
use the HMMER candidate lists from Step 9.

---

## Step 9 — HMMER CYP/OMT gene finder (Script 02)

**Script:** `scripts/07_domains/run_hmmer_genefinder.sbatch`

**Purpose:** Find cytochrome P450 (CYP) and O-methyltransferase (OMT) candidate
transcripts in your Trinity assembly by searching TransDecoder proteins for
conserved Pfam domains.

### Why this step matters for your nutmeg work

Your project asks: *which CYP and OMT genes are expressed in nutmeg fruit vs
leaf, and how do they change?* Trinity gives you thousands of transcript IDs,
but most are unrelated to terpene biosynthesis. HMMER is your **first sensitive
filter** — it finds transcripts whose predicted proteins contain the structural
domains that define each family:

| Family | Pfam domain | What it means biologically |
|--------|-------------|----------------------------|
| CYP450 | PF00067 (Cytochrome_P450) | Enzymes that oxidize/modify terpenoid backbones |
| OMT | PF00891, PF08100, PF01596 | Enzymes that add methyl groups to terpenoids |

HMMER profiles are trained on hundreds of known family members, so they catch
divergent nutmeg genes that a simple keyword search would miss. The output
transcript IDs match your RSEM/PyDESeq2 count matrix (same Trinity IDs).

### How this relates to other HMMER scripts in the repo

| Script | When to use | Input |
|--------|-------------|-------|
| **`run_hmmer_genefinder.sbatch`** (this step) | Trinity de novo (MF nutmeg) | `.transdecoder.pep` |
| `run_hmmer.sbatch` | Carrot reference genome (DC/DG) | `all_genes_protein.fasta` — scans **all** Pfam domains |
| `gene_hmmer_scan.py` | Building carrot CYP/OMT master lists | One Pfam at a time, outputs CSV |

For MF nutmeg, use **`run_hmmer_genefinder.sbatch`** — it runs four targeted
`hmmsearch` commands with Pfam's curated threshold (`--cut_ga`) and writes
simple ID lists you can intersect with DESeq2 results.

### Inputs required

Everything below must be in place before you submit `run_hmmer_genefinder.sbatch`.

#### 1. Upstream pipeline (MF nutmeg)

| Step | Script | What it produces |
|------|--------|------------------|
| Trinity assembly | `trinity_pooled.sbatch` | Pooled MF transcriptome |
| CD-HIT dedup | `run_cdhit.sbatch MF` | `MF_trinity_cdhit95.fasta` |
| RSEM quantification | `run_rsem.sbatch MF` | Per-sample counts (same transcript IDs) |
| **TransDecoder** | **`run_transdecoder.sbatch MF`** | **`.transdecoder.pep` ← main input for this step** |

#### 2. Required input file

| Item | MF path |
|------|---------|
| **Protein FASTA** | `01_processed/00_6_cdhit/MF_trinity_cdhit95.fasta.transdecoder.pep` |

Check it exists:

```bash
ls -lh /projects/tholl_lab_1/daisy_analysis/01_processed/00_6_cdhit/MF_trinity_cdhit95.fasta.transdecoder.pep
grep -c '^>' /projects/tholl_lab_1/daisy_analysis/01_processed/00_6_cdhit/MF_trinity_cdhit95.fasta.transdecoder.pep
```

Example header (TransDecoder format):

```
>TRINITY_DN0_c0_g1_i1.p1 type:complete len:312 ...
MKVLFLLLIA...
```

#### 3. Software (rnaseq conda env)

| Tool | Check | Install if missing |
|------|-------|-------------------|
| `hmmsearch` | `hmmsearch -h` | `conda install -c bioconda hmmer` |
| `hmmpress` | `hmmpress -h` | (same package as hmmsearch) |

```bash
conda activate rnaseq
hmmsearch -h    # should print HMMER version, not "command not found"
```

#### 4. Pfam HMM profiles (auto-downloaded on first run)

Saved to `07_NRdatabase/hmmerdb/pfam_profiles/` if not already present:

| File | Family |
|------|--------|
| `PF00067.hmm` | Cytochrome P450 (CYP) |
| `PF00891.hmm` | Methyltransf_2 (OMT) |
| `PF08100.hmm` | OMT dimerisation domain |
| `PF01596.hmm` | CCoA-type OMT |

> **Note on OMT domains:** PF00067 (P450) is well verified. Confirm the three
> OMT accessions on [InterPro](https://www.ebi.ac.uk/interpro/) for your specific
> OMT subtype before you treat them as definitive.

#### 5. Command-line arguments

| Argument | Required? | Meaning |
|----------|-----------|---------|
| `MF` (species code) | **Yes** | Must match a code in `scripts/config.sh` (MF = nutmeg) |
| `--force` | No | Re-run even if output files already exist |

### Submit on ARC (MF nutmeg)

```bash
cd /projects/tholl_lab_1/daisy_analysis/05_rnaseq-code

# Default — skips if outputs already exist
sbatch scripts/07_domains/run_hmmer_genefinder.sbatch MF

# Re-run from scratch
sbatch scripts/07_domains/run_hmmer_genefinder.sbatch MF --force
```

**Runtime:** Usually a few minutes (much faster than full Pfam `hmmscan`).

**Monitor:**

```bash
squeue -u $USER
tail -f /projects/tholl_lab_1/daisy_analysis/06_analysis/hmmer_genefinder_<JOBID>.out
```

### Output (MF example)

| Item | Path |
|------|------|
| **CYP candidates** | `06_analysis/hmmer_genefinder_MF/cyp450_hmmer_ids.txt` |
| **OMT candidates** | `06_analysis/hmmer_genefinder_MF/omt_hmmer_ids.txt` |
| Raw domain tables | `06_analysis/hmmer_genefinder_MF/*.domtbl` (keep for e-value filtering) |

**First few lines of `cyp450_hmmer_ids.txt`:**

```
TRINITY_DN123_c0_g1_i1
TRINITY_DN456_c0_g2_i3
TRINITY_DN789_c1_g1_i1
```

One Trinity transcript ID per line — no `.p1` suffix (stripped automatically).

### What the script does internally

1. Checks `hmmsearch` is installed and prints the HMMER version.
2. Verifies the `.pep` file and all four Pfam HMM files exist (downloads if needed).
3. Runs four `hmmsearch --cut_ga` searches against the protein FASTA.
4. Parses each `.domtbl` file — column 1 is the protein ID.
5. Strips trailing `.p1`, `.p2`, etc. to recover transcript IDs.
6. Merges the three OMT domain results into one OMT candidate set.
7. Writes `cyp450_hmmer_ids.txt` and `omt_hmmer_ids.txt`.

### Troubleshooting

**`TransDecoder protein file not found`**

Run Step 7 first:

```bash
sbatch scripts/03_assembly/run_transdecoder.sbatch MF
```

**`hmmsearch not found`**

```bash
conda activate rnaseq
conda install -c bioconda hmmer
```

**Could not download Pfam HMM**

Download manually from InterPro and save to
`07_NRdatabase/hmmerdb/pfam_profiles/PF00067.hmm` (repeat for PF00891, PF08100, PF01596).

**Zero candidates found**

- Check the `.pep` file has sequences: `grep -c '^>' ...transdecoder.pep`
- Inspect raw `.domtbl` files — maybe `--cut_ga` is too strict for very divergent nutmeg genes.
- Confirm TransDecoder ran on the CD-HIT assembly (same IDs as RSEM).

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

# 5. RSEM + Bowtie2 per sample (~hours — all 6 samples in one job)
sbatch scripts/04_rsem/run_rsem.sbatch MF

# 6. Build count matrix + metadata for PyDESeq2 (~minutes)
sbatch scripts/04_counting/run_count_matrix.sbatch MF --type rsem

# 7. TransDecoder ORF → protein (slow — run once; after RSEM)
sbatch scripts/03_assembly/run_transdecoder.sbatch MF

# 8. PyDESeq2 differential expression (Fruit vs Leaf)
#    Can submit while step 7 is still running
CONTRAST_A=F CONTRAST_B=L sbatch scripts/05_pydeseq2/run_step1_analysis.sbatch MF

# 9. HMMER CYP/OMT gene finder (~minutes — after TransDecoder finishes)
sbatch scripts/07_domains/run_hmmer_genefinder.sbatch MF
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

### RSEM fails with Perl mismatch (`Cwd.c: loadable library and perl binaries are mismatched`)

This happened when the script used `module load RSEM` + `module load Bowtie2` — ARC
has no RSEM module, and unpinned Bowtie2 swaps GCC/Perl. The script now uses the
`rnaseq` conda env only (RSEM + Bowtie2 installed together).

Verify on the login node:

```bash
conda activate rnaseq
which rsem-calculate-expression bowtie2
rsem-calculate-expression --version
```

If missing, install once (login node, needs internet):

```bash
conda activate rnaseq
conda install -y -c bioconda -c conda-forge rsem bowtie2 trinity
```

If a failed job left a partial index, remove it before resubmitting:

```bash
rm -rf /projects/tholl_lab_1/daisy_analysis/01_processed/00_7_RSEM/MF_rsem_ref
```

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

### Count matrix step says "RSEM count matrix not found"

RSEM must finish first:

```bash
ls -lh /projects/tholl_lab_1/daisy_analysis/01_processed/00_7_RSEM/RSEM.gene.counts.matrix
```

If missing, re-run:

```bash
sbatch scripts/04_rsem/run_rsem.sbatch MF
```

### PyDESeq2 says "Count matrix not found"

You need Step 6 after RSEM:

```bash
sbatch scripts/04_counting/run_count_matrix.sbatch MF --type rsem
```

Check output:

```bash
ls -lh /projects/tholl_lab_1/daisy_analysis/03_count_tables/00_5_MF_trinity/gene_count_matrix.tsv
ls -lh /projects/tholl_lab_1/daisy_analysis/03_count_tables/00_5_MF_trinity/sample_metadata.tsv
```

### TransDecoder not found on PATH

Install in the `rnaseq` conda env (login node):

```bash
conda activate rnaseq
conda install -y -c bioconda transdecoder
which TransDecoder.LongOrfs TransDecoder.Predict
```

Then resubmit:

```bash
sbatch scripts/03_assembly/run_transdecoder.sbatch MF
```

### TransDecoder says RSEM count matrix not found

You must finish Steps 5 and 6 first:

```bash
ls -lh /projects/tholl_lab_1/daisy_analysis/01_processed/00_7_RSEM/RSEM.gene.counts.matrix
sbatch scripts/04_rsem/run_rsem.sbatch MF
sbatch scripts/04_counting/run_count_matrix.sbatch MF --type rsem
```

### TransDecoder says assembly not found

CD-HIT must finish first (default input):

```bash
ls -lh /projects/tholl_lab_1/daisy_analysis/01_processed/00_6_cdhit/MF_trinity_cdhit95.fasta
sbatch scripts/04_busco/run_cdhit.sbatch MF
```

### TransDecoder skipped but you want to re-run

The script skips if `.transdecoder.pep` already exists. Force a fresh run:

```bash
sbatch scripts/03_assembly/run_transdecoder.sbatch MF cdhit --force
```

---

## What comes after this workflow

Once you are happy with assembly quality (BUSCO C% and acceptable D% after
CD-HIT) and RSEM alignment rates look good:

1. **Count matrix + metadata** — Step 6 above
2. **TransDecoder** — Step 7 above (proteins for BLAST/HMMER)
3. **PyDESeq2** — Step 8 above (can run while TransDecoder is going)
4. **HMMER CYP/OMT scan** — Step 9 above (targeted domain search on Trinity proteins)
5. **BLAST + combine** — Scripts 03–04 in your CYP/OMT pipeline (see main [README.md](../README.md))

For Kraken2 decontamination details, see
[README_kraken.md](../scripts/01_qc/README_kraken.md).

For full pipeline context, see the main [README.md](../README.md).
