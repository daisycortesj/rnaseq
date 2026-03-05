# Count Matrix Pipeline

How raw read counts go from BAM files to a DESeq2-ready gene count matrix.

There are two supported paths. Both produce the same final output files that
PyDESeq2 needs.

---

## Pipeline Overview

```
                      ┌─────────────────────────────────────────┐
                      │        STAR Alignment (already done)    │
                      │  *_Aligned.sortedByCoord.out.bam        │
                      │  *_ReadsPerGene.out.tab (if --quantMode)│
                      └──────────┬──────────────┬───────────────┘
                                 │              │
                    PATH A (STAR)│              │PATH B (featureCounts)
                                 │              │
                                 │    ┌─────────▼──────────┐
                                 │    │  samtools sort      │
                                 │    │  Re-sort BAMs       │
                                 │    └─────────┬──────────┘
                                 │              │
                                 │    ┌─────────▼──────────┐
                                 │    │  featureCounts      │
                                 │    │  Count reads/gene   │
                                 │    └─────────┬──────────┘
                                 │              │
                      ┌──────────▼──────────────▼───────────────┐
                      │       build_count_matrix.py             │
                      │  Parses either source into clean matrix │
                      └──────────────────┬──────────────────────┘
                                         │
                              ┌───────────▼────────────┐
                              │  gene_count_matrix.tsv  │
                              │  sample_metadata.tsv    │
                              │  count_summary.txt      │
                              └────────────────────────┘
                                         │
                              ┌───────────▼────────────┐
                              │  PyDESeq2 analysis     │
                              └────────────────────────┘
```

---

## PATH A — STAR GeneCounts (your original pipeline)

STAR counts reads per gene **during** alignment when you pass
`--quantMode GeneCounts`. It outputs one `ReadsPerGene.out.tab` file per
sample. `build_count_matrix.py` merges those individual files into a single
matrix.

**When to use:** Any species where STAR was run with `--quantMode GeneCounts`
(DC and DG with carrot genome). Also the only option for MF (nutmeg) since
there is no GTF for featureCounts.

### Command

```bash
sbatch scripts/build_count_matrix.sbatch DC --type star
```

That single command is all you need — STAR already produced the count files
during alignment.

---

## PATH B — featureCounts (previous student's workflow)

This follows the traditional approach: align first, then count separately.
featureCounts (from the Subread package) reads sorted BAM files and a GTF
annotation, counts how many reads overlap each exon, and summarizes by
`gene_id`. It produces a single combined output file with all samples.

**When to use:** DC and DG (requires a GTF). Cannot be used for MF (no GTF).

### Step 1 — samtools sort

Re-sorts STAR BAMs with samtools to ensure a standard sort header.

| | |
|---|---|
| **Script** | `scripts/samtools_sort.sbatch` |
| **Tool** | samtools sort + samtools index |
| **Input** | `02_mapped/{species}/*_Aligned.sortedByCoord.out.bam` |
| **Output** | `02_mapped/{species}/*_Aligned.sortedBySamtools.bam` (+ `.bai` index) |
| **SLURM** | 1 node, 8 CPUs, 16 GB, 3 hours |

```bash
sbatch scripts/samtools_sort.sbatch DC
```

**What happens:** For each BAM file in the species directory, samtools
re-sorts it by genomic coordinate and creates an index (`.bai`). If a
sorted file already exists, it skips that sample.

**Why this step?** STAR sorts during alignment, but samtools sort writes a
standard `@HD SO:coordinate` header that some tools are strict about. This
mirrors the previous student's workflow.

### Step 2 — featureCounts

Counts how many read pairs map to each gene across all samples at once.

| | |
|---|---|
| **Script** | `scripts/featurecounts.sbatch` |
| **Tool** | featureCounts (Subread package) |
| **Input** | `02_mapped/{species}/*_Aligned.sortedBySamtools.bam` + GTF |
| **Output** | `03_count_tables/{species}/featurecounts.txt` (+ `.summary`) |
| **GTF** | `04_reference/dc_genomic.gtf` (carrot genome, used for DC and DG) |
| **SLURM** | 1 node, 8 CPUs, 16 GB, 3 hours |

```bash
sbatch scripts/featurecounts.sbatch DC
```

**What happens:** featureCounts reads all sorted BAMs plus the GTF annotation.
For each gene, it counts how many read pairs overlap its exons. When a read
overlaps multiple genes, it is assigned to the one with the largest overlap
(`--largestOverlap`). Output is a single tab-delimited file with annotation
columns (Geneid, Chr, Start, End, Strand, Length) followed by one count
column per BAM.

**Key flags used:**

| Flag | Purpose |
|------|---------|
| `-a <GTF>` | Annotation file with gene/exon coordinates |
| `-t exon` | Count at exon level (standard for RNA-seq) |
| `-g gene_id` | Summarize exon counts by gene_id attribute |
| `--largestOverlap` | Ambiguous reads go to the gene with most overlap |
| `-p --countReadPairs` | Paired-end mode: count fragments, not individual reads |
| `-T <threads>` | Parallelize across CPUs |

### Automatic — Build count matrix + sample metadata

After featureCounts finishes, the same sbatch job automatically calls
`build_count_matrix.py` to produce the PyDESeq2-ready files. No extra
step needed.

**What happens automatically:** The Python script strips featureCounts
annotation columns, converts full BAM paths to clean sample names (e.g.,
`/path/DC1L1_Aligned.sortedBySamtools.bam` becomes `DC1L1`), extracts
biological metadata from sample names (species, condition, replicate), and
writes `gene_count_matrix.tsv` + `sample_metadata.tsv` + `count_summary.txt`.

---

## Output Files

All three files land in `03_count_tables/{species}/` (e.g., `03_count_tables/00_1_DC/`).

### `gene_count_matrix.tsv`

The main input for PyDESeq2. Tab-delimited, genes as rows, samples as columns.

```
            DC1L1   DC1L2   DC1R1   DC1R2   ...
gene_001      150     130     200     180
gene_002       42      38      10      12
...
```

### `sample_metadata.tsv`

Tells PyDESeq2 which sample belongs to which experimental group.

```
sample  group  group_number  condition  replicate  treatment  full_condition
DC1L1   DC     1             L          1          DC1        DC1_L
DC1R1   DC     1             R          1          DC1        DC1_R
...
```

### `count_summary.txt`

Human-readable stats: total samples, total genes, per-sample read totals,
and gene filtering thresholds (how many genes have zero/low counts).

---

## Quick Reference — Full Command Sequences

### PATH A (STAR — one command)

```bash
sbatch scripts/build_count_matrix.sbatch DC --type star
```

### PATH B (featureCounts — two commands, run sequentially)

```bash
# Wait for each job to finish before submitting the next
sbatch scripts/samtools_sort.sbatch DC
sbatch scripts/featurecounts.sbatch DC
```

The second command runs featureCounts **and** builds the count matrix +
sample metadata automatically — no third step needed.

### After either path — run PyDESeq2

```bash
sbatch scripts/run_pydeseq2_step1_analysis.sbatch DC
```

---

## Species Support

| Species | Code | PATH A (STAR) | PATH B (featureCounts) | Reason |
|---------|------|:---:|:---:|--------|
| D. carota ssp. maximus | DC | Yes | Yes | Has GTF (`dc_genomic.gtf`) |
| Daucus glaber | DG | Yes | Yes | Aligned to carrot genome (same GTF) |
| Myristica fragrans | MF | Yes | No | No GTF available |

---

## When Do You Need `build_count_matrix.sbatch`?

| Path | Do you need `build_count_matrix.sbatch`? | Why? |
|------|:---:|--------|
| PATH A (STAR) | **Yes** | Merges per-sample STAR files into one matrix + generates metadata |
| PATH B (featureCounts) | **No** | `featurecounts.sbatch` does this automatically at the end |

For PATH B, `featurecounts.sbatch` runs featureCounts and then
immediately calls `build_count_matrix.py` to produce the clean matrix +
metadata. If that second step fails for some reason, you can retry it
manually with:

```bash
sbatch scripts/build_count_matrix.sbatch DC --type featurecounts
```

---

## Differences Between the Two Paths

| Aspect | PATH A (STAR GeneCounts) | PATH B (featureCounts) |
|--------|------------------------|----------------------|
| When counting happens | During alignment | After alignment |
| Number of steps | 1 | 3 |
| Overlap handling | Read must map unambiguously to one gene | Configurable (`--largestOverlap`) |
| Counting unit | Whole gene body | Exon-level, summarized to gene |
| Output format | One file per sample | One combined file |
| GTF required at count time | No (baked into STAR index) | Yes (passed to featureCounts) |
| Paired-end awareness | Counts individual reads | Counts read pairs (`--countReadPairs`) |
