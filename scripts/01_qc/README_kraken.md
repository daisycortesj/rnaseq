# Kraken2 Pre-Assembly Decontamination

**Script:** `scripts/01_qc/kraken2_prefilter.sbatch`
**Where it fits:** fastp → RiboDetector → **Kraken2 (here)** → Trinity

---

## Why we do this

MultiQC flagged a GC content warning for all 6 samples. The fruit samples
(MFF) averaged ~50% GC while the leaf samples (MFL) averaged ~41% — a 9-point
gap that strongly suggests microbial contamination in the fruit tissue. Plants
run 40–45% GC; fungi and bacteria often run higher.

Kraken2 classifies every read by comparing short sequences (k-mers) against a
reference database. Reads that **do not match** anything are likely your nutmeg
reads and are saved as `*_filtered.fq.gz` for Trinity. Reads that **do match**
are likely contamination and are saved separately for your records.

---

## Step 0 — Check if the database already exists

Before submitting the script, run this on a login node:

```bash
ls /projects/tholl_lab_1/daisy_analysis/07_NRdatabase/00_1_kraken2_pluspf/
```

You need to see all three of these files:
```
hash.k2d
opts.k2d
taxo.k2d
```

If those files are there, skip to **Step 2**.

Also worth checking if ARC has a shared copy (saves you 77 GB of downloading):
```bash
ls /projects/arcsingularity/kraken2/ 2>/dev/null
ls /apps/ 2>/dev/null | grep -i kraken
```
Or email arc@vt.edu and ask: *"Do you have a shared Kraken2 PlusPF database?"*

---

## Step 1 — Download the database (only if it does not exist)

Run these commands on a **login node** (not in a job). This takes a while due
to the 77 GB download size.

```bash
# Create the folder
mkdir -p /projects/tholl_lab_1/daisy_analysis/07_NRdatabase/00_1_kraken2_pluspf

# Move into it
cd /projects/tholl_lab_1/daisy_analysis/07_NRdatabase/00_1_kraken2_pluspf

# Download the PlusPF database (check for a newer date at the URL below)
# Latest releases: https://benlangmead.github.io/aws-indexes/k2
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20240904.tar.gz

# Extract it (~100 GB after extraction)
tar -xzvf k2_pluspf_20240904.tar.gz

# Delete the compressed archive to save space
rm k2_pluspf_20240904.tar.gz

# Confirm the 3 required files are present
ls -lh /projects/tholl_lab_1/daisy_analysis/07_NRdatabase/00_1_kraken2_pluspf/hash.k2d \
        /projects/tholl_lab_1/daisy_analysis/07_NRdatabase/00_1_kraken2_pluspf/opts.k2d \
        /projects/tholl_lab_1/daisy_analysis/07_NRdatabase/00_1_kraken2_pluspf/taxo.k2d
```

**What PlusPF contains:** archaea, bacteria, viruses, plasmids, human,
UniVec_Core, protozoa, and fungi. The fungi component is the most important
for your nutmeg samples — endophytic fungi are the most common plant RNA-seq
contaminant.

---

## Step 2 — Confirm Kraken2 is installed in your conda environment

Kraken2 is installed via conda into the `rnaseq` environment. If you have not
done this yet, run once on a login node:

```bash
conda activate rnaseq
conda install -c bioconda kraken2

# Confirm it installed — you should see a version number printed
kraken2 --version
```

You only need to do this **once**. The script activates the `rnaseq`
environment automatically at the start of every job.

---

## Step 3 — Submit the job

From the root of your `rnaseq` repo:

```bash
sbatch scripts/01_qc/kraken2_prefilter.sbatch
```

You should see: `Submitted batch job 12345678` (your job ID will differ).

---

## Step 4 — Monitor the job

Check if the job is running or queued:
```bash
squeue -u daisycortesj
```

Watch the live log as it runs:
```bash
tail -f /projects/tholl_lab_1/daisy_analysis/02_kraken_filtered/kraken2_prefilter_<JOBID>.out
```
Replace `<JOBID>` with the number from the `Submitted batch job` message.

Check CPU and memory usage on the node:
```bash
jobload <JOBID>
```

---

## Step 5 — Check the output

After the job finishes, your output directory should contain:

```
02_kraken_filtered/
├── MFF1_1_filtered.fq.gz        ← unclassified R1 (feed to Trinity)
├── MFF1_2_filtered.fq.gz        ← unclassified R2 (feed to Trinity)
├── MFF1_classified_1.fq.gz      ← classified R1 (contamination, saved)
├── MFF1_classified_2.fq.gz      ← classified R2 (contamination, saved)
├── MFF1_kraken_output.txt.gz    ← per-read classification (compressed)
├── MFF1_kraken_report.txt       ← per-taxon summary (human-readable)
├── ... (same for MFF2, MFF3, MFL1, MFL2, MFL3)
└── contamination_summary.tsv    ← comparison table across all samples
```

---

## Step 6 — Interpret contamination_summary.tsv

View the table:
```bash
cat /projects/tholl_lab_1/daisy_analysis/02_kraken_filtered/contamination_summary.tsv
```

It looks like this:

| Sample | Total_Reads | Unclassified_pct | Classified_pct | Bacterial_pct | Fungal_pct | Viral_pct |
|--------|-------------|-----------------|----------------|---------------|------------|-----------|
| MFF1   | 42000000    | 78.5            | 21.5           | 8.2           | 12.1       | 0.4       |
| MFL1   | 38000000    | 93.1            | 6.9            | 3.4           | 2.8        | 0.2       |
| ...    | ...         | ...             | ...            | ...           | ...        | ...       |

**What to look for:**
- **Unclassified_pct** — this is the fraction of reads you *keep* for Trinity.
  A healthy plant sample should be 85–98% unclassified with PlusPF (PlusPF
  does not include any plant genomes, so nutmeg reads will not be classified).
- **Classified_pct** — this is removed contamination. If MFF is significantly
  higher than MFL here, it confirms the GC% hypothesis.
- **Fungal_pct** — likely the biggest contributor to contamination in fruit
  tissue. Endophytic fungi are common in nutmeg (Myristica fragrans).
- Compare MFF rows vs MFL rows — you expect MFF to show more contamination.

To look deeper at what specific organisms were found in a sample:
```bash
# Show the top 20 classified taxa for MFF1 sorted by read count
sort -t$'\t' -k2 -nr /projects/tholl_lab_1/daisy_analysis/02_kraken_filtered/MFF1_kraken_report.txt | head -20
```

---

## Step 7 — Update trinity_pooled.sbatch and rerun Trinity

Open `scripts/03_assembly/trinity_pooled.sbatch` and change **two things**:

1. **Input directory** — point to the filtered reads:
```bash
# Before:
INPUT_DIR="/projects/tholl_lab_1/daisy_analysis/01_processed/00_3_MF"

# After:
INPUT_DIR="/projects/tholl_lab_1/daisy_analysis/02_kraken_filtered"
```

2. **File suffix** — update the suffix from `_clean` to `_filtered`:
```bash
# Before:
${INPUT_DIR}/MFF1_1_clean.fq.gz

# After:
${INPUT_DIR}/MFF1_1_filtered.fq.gz
```
Do this for all 12 file paths in the `LEFT_READS` and `RIGHT_READS` blocks.

Then submit Trinity when you're ready:
```bash
sbatch scripts/03_assembly/trinity_pooled.sbatch
```

---

## Important notes

- **Do not delete your `_clean.fq.gz` files.** Even after this succeeds, keep
  the original fastp-cleaned reads as a backup. Disk space permitting, always
  keep one copy of each processing stage.

- **"Classified" does not always mean contamination.** If a read happened to
  match a plant genome that *is* in PlusPF (e.g., human sequences that are
  present for decontamination purposes), it would appear classified even though
  it's not a contaminant. Because PlusPF does not include nutmeg or any closely
  related plant, this effect should be very small — but it's worth knowing.

- **Expected runtime:** ~1–2 hours for all 6 samples with 16 threads. The
  12-hour SLURM time request is a generous safety buffer.

---

## Reference

Wood DE, Lu J, Langmead B. 2019. Improved metagenomic analysis with Kraken 2.
*Genome Biology* 20:257. [DOI 10.1186/s13059-019-1891-0](https://doi.org/10.1186/s13059-019-1891-0)

Pre-built Kraken2 databases: https://benlangmead.github.io/aws-indexes/k2
Kraken2 manual: https://github.com/DerrickWood/kraken2/wiki/Manual
