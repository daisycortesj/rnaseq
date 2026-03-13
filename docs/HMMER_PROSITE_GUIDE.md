# HMMER & PROSITE Database Setup and Usage Guide

## Overview

After BLAST analysis, you can run additional protein domain and motif analyses to better characterize your candidate genes:

- **HMMER (Pfam)**: Identifies protein domains and families (structural/functional units)
- **PROSITE**: Identifies short functional motifs and signatures (active sites, binding sites)

---

## ⚠️ PREREQUISITES: Install Software Tools First

Before you can use HMMER and PROSITE, you need to **install the software** in your conda environment.

### Quick Check: Are Tools Installed?

```bash
conda activate rnaseq
which hmmscan
which hmmpress
```

If these show paths, tools are installed. If "command not found", you need to install them.

### Install Missing Tools (5-10 minutes)

**Option 1: Update from environment.yml (Recommended)**

Your `environment.yml` now includes HMMER and EMBOSS:

```bash
conda activate rnaseq
conda env update -f environment.yml --prune
```

**Option 2: Install Only New Tools (Faster)**
```bash
conda activate rnaseq
conda install -c bioconda hmmer emboss
```

📖 **Detailed update instructions**: See `UPDATE_ENVIRONMENT.md`

### Verify Installation

After installation, check that tools are available:

```bash
# Check HMMER
hmmscan -h
hmmpress -h

# Check EMBOSS (for PROSITE)
ps_scan.pl -h
```

✓ Once tools are installed, proceed to download databases below.

---

## When to Use Each Tool

### Use HMMER (Pfam) when:
- ✅ You want to identify **protein domains** (larger functional units, 50-500 aa)
- ✅ You want to classify proteins into **families** (e.g., P450, kinase, etc.)
- ✅ You need **sensitive** domain detection (more sensitive than BLAST)
- ✅ You want to understand **protein architecture**

### Use PROSITE when:
- ✅ You want to find **short functional motifs** (5-50 aa)
- ✅ You're looking for **active sites** or **binding sites**
- ✅ You want to identify **post-translational modification sites**
- ✅ You need **highly specific** pattern matching

### Recommendation:
**Run HMMER first** - it's more commonly used and provides better overall functional annotation. PROSITE is complementary but optional.

---

## Step 1: Download Databases

Run the automated download script:

```bash
sbatch scripts/setup/download_protein_databases.sbatch
```

### What it downloads:

1. **Pfam-A.hmm** (~600 MB compressed → ~2.3 GB uncompressed)
   - Location: `/projects/tholl_lab_1/daisy_analysis/07_NRdatabase/Pfam-A.hmm`
   - Contains ~20,000 protein families
   - Automatically indexed with `hmmpress`

2. **PROSITE** (~30 MB)
   - Location: `/projects/tholl_lab_1/daisy_analysis/07_NRdatabase/prosite/`
   - Contains ~2,000 patterns and profiles
   - Files: `prosite.dat` and `prosite.psa`

### Runtime:
- **10-60 minutes** depending on download speed
- Much faster than BLAST database downloads!

---

## Step 2: Run HMMER (Pfam) Analysis

### Input Requirements:
You need **protein sequences** (FASTA format). Generate them with:

```bash
# Translate your CDS to proteins first
sbatch scripts/06_blast/run_translate_cds.sbatch DC
```

This creates: `/projects/tholl_lab_1/daisy_analysis/06_analysis/blast_input_DC/all_genes_protein.fasta`

### Run HMMER on all proteins:

```bash
sbatch scripts/07_domains/run_hmmer.sbatch DC all_genes_protein.fasta
```

### Or run on filtered candidates only:

```bash
# First filter your BLAST results to candidates
python scripts/filter_combined_results.py \
    --input combined_blast_results_DC.tsv \
    --output filtered_candidates.txt

# Then run HMMER on just those proteins
sbatch scripts/07_domains/run_hmmer.sbatch DC filtered_proteins.fasta
```

### Output Files:

1. **`all_genes_protein_pfam_domains.txt`** - Domain table (tab-separated)
   - Columns: target name, accession, query name, E-value, score, bias, etc.
   - Use this for downstream analysis

2. **`all_genes_protein_pfam_domains.out`** - Human-readable output
   - Detailed alignments and descriptions
   - Use this for manual inspection

### Runtime:
- **30 min - 6 hours** depending on number of sequences
- ~50,000 proteins: 2-4 hours with 16 CPUs

---

## Step 3: (Optional) Run PROSITE Analysis

### Run PROSITE scan:

```bash
sbatch scripts/07_domains/run_prosite.sbatch DC all_genes_protein.fasta
```

### Output:
- **`all_genes_protein_prosite_motifs.txt`** - Motif table

### Runtime:
- **30 min - 4 hours** depending on number of sequences

---

## Understanding the Results

### HMMER Results Interpretation

Example line from `*_pfam_domains.txt`:

```
PF00067.23  Cytochrome_P450  gene_12345  1.2e-45  150.3  domain:125-487
```

This means:
- Gene `gene_12345` contains a **Cytochrome P450 domain** (PF00067)
- **E-value: 1.2e-45** (very significant, <1e-5 is good)
- **Score: 150.3** (higher = better match)
- Domain location: amino acids 125-487

### Key Pfam Domains for Plant Terpene Biosynthesis:

| Pfam ID | Domain Name | Description |
|---------|-------------|-------------|
| PF00067 | Cytochrome_P450 | P450 enzymes (terpene modification) |
| PF01397 | Terpene_synth | Terpene synthase (TPS) |
| PF03936 | Terpene_synth_C | TPS C-terminal domain |
| PF00432 | Prenyltransferase | Prenyl transferase (terpene backbone) |
| PF08240 | ADH_N | Alcohol dehydrogenase (terpene modification) |
| PF00155 | ADP_ribosylation | NAD-binding domain |

### PROSITE Results Interpretation

Look for functional motifs like:
- **PS00101**: Active sites
- **PS00102**: Binding sites
- **PS00103**: Phosphorylation sites

---

## Workflow Comparison

### Option A: Full Analysis (Recommended)

```bash
# 0. Install tools first (if not already installed)
conda env update -f environment.yml --prune       # 5-10 min

# 1. Download databases
sbatch scripts/setup/download_protein_databases.sbatch  # 10-60 min

# 2. Translate CDS to protein (if not done)
sbatch scripts/06_blast/run_translate_cds.sbatch DC       # 5-10 min

# 3. Run BLAST (you already did this)
sbatch scripts/06_blast/run_blastp_strict.sbatch DC     # 2-12 hours

# 4. Run HMMER for domain identification
sbatch scripts/07_domains/run_hmmer.sbatch DC all_genes_protein.fasta  # 2-4 hours

# 5. (Optional) Run PROSITE
sbatch scripts/07_domains/run_prosite.sbatch DC all_genes_protein.fasta  # 1-4 hours
```

### Option B: Filtered Candidates Only (Faster)

```bash
# 0. Install tools first (if not already installed)
conda env update -f environment.yml --prune

# 1. Download databases
sbatch scripts/setup/download_protein_databases.sbatch

# 2. Filter BLAST results to top candidates
python scripts/filter_combined_results.py \
    --input combined_blast_results_DC.tsv \
    --output filtered_candidates.txt \
    --top_n 500

# 3. Extract those proteins to new FASTA
# (you'll need to create this, or run HMMER on filtered list)

# 4. Run HMMER on filtered set (much faster)
sbatch scripts/07_domains/run_hmmer.sbatch DC filtered_proteins.fasta  # 10-30 min
```

---

## Combining Results

To get a complete picture, combine BLAST + HMMER results:

### Method: Merge by gene ID

1. BLAST tells you: **"This gene is similar to Protein X"**
2. HMMER tells you: **"This gene contains Domain Y"**

Example combined interpretation:
```
gene_12345:
  - BLAST hit: "Limonene synthase" (E-value: 1e-50)
  - Pfam domain: "Terpene_synth" (PF01397, E-value: 1e-45)
  → Conclusion: High-confidence limonene synthase candidate
```

### Python script to merge (you can create this):

```python
# merge_blast_hmmer.py
import pandas as pd

# Load results
blast = pd.read_csv("combined_blast_results_DC.tsv", sep="\t")
hmmer = pd.read_csv("all_genes_protein_pfam_domains.txt", sep="\t", comment="#")

# Merge on gene ID
merged = blast.merge(hmmer, left_on="qseqid", right_on="query_name", how="left")

# Save
merged.to_csv("combined_blast_hmmer_results.tsv", sep="\t", index=False)
```

---

## Best Practices

### 1. Run on Non-Filtered List First
- ✅ Run HMMER on **all proteins** to avoid missing domains
- ✅ Filter after combining BLAST + HMMER results
- ❌ Don't filter before HMMER - you might miss important domains

### 2. Use Appropriate E-value Thresholds
- HMMER: **E < 1e-5** (standard for domains)
- PROSITE: **Score > 10** (pattern-dependent)

### 3. Prioritize HMMER over PROSITE
- HMMER is more comprehensive and sensitive
- PROSITE adds detail but isn't required

### 4. Verify Important Hits Manually
- Check domain architecture makes sense
- Look for expected domain combinations
- Compare to known proteins in literature

---

## Troubleshooting

### Installation Issues

**Problem**: Can't activate conda environment
```bash
# Reinitialize conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq
```

**Problem**: "conda: command not found"
```bash
# Check if conda is in PATH
which conda

# If not found, add to PATH (adjust path if needed)
export PATH="$HOME/miniconda3/bin:$PATH"
source ~/miniconda3/etc/profile.d/conda.sh
```

**Problem**: Installation fails with "package not found"
```bash
# Update conda channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# Try installation again
conda install -c bioconda hmmer
```

**Problem**: HMMER tools installed but not found
```bash
# Check if in correct environment
conda activate rnaseq
which hmmscan

# If still not found, try reinstalling
conda uninstall hmmer
conda install -c bioconda hmmer
```

**Problem**: Permission denied when running install script
```bash
# Make script executable
chmod +x scripts/install_protein_tools.sh

# Run again
bash scripts/install_protein_tools.sh
```

### Database Download Issues

**Problem**: wget/curl fails
```bash
# Manual download
cd /projects/tholl_lab_1/daisy_analysis/07_NRdatabase
wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
```

**Problem**: HMMER tools not found
```bash
# Install HMMER
conda activate rnaseq
conda install -c bioconda hmmer
```

### HMMER Run Issues

**Problem**: "Pfam database not found"
```bash
# Check if database exists
ls -lh /projects/tholl_lab_1/daisy_analysis/07_NRdatabase/Pfam-A.hmm*

# Should see:
# Pfam-A.hmm
# Pfam-A.hmm.h3m
# Pfam-A.hmm.h3i
# Pfam-A.hmm.h3f
# Pfam-A.hmm.h3p
```

**Problem**: HMMER runs but no results
- Check E-value threshold (default: 1e-5)
- Some proteins genuinely have no known domains
- Check input is protein sequence (not nucleotide)

### PROSITE Issues

**Problem**: ps_scan.pl not found
```bash
# Install EMBOSS
conda activate rnaseq
conda install -c bioconda emboss
```

---

## Quick Reference Commands

```bash
# Check if tools are installed
conda activate rnaseq
which hmmscan hmmpress

# Install tools (if needed)
conda env update -f environment.yml --prune

# Download databases
sbatch scripts/setup/download_protein_databases.sbatch

# Check database status
ls -lh /projects/tholl_lab_1/daisy_analysis/07_NRdatabase/Pfam-A.hmm*
ls -lh /projects/tholl_lab_1/daisy_analysis/07_NRdatabase/prosite/

# Run HMMER
sbatch scripts/07_domains/run_hmmer.sbatch DC all_genes_protein.fasta

# Run PROSITE
sbatch scripts/07_domains/run_prosite.sbatch DC all_genes_protein.fasta

# Check job status
squeue -u $USER

# View results
less all_genes_protein_pfam_domains.txt
less all_genes_protein_prosite_motifs.txt
```

---

## Summary

1. **Download databases**: `sbatch scripts/setup/download_protein_databases.sbatch`
2. **Run HMMER first**: Most useful for functional annotation
3. **Run PROSITE optionally**: Adds detail on specific motifs
4. **Combine with BLAST**: Use both to get complete picture
5. **Prioritize domains**: Use Pfam domains to filter/prioritize candidates

Good luck with your domain analysis! 🧬
