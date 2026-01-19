# RNA-seq Pipeline: Beginner's Guide

## What is RNA-seq? (The Big Picture)

**RNA-seq** stands for "RNA sequencing." Think of it like this:

- Your **DNA** is like a library of instruction books (genes) that tell your cells how to work
- **RNA** is like photocopies of specific pages from those books that your cells are actually using right now
- **RNA-seq** is a technology that reads these "photocopies" to see which genes are active (turned on) and how much they're being used

**Why is this useful?**
- Compare healthy vs. diseased cells to see which genes are different
- Understand how cells respond to treatments
- Discover new genes or gene variants
- Study development, stress responses, or any biological process

---

## What This Pipeline Does (Simple Explanation)

This pipeline takes raw RNA-seq data (millions of tiny DNA/RNA fragments) and processes them to answer biological questions. It's like having a smart assistant that:

1. **Organizes** the data (builds an index)
2. **Matches** your data to a reference (aligns reads)
3. **Assembles** new sequences when no reference exists (Trinity)

The pipeline has **three main jobs** that work together:

---

## The Three Workflows Explained

### Workflow 1: Building a Genome Index (STAR Index)

**What it does:** Creates a searchable "phone book" from a reference genome

**Real-world analogy:** 
Imagine you have a huge dictionary (the reference genome) and you want to find words quickly. Instead of reading the whole dictionary every time, you create an index (like the alphabetical guide at the top of dictionary pages). This index tells you where everything is located.

**What you need:**
- A reference genome file (`.fa` or `.fasta`) - the "dictionary"
- Optionally: A gene annotation file (`.gtf`) - like a table of contents showing where genes are

**What it produces:**
- An index directory that STAR can use to quickly find where reads belong

**When to use it:**
- First step when you have a reference genome
- Only needs to be done once per genome
- Takes a while (hours), but you only do it once

**Example:**
```
"Take this carrot genome file and create a searchable index so I can quickly find where any DNA sequence belongs"
```

---

### Workflow 2: Aligning Reads (STAR Alignment)

**What it does:** Matches your RNA-seq reads to locations in the reference genome

**Real-world analogy:**
You have thousands of puzzle pieces (your RNA-seq reads) and a completed puzzle picture (the reference genome). The alignment process figures out where each piece belongs in the picture.

**What you need:**
- A pre-built genome index (from Workflow 1)
- Your RNA-seq data files (`.fq` or `.fastq` files, usually compressed as `.gz`)
  - `reads-left` (R1): Forward reads
  - `reads-right` (R2): Reverse reads (for paired-end sequencing)

**What it produces:**
- **BAM file**: A file showing where each read mapped to the genome
- **Gene counts table**: How many reads mapped to each gene (if you have gene annotations)
- **Statistics**: How well the alignment worked

**When to use it:**
- When you have a reference genome AND your RNA-seq data
- This is the main analysis step for most RNA-seq projects
- Run this for each sample you have

**Example:**
```
"Take my RNA-seq data from this carrot sample and figure out which genes are being expressed and how much"
```

---

### Workflow 3: De Novo Assembly (Trinity)

**What it does:** Builds a transcriptome (collection of all RNA sequences) from scratch when you don't have a reference genome

**Real-world analogy:**
Imagine you have a box of puzzle pieces but NO picture to guide you. Trinity is like a smart algorithm that figures out how the pieces fit together to reconstruct the picture, even though you've never seen it before.

**What you need:**
- Your RNA-seq data files (paired-end reads)
  - `reads-left` (R1): Forward reads
  - `reads-right` (R2): Reverse reads

**What it produces:**
- **Trinity.fasta**: A file containing all the assembled transcripts (RNA sequences)
- **Gene-to-transcript mapping**: Which transcripts belong to which genes

**When to use it:**
- When you DON'T have a reference genome for your organism
- For non-model organisms (organisms that haven't been fully sequenced)
- This is computationally intensive and can take days

**Example:**
```
"I have RNA-seq data from a rare plant that has no genome sequence. Build me a transcriptome from scratch"
```

---

## How the Pipeline Works (Step-by-Step)

### The Big Picture Flow

```
┌─────────────────────────────────────────────────────────┐
│                    YOUR RNA-SEQ DATA                     │
│         (Millions of short DNA/RNA fragments)            │
└────────────────────┬──────────────────────────────────────┘
                     │
                     ▼
        ┌────────────────────────┐
        │  Do you have a         │
        │  reference genome?     │
        └───┬──────────────┬──────┘
            │              │
         YES│              │NO
            │              │
            ▼              ▼
    ┌──────────────┐  ┌──────────────┐
    │  Workflow 1: │  │ Workflow 3:  │
    │  Build Index │  │   Trinity    │
    │   (STAR)     │  │  Assembly    │
    └──────┬───────┘  └──────┬───────┘
           │                  │
           ▼                  │
    ┌──────────────┐          │
    │  Workflow 2: │          │
    │   Align      │          │
    │   (STAR)     │          │
    └──────┬───────┘          │
           │                  │
           └────────┬─────────┘
                    │
                    ▼
            ┌──────────────┐
            │   RESULTS    │
            │  (BAM files, │
            │  gene counts,│
            │  transcripts)│
            └──────────────┘
```

### Detailed Step-by-Step

#### Scenario A: You Have a Reference Genome

1. **Step 1: Build Index** (One time, takes hours)
   - Input: Reference genome FASTA file
   - Process: STAR creates a searchable index
   - Output: Index directory
   - Command example: `--mode index --fasta genome.fa --outdir index/`

2. **Step 2: Align Reads** (Repeat for each sample)
   - Input: Index + your RNA-seq reads
   - Process: STAR matches reads to genome locations
   - Output: BAM files + gene counts
   - Command example: `--mode align --genome-index index/ --reads-left R1.fq.gz --reads-right R2.fq.gz`

#### Scenario B: No Reference Genome

1. **Trinity Assembly** (One time per sample, takes days)
   - Input: Your RNA-seq reads
   - Process: Trinity assembles transcripts from scratch
   - Output: Assembled transcriptome (Trinity.fasta)
   - Command example: `--mode trinity --reads-left R1.fq.gz --reads-right R2.fq.gz --outdir trinity/`

---

## Key Concepts Explained Simply

### What are "reads"?
- When you sequence RNA, the machine breaks it into millions of tiny pieces (usually 50-300 letters long)
- Each piece is called a "read"
- You get two files: R1 (forward) and R2 (reverse) because sequencing reads both ends

### What is a "reference genome"?
- A complete DNA sequence of an organism that scientists have already figured out
- Like having the answer key to a test
- Examples: Human genome, mouse genome, Arabidopsis (plant) genome

### What is "alignment"?
- The process of figuring out where each read came from in the genome
- Like matching puzzle pieces to their location in the picture

### What is "gene expression"?
- How much a gene is being used (turned on)
- Measured by counting how many reads map to that gene
- More reads = gene is more active

### What is "de novo assembly"?
- Building something from scratch without a reference
- "De novo" means "from the beginning" in Latin
- Like solving a puzzle without the box picture

---

## How to Explain This to Someone Else

### The Elevator Pitch (30 seconds)
"This pipeline processes RNA-seq data to understand which genes are active in biological samples. It can either match reads to a known genome (like using GPS with a map) or assemble transcripts from scratch when no genome exists (like building a map from scratch)."

### The Detailed Explanation (5 minutes)

**Start with the biological question:**
"We want to know which genes are turned on in our samples. For example, comparing healthy vs. diseased tissue, or before vs. after treatment."

**Explain the data:**
"RNA-seq gives us millions of short DNA sequences (reads) that represent active genes. But these are just fragments - we need to figure out what they mean."

**Explain the two paths:**

**Path 1 (with reference):**
"If we have a reference genome (like a complete map), we:
1. First build an index (create a searchable database) - this is like creating a GPS system
2. Then align our reads (match each fragment to its location) - this is like using GPS to find where you are on the map
3. Count reads per gene to measure expression"

**Path 2 (without reference):**
"If we don't have a reference genome, we use Trinity to:
1. Assemble the fragments into complete transcripts (like solving a jigsaw puzzle without the box)
2. This creates a new transcriptome we can use for analysis"

**End with the output:**
"Either way, we end up with data showing which genes are expressed and how much, which we can then use for statistical analysis and visualization."

---

## Common Questions

### Q: Why do I need to build an index first?
**A:** The index makes alignment much faster. Without it, STAR would have to search through the entire genome for every single read, which would take forever. The index is like a phone book - instead of calling every number to find someone, you look them up alphabetically.

### Q: Can I skip the index step?
**A:** No, if you're using STAR alignment, you need the index. However, you only build it once per genome, and then reuse it for all your samples.

### Q: How long does each step take?
**A:** 
- Index building: Hours (depends on genome size)
- Alignment: Minutes to hours per sample (depends on data size)
- Trinity: Days (very computationally intensive)

### Q: What if my organism has no reference genome?
**A:** Use Trinity (Workflow 3). It will assemble a transcriptome from your reads. You can then use this transcriptome as a "reference" for future analyses.

### Q: What's the difference between single-end and paired-end reads?
**A:** 
- **Single-end**: You sequence one end of each fragment
- **Paired-end**: You sequence both ends (R1 and R2)
- Paired-end is better because you get more information and can assemble more accurately

---

## File Types You'll Encounter

- **FASTA (.fa, .fasta)**: Text file containing DNA/RNA sequences
- **FASTQ (.fq, .fastq)**: Text file containing sequences with quality scores (usually compressed as .gz)
- **GTF/GFF**: Annotation file showing where genes are located in the genome
- **BAM**: Binary file containing aligned reads (compressed version of SAM)
- **TAB/TSV**: Tab-separated text files (like spreadsheets) containing gene counts

---

## The Pipeline in Action: Real Example

Let's say you're studying carrots and want to compare gene expression in roots vs. leaves:

**Step 1: Build index (once)**
```bash
python -m rna_pipeline.cli \
  --mode index \
  --fasta carrot_genome.fa \
  --gtf carrot_genes.gtf \
  --outdir carrot_index
```

**Step 2: Align each sample**
```bash
# Root sample
python -m rna_pipeline.cli \
  --mode align \
  --genome-index carrot_index \
  --reads-left root_R1.fq.gz \
  --reads-right root_R2.fq.gz \
  --outdir root_alignment

# Leaf sample
python -m rna_pipeline.cli \
  --mode align \
  --genome-index carrot_index \
  --reads-left leaf_R1.fq.gz \
  --reads-right leaf_R2.fq.gz \
  --outdir leaf_alignment
```

**Step 3: Analyze results**
- Compare gene counts between root and leaf
- Find genes that are more expressed in one tissue
- Do statistical analysis to find significant differences

---

## Summary

This pipeline is a tool that:
1. **Organizes** reference genomes for fast searching (index)
2. **Matches** your data to known genomes (alignment)
3. **Builds** transcriptomes when no reference exists (Trinity)

It takes raw sequencing data and turns it into interpretable biological information about which genes are active and how much.

The beauty of this pipeline is that it automatically figures out what to do based on what files you give it - you don't need to be a bioinformatics expert to use it!


