# RNA-seq Pipeline: STAR vs Trinity - Complete Teaching Guide
*For beginners: Building a solid foundation in RNA-seq differential expression*

---

## Table of Contents
1. [Starting Point: FASTQ Files](#1-starting-point-fastq-files)
2. [Quality Control (QC)](#2-quality-control-qc)
3. [The Great Split: STAR vs Trinity](#3-the-great-split-star-vs-trinity)
4. [STAR Pipeline (Reference-Based)](#4-star-pipeline-reference-based)
5. [Trinity Pipeline (De Novo Assembly)](#5-trinity-pipeline-de-novo-assembly)
6. [Differential Expression Analysis](#6-differential-expression-analysis)
7. [Key Concepts Summary](#7-key-concepts-summary)

---

## 1. Starting Point: FASTQ Files

### What is a FASTQ file?

A **FASTQ file** is the raw output from your sequencing machine. It contains millions of short DNA sequences (called "reads") that represent pieces of your RNA molecules.

**Think of it like this:** Imagine you took a book, shredded it into millions of tiny strips (each 50-150 letters long), and now you need to figure out what the book originally said. Each strip is a "read."

### FASTQ Format (Toy Example)

```
@READ_001                           ← Read name (unique ID)
ACTGACTGACTGACTGACTG                ← The actual sequence (A, T, G, C bases)
+                                   ← Separator line
IIIIIIIIIIIIIIIIIII                ← Quality scores (one per base)
```

**Quality scores:** Each letter represents how confident the sequencer was about that base call.
- `I` = very high quality (Phred score ~40, 99.99% accuracy)
- `#` = low quality (Phred score ~2, only 60% accuracy)

### Paired-End Sequencing (Most Common for RNA-seq)

When you do RNA-seq, you usually get **paired-end reads**:
- **R1 (Read 1)**: Sequence from the left/forward end of the fragment
- **R2 (Read 2)**: Sequence from the right/reverse end of the same fragment

```
Your RNA fragment:
5'═══════════════════════════════════3'
   R1 →→→→          ←←←← R2
```

**Why paired-end?** It gives you more information! You know these two reads came from the same RNA molecule and roughly how far apart they are. This helps with alignment accuracy.

**Your files will look like:**
```
sample1_R1.fastq.gz   (forward reads)
sample1_R2.fastq.gz   (reverse reads)
```

---

## 2. Quality Control (QC)

Before you do anything, you need to check if your sequencing data is good quality. Bad data = bad results, no matter how fancy your analysis is!

### Step 2a: FastQC - Initial Quality Check

**What is FastQC?** A tool that generates a visual HTML report showing the quality of your reads.

**What it checks:**
- **Per-base quality scores:** Are most bases high quality?
- **Sequence length distribution:** Are reads the expected length?
- **GC content:** Does it match what you expect for your organism?
- **Adapter content:** Are there leftover sequencing adapters?
- **Duplicate sequences:** Too many duplicates might indicate PCR bias

**Your script reference:** `rna_pipeline/tools/qc.py` lines 74-96
```python
def build_fastqc_cmd(fastq_files: List[str], outdir: str, threads: int):
    """Run FastQC on your FASTQ files"""
    cmd = [
        "fastqc",
        "--outdir", outdir,      # Where to save reports
        "--threads", str(threads),
        "--quiet"
    ]
    cmd.extend(fastq_files)      # Add all your FASTQ files
    return cmd
```

**How to interpret FastQC results:**
- ✅ **PASS** (green): This metric looks good!
- ⚠️ **WARN** (orange): Not ideal but might be okay
- ❌ **FAIL** (red): This is problematic - you likely need to fix it

### Step 2b: MultiQC - Aggregate QC Reports

**What is MultiQC?** After you run FastQC on 10, 20, or 100 samples, you don't want to open 100 separate HTML files! MultiQC combines all FastQC reports into one easy-to-read summary.

**Your script reference:** `rna_pipeline/tools/qc.py` lines 99-121

**What to look for in MultiQC:**
1. **Are all samples similar quality?** If one sample is much worse, it might be a failed library.
2. **Overall read quality:** Mean quality score should be >28-30 across most of the read.
3. **Adapter contamination:** Should be <1% (ideally 0%)

### Step 2c: fastp - Trimming and Cleaning (if needed)

**What is trimming?** If your QC shows problems (low quality bases, adapters), you need to "clean" your reads before using them.

**fastp does:**
1. **Adapter removal:** Cuts off sequencing adapters (leftover junk from library prep)
2. **Quality trimming:** Removes low-quality bases from read ends
3. **Length filtering:** Removes reads that are too short after trimming

**Toy Example - Before and After Trimming:**

```
BEFORE (raw read):
Read: ACTGACTGACTGACTGACTAGATCGGAAGAGC
Quality: IIIIIIIIIIIIIIIII#################
                         ↑ quality drops
                         ↑ adapter sequence

AFTER fastp (trimmed):
Read: ACTGACTGACTGACTGACTGA
Quality: IIIIIIIIIIIIIIIII###
                          ↑ still removed low quality end
                          
Final kept read:
Read: ACTGACTGACTGACTGACTGA
Quality: IIIIIIIIIIIIIIIII
```

**Your script reference:** `rna_pipeline/tools/fastp.py` lines 8-61

**Key parameters:**
- `min_quality`: Default 20 (99% base accuracy)
- `min_length`: Default 50 bp (reads shorter than this are discarded)
- `--detect_adapter_for_pe`: fastp automatically finds and removes adapters!

**When to use fastp:**
- ❌ Your FastQC shows adapter contamination >5%
- ❌ Quality drops below Q20 in the last 10-20 bases
- ❌ You see lots of "overrepresented sequences" that look like adapters
- ✅ **If your FastQC looks good, you can skip trimming!**

---

## 3. The Great Split: STAR vs Trinity

Now you have clean, high-quality reads. Here's where the path diverges based on one question:

### Do You Have a Reference Genome?

```
┌────────────────────────────────────────────┐
│ Do you have a reference genome/            │
│ transcriptome for your organism?           │
└────────────────┬───────────────────────────┘
                 │
         ┌───────┴───────┐
         ↓               ↓
    YES: STAR        NO: Trinity
 (Reference-based)  (De novo assembly)
```

**Reference genome:** A complete, annotated genome sequence for your organism (like human, mouse, Arabidopsis, etc.).

**When to use each:**

| Use STAR if... | Use Trinity if... |
|----------------|-------------------|
| ✅ Well-studied organism (human, mouse, Arabidopsis) | ✅ Non-model organism with no reference |
| ✅ You want faster analysis | ✅ You want to discover novel transcripts |
| ✅ You have a genome annotation (GTF file) | ✅ The reference genome is poor quality |
| ✅ You're studying known genes | ✅ You're studying a new species |

**Your organism:** It looks like you're working with a plant (based on `DC1L1`, `DC1R1` - leaf samples, and `DGR1` - possibly root?). If this is a well-studied plant like Arabidopsis or tomato, use STAR. If it's a wild species or understudied plant, Trinity might be better.

---

## 4. STAR Pipeline (Reference-Based)

**STAR** = **S**pliced **T**ranscripts **A**lignment to a **R**eference

This is the most common approach for RNA-seq if you have a reference genome.

### The STAR Pipeline Overview

```
Clean FASTQ files
       ↓
[Step 1: Build STAR Index]  ← One-time setup per genome
       ↓
   STAR Index
       ↓
[Step 2: Map/Align Reads]    ← Do this for each sample
       ↓
    BAM files (aligned reads)
       ↓
[Step 3: Count Reads per Gene]
       ↓
   Count Matrix (genes × samples)
       ↓
[Step 4: Differential Expression (PyDESeq2)]
       ↓
   DEG Results (which genes are up/down)
```

### Step 4.1: Build STAR Index (One-Time Setup)

**What is indexing?** Think of it like building an index for a book. Instead of reading page-by-page to find a word, you look in the index and jump straight to the right page. STAR builds an index of your genome so it can quickly find where reads align.

**What you need:**
1. **Genome FASTA file** (`.fa` or `.fasta`): The complete DNA sequence
2. **Genome annotation GTF file** (`.gtf`): Tells where genes are located

**Example file names:**
```
Athaliana.genome.fa       ← DNA sequence (all chromosomes)
Athaliana.annotation.gtf  ← Gene locations and structures
```

**What's in a GTF file?** (Toy Example)

```
Chr1  gene  1000  2000  GENE123  ← Gene GENE123 is on Chr1, positions 1000-2000
Chr1  exon  1000  1200  GENE123  ← First exon
Chr1  exon  1800  2000  GENE123  ← Second exon
                 ↑
               intron (not shown, it's the gap)
```

**Your script reference:** `rna_pipeline/tools/star.py` lines 34-69

```python
def build_star_index_cmd(fasta: str, gtf: str, outdir: str, threads: int, readlen: int):
    """
    Build STAR genome index
    
    Key parameters:
    - readlen: Your read length (usually 50-150bp)
    - sjdbOverhang: readlen - 1 (for splice junction detection)
    - genomeSAindexNbases: Adjusted based on genome size
    """
```

**What happens during indexing:**
- STAR reads the genome and creates a special compressed index
- It uses the GTF to mark splice junctions (where exons connect)
- This takes 20-60 minutes and uses 10-30GB RAM
- **You only do this ONCE per genome**

**Output:** A directory with index files (e.g., `genome_index/`)

### Step 4.2: Align/Map Reads to Genome

**What is alignment/mapping?** Taking each of your short reads (50-150 bp) and figuring out where in the genome it came from.

**Why is RNA-seq alignment special?** Unlike DNA-seq, RNA has **introns removed** (splicing). So a read might span an exon-exon junction:

```
Genome DNA:
EXON1═════════INTRON═════════EXON2═════════
5'---1000bp---+----5000bp----+---1000bp---3'

Your RNA read (after splicing):
EXON1════|════EXON2
     Read spans junction!

Regular aligner (like BWA): "This read doesn't match the genome!"
STAR (splice-aware): "Aha! This read spans from exon1 to exon2."
```

**Your script reference:** `rna_pipeline/tools/star.py` lines 72-138

```python
def build_star_align_cmd(genome_index: str, reads_left: str, reads_right: str, 
                         outprefix: str, threads: int, quant_mode: bool = True):
    """
    Align reads to indexed genome
    
    Key outputs:
    - BAM file: Binary file with aligned reads
    - ReadsPerGene.out.tab: Gene counts (if quant_mode=True)
    """
```

**Key parameters:**
- `--genomeDir`: Your STAR index directory
- `--readFilesIn`: Your R1 and R2 FASTQ files
- `--outSAMtype BAM SortedByCoordinate`: Output sorted BAM files
- `--quantMode GeneCounts`: Count reads per gene (we want this!)

**What is a BAM file?** Binary Alignment Map - it contains:
```
Read_001: Aligned to Chr1, position 12345, MAPQ=60 (high confidence)
Read_002: Aligned to Chr3, position 98765, MAPQ=30 (medium confidence)
Read_003: Unmapped (couldn't find where it belongs)
```

**STAR outputs:**
```
sample1_Aligned.sortedByCoord.out.bam  ← Aligned reads
sample1_ReadsPerGene.out.tab           ← Gene counts (THIS IS GOLD!)
sample1_Log.final.out                  ← Alignment statistics
```

### Step 4.3: Build Count Matrix

**What is a count matrix?** A table where:
- **Rows** = genes
- **Columns** = samples  
- **Values** = number of reads that mapped to each gene in each sample

**Toy Example:**

```
            sample1_leaf  sample2_leaf  sample3_root  sample4_root
GENE_A            450           523           892          1045
GENE_B           1023          1156            23            45
GENE_C             12            18          3456          3890
```

**What this tells us:**
- GENE_A: ~500 reads in leaf, ~900 in root → probably higher in root
- GENE_B: ~1100 in leaf, ~35 in root → probably higher in leaf  
- GENE_C: ~15 in leaf, ~3700 in root → MUCH higher in root!

**Where do counts come from?** STAR's `ReadsPerGene.out.tab` file:

```
# STAR output for sample1_leaf:
GENE_A    450    445    5      ← unstranded, first-strand, second-strand
GENE_B    1023   1020   3
GENE_C    12     11     1
```

**Your script reference:** `rna_pipeline/tools/build_count_matrix.py` lines 59-85

```python
def read_star_counts(filepath):
    """
    Read STAR ReadsPerGene.out.tab file
    
    STAR output columns:
    1. gene_id
    2. unstranded counts  ← We use this!
    3. first strand counts
    4. second strand counts
    """
```

**How the count matrix is built:** Lines 86-186

1. Find all `*ReadsPerGene.out.tab` files
2. Parse sample names (e.g., `DC1L1_ReadsPerGene.out.tab` → sample name `DC1L1`)
3. Extract counts for each gene
4. Merge into one big matrix
5. Save as `gene_count_matrix.tsv`

**What you get:**
```
count_matrices/
├── gene_count_matrix.tsv      ← Main output: genes × samples
├── sample_metadata.tsv         ← Sample info (condition, replicate, etc.)
└── count_summary.txt           ← Statistics
```

### Step 4.4: Differential Expression Analysis (PyDESeq2)

**Now the real science begins!** You have counts, but raw counts are hard to compare. You need statistics!

**PyDESeq2** is a Python implementation of DESeq2, a statistical method for finding differentially expressed genes (DEGs).

**The core question:** Is GENE_A **truly** more expressed in root vs leaf, or is the difference just random noise?

#### Key Concepts You MUST Understand:

##### A. Normalization

**Problem:** Sample 1 might have 10 million total reads, Sample 2 might have 20 million. If we compare raw counts, Sample 2 will always look "higher" just because we sequenced it deeper!

**Solution: Size factor normalization**

**Toy Example:**

```
Raw counts:
            sample1  sample2  sample3
GENE_A        100      200      150
GENE_B        200      400      300
GENE_C        300      600      450
Total reads  10M      20M      15M
             ↑        ↑        ↑
         sequenced at different depths!

After normalization (divide by size factors):
Size factors:   1.0      2.0      1.5

Normalized:
            sample1  sample2  sample3
GENE_A        100      100      100   ← Now comparable!
GENE_B        200      200      200
GENE_C        300      300      300
```

**DESeq2's size factor** (simplified): Median of ratios method
- Calculate geometric mean for each gene across samples
- For each sample, divide counts by geometric mean
- Take median of these ratios → size factor

**In your code:** `pydeseq2_analysis.py` lines 464-476

##### B. Dispersion

**Dispersion** = biological variability (how much does the same gene vary between replicates?)

**Why it matters:** With only 2-3 replicates per condition, we need to estimate variance accurately.

**Toy Example:**

```
GENE_A in leaf samples:
Rep1: 100 reads
Rep2: 110 reads
Rep3: 105 reads
Mean: 105, Variance: small → LOW dispersion (reliable!)

GENE_B in leaf samples:
Rep1: 50 reads
Rep2: 200 reads
Rep3: 30 reads
Mean: 93, Variance: HUGE → HIGH dispersion (noisy!)
```

**DESeq2's approach:**
1. Estimate dispersion for each gene
2. Shrink dispersions toward a fitted curve (borrowing information across genes)
3. Use shrunk dispersions for statistical tests

**Visual concept:**

```
Dispersion vs Mean Expression
    ↑
Disp|     •     Gene-wise estimates
    |   • • •
    | • • • • •  ← Fitted trend (red line)
    |• • • • • • •
    |─────────────→
           Mean expression
```

##### C. Statistical Testing

**For each gene, DESeq2 asks:** Is the difference in expression between conditions **statistically significant**?

**It uses:** Negative binomial test (appropriate for count data)

**Outputs for each gene:**

1. **baseMean**: Average normalized expression across all samples
   - Example: baseMean = 1000 → moderately expressed gene
   - Example: baseMean = 10 → lowly expressed gene

2. **log2FoldChange (log2FC)**: Effect size
   - **log2FC = +1** → 2× higher in condition A (root) than B (leaf)
   - **log2FC = +2** → 4× higher in A than B
   - **log2FC = -1** → 2× lower in A than B (i.e., 2× higher in B)

   **Why log2?** It's symmetric:
   ```
   Raw fold changes: 0.5×  1×  2×  4×  8×
   log2 fold changes: -1   0  +1  +2  +3
                       ↑       ↑
                    symmetric around 0
   ```

3. **pvalue**: Probability that difference is due to chance
   - p = 0.001 → 0.1% chance this is random (strong evidence!)
   - p = 0.05 → 5% chance random (conventional cutoff)
   - p = 0.3 → 30% chance random (not significant)

4. **padj (adjusted p-value / FDR)**: Corrected for multiple testing
   - **Why adjust?** You're testing 20,000 genes. By chance, 5% (1000 genes) will have p < 0.05 even if nothing is real!
   - **FDR (False Discovery Rate)**: Controls the proportion of false positives
   - **padj < 0.05** → Less than 5% of your "significant" genes are false positives

**Example Results Table:**

```
gene_id    baseMean  log2FC  pvalue   padj     Interpretation
GENE_A     1000      +3.2    1e-50   1e-47    *** Strongly up in root
GENE_B     500       -2.1    1e-20   1e-18    *** Strongly up in leaf  
GENE_C     200       +0.5    0.01    0.03     * Weakly up in root (but significant)
GENE_D     100       +1.0    0.06    0.15     Not significant (padj > 0.05)
GENE_E     10        +5.0    0.001   0.10     High fold change but noisy (low counts)
```

##### D. Filtering DEGs (Differentially Expressed Genes)

**Common criteria:**
- **padj < 0.05** (statistically significant)
- **|log2FC| > 1** (at least 2-fold change) or **|log2FC| > 2** (4-fold change)

**Why filter by fold change too?** A gene might be statistically significant (padj < 0.05) but biologically boring (only 1.1× difference).

**Your script:** `pydeseq2_analysis.py` lines 491-519

```python
# Filter DEGs
deg_df = results_df[results_df['padj'] < padj_cutoff]  # Significant
deg_df = deg_df[np.abs(deg_df['log2FoldChange']) > lfc_cutoff]  # Biologically meaningful
```

##### E. Heatmaps and Visualization

**Heatmap:** Visual representation of gene expression across samples

**Toy Example (3 genes, 4 samples):**

```
            L1    L2    R1    R2    ← Samples (L=leaf, R=root)
GENE_A     low   low   high  high   
GENE_B     high  high  low   low    
GENE_C     med   med   med   med    

Heatmap:
         L1   L2 | R1   R2
GENE_A  🔵  🔵 | 🔴  🔴     🔵 = low expression (blue)
GENE_B  🔴  🔴 | 🔵  🔵     🔴 = high expression (red)
GENE_C  ⚪  ⚪ | ⚪  ⚪     ⚪ = medium expression (white)
```

**Scaling:** To compare patterns across genes, we often **center** or **z-score** each row (gene):

```
Before centering (raw log2 values):
            L1    L2    R1    R2   Mean
GENE_A      2     3     8     9    5.5
GENE_B     10    11     4     5    7.5

After centering (subtract row mean):
            L1    L2    R1    R2
GENE_A     -3.5  -2.5  +2.5  +3.5  ← Pattern: low in L, high in R
GENE_B     +2.5  +3.5  -3.5  -2.5  ← Pattern: high in L, low in R
```

**Z-score** (one step further): Divide by standard deviation too, so all genes are on the same scale.

**Your script:** `pydeseq2_analysis.py` lines 584-600

---

## 5. Trinity Pipeline (De Novo Assembly)

**Trinity** is used when you **don't have a reference genome**. Instead of aligning reads to a genome, Trinity **assembles** the reads into transcripts.

### The Trinity Pipeline Overview

```
Clean FASTQ files
       ↓
[Step 1: Trinity Assembly]  ← Computationally intensive!
       ↓
   Trinity.fasta (assembled transcripts)
       ↓
[Step 2: Abundance Estimation (RSEM/Kallisto)]
       ↓
   Count Matrix (transcripts × samples)
       ↓
[Step 3: Differential Expression (DESeq2/edgeR)]
       ↓
   DEG Results
```

### What is De Novo Assembly?

**Analogy:** Remember the shredded book? With STAR, you had the original book to compare against. With Trinity, you don't! You have to reconstruct the book from overlapping strips.

**Example:**

```
Your reads (short fragments):
Read1: ACTGACTGACTGA
Read2:      CTGACTGACCCTA
Read3:           TGACCCTAGGATA
Read4:                  CTAGGATAGGGG

Trinity finds overlaps and assembles:
ACTGACTGACTGACCCTAGGATAGGGG
       └─overlap─┘
```

**Challenges:**
- **Isoforms:** One gene can produce multiple transcript variants (alternative splicing)
- **Repeats:** Similar sequences from different genes can be confusing
- **Errors:** Sequencing errors can create false branches

### Step 5.1: Trinity Assembly

**Your script reference:** `rna_pipeline/tools/trinity.py` lines 2-35

```python
def build_trinity_cmd(left: str, right: str, outdir: str, threads: int, mem_gb: int):
    """
    Build Trinity de novo assembly command
    
    Trinity does three main steps internally:
    1. Inchworm: Build initial contigs from high-coverage reads
    2. Chrysalis: Cluster contigs and build de Bruijn graphs
    3. Butterfly: Resolve alternative splicing and output transcripts
    """
```

**Input:** R1 and R2 FASTQ files (all samples can be combined or assembled separately)

**Output:** `Trinity.fasta` - assembled transcripts

**Example Trinity.fasta:**

```
>TRINITY_DN1000_c0_g1_i1 len=1523
ATGACTGACTGACTG...TAATAG  (1523 bases)
>TRINITY_DN1000_c0_g1_i2 len=1350
ATGACTGACTGACTG...TAATAG  (1350 bases, isoform 2)
>TRINITY_DN1001_c0_g1_i1 len=892
CCCTAGGATACCCTA...AAAAAA  (892 bases)
```

**Trinity naming:**
- `DN1000`: Component number
- `g1`: Gene within component
- `i1`, `i2`: Isoform 1, isoform 2

**What are isoforms?** Different versions of the same gene due to alternative splicing:

```
Gene structure:
EXON1═══INTRON1═══EXON2═══INTRON2═══EXON3═══INTRON3═══EXON4

Isoform 1 (full length):
EXON1═══EXON2═══EXON3═══EXON4

Isoform 2 (skips exon 3):
EXON1═══EXON2═══EXON4

Isoform 3 (alternative exon 2):
EXON1═══EXON2b═══EXON3═══EXON4
```

**Trinity Resources:**
- **Time:** 6-48 hours (depends on data size)
- **RAM:** 30-100+ GB
- **CPU:** More threads = faster

**Trinity is expensive!** This is why you use STAR if you can.

### Step 5.2: Quantification (After Assembly)

After Trinity gives you assembled transcripts, you need to count how many reads from each sample map to each transcript.

**Common tools:**
1. **RSEM** (RNA-Seq by Expectation-Maximization)
2. **Kallisto** (fast pseudoalignment)
3. **Salmon** (fast quasi-mapping)

**Process:**
```
Trinity.fasta (reference transcriptome)
     +
Sample1_R1.fq, Sample1_R2.fq
     ↓ [RSEM align_and_estimate_abundance]
Sample1.isoforms.results
```

**Output:** Counts for each transcript in each sample

### Step 5.3: Build Count Matrix (Similar to STAR)

Same concept as STAR, but with transcripts instead of genes.

### Trinity vs STAR Summary

| Feature | STAR | Trinity |
|---------|------|---------|
| **Needs reference?** | Yes | No |
| **Speed** | Fast (~30 min/sample) | Slow (6-48 hours total) |
| **Memory** | 30GB | 50-100GB |
| **Output** | Aligned reads (BAM) + counts | Assembled transcripts (FASTA) |
| **Gene structure** | Uses known genes (GTF) | Discovers transcripts de novo |
| **Isoforms** | Can detect with GTF | Discovers novel isoforms |
| **Best for** | Model organisms | Non-model organisms |
| **Difficulty** | Easier | More complex |

---

## 6. Differential Expression Analysis

**This is the same for both STAR and Trinity!** Once you have a count matrix, the analysis is identical.

### Workflow

```
Count Matrix (genes/transcripts × samples)
       ↓
Sample Metadata (which samples are leaf? which are root?)
       ↓
[PyDESeq2 / DESeq2 / edgeR]
       ↓
Statistical Analysis:
  - Normalize counts
  - Estimate dispersion
  - Test for differential expression
       ↓
Results:
  - Which genes are up/down?
  - By how much (fold change)?
  - How confident are we (p-value)?
       ↓
Biological Interpretation:
  - What pathways are enriched?
  - What does this mean for my biological question?
```

### Your PyDESeq2 Script

**Reference:** `pydeseq2_analysis.py` (the whole file!)

**Key steps in your script:**

1. **Read data** (lines 78-99)
   - Load count matrix
   - Load sample metadata

2. **Create DESeq dataset** (lines 212-270)
   - Filter low-count genes (< 10 reads total)
   - Set up statistical model

3. **Run DESeq2** (lines 280-309)
   - Estimate size factors (normalization)
   - Estimate dispersions (variability)
   - Fit negative binomial model

4. **Statistical testing** (lines 312-380)
   - Compute p-values for each gene
   - Adjust for multiple testing (FDR)

5. **Generate visualizations** (lines 383-826)
   - MA plot: log fold change vs mean expression
   - Volcano plot: log fold change vs significance
   - Heatmaps: Gene expression patterns

---

## 7. Key Concepts Summary

### Glossary

**FASTQ:** Raw sequencing data format containing reads and quality scores

**Read:** A short DNA sequence (50-150 bp) from one fragment

**Paired-end:** Sequencing both ends of each fragment (R1 and R2)

**Quality score (Phred):** Confidence in each base call (Q20 = 99%, Q30 = 99.9%)

**Adapter:** Synthetic DNA sequences added during library prep (should be removed)

**Trimming:** Removing low-quality bases and adapters from reads

**Reference genome:** Complete DNA sequence of an organism

**GTF (Gene Transfer Format):** File describing gene locations and structures

**Exon:** Coding region of a gene (kept in mRNA)

**Intron:** Non-coding region (removed during splicing)

**Splice junction:** Where two exons connect in mRNA

**Index (STAR):** Preprocessed genome for fast alignment

**Alignment/Mapping:** Finding where reads came from in the genome

**BAM:** Binary file containing aligned reads

**Count matrix:** Table of read counts (genes × samples)

**De novo assembly:** Reconstructing transcripts without a reference

**Isoform:** Alternative version of a gene's transcript

**Normalization:** Adjusting for sequencing depth differences

**Size factor:** Normalization factor for each sample

**Dispersion:** Measure of biological variability

**log2FoldChange:** Expression difference on log2 scale (log2FC = 1 means 2× difference)

**p-value:** Probability that difference is due to chance

**padj / FDR:** p-value adjusted for multiple testing (False Discovery Rate)

**DEG:** Differentially Expressed Gene (statistically significant)

**Heatmap:** Visual representation of expression patterns

**z-score:** Standardized expression value (mean = 0, SD = 1)

### Decision Tree: Which Pipeline Should I Use?

```
START HERE
    ↓
Do you have a high-quality reference genome?
    ├─ YES → Do you have a GTF annotation?
    │         ├─ YES → Use STAR ✅ (easiest, fastest, most common)
    │         └─ NO  → Use STAR without GTF or Trinity
    │
    └─ NO  → Is your organism closely related to one with a reference?
              ├─ YES → Try STAR with related species genome
              └─ NO  → Use Trinity (de novo assembly)
```

### Typical Analysis Timeline

**STAR pipeline (per sample):**
1. FastQC: 5-15 minutes
2. fastp (if needed): 10-30 minutes
3. STAR indexing: 30-60 minutes (once per genome)
4. STAR alignment: 15-45 minutes per sample
5. Build count matrix: 5 minutes
6. PyDESeq2 analysis: 5-20 minutes

**Total for 12 samples: ~1-2 days (mostly hands-off computer time)**

**Trinity pipeline:**
1. FastQC: 5-15 minutes
2. fastp: 10-30 minutes
3. Trinity assembly: **6-48 hours** (bottleneck!)
4. Quantification: 1-4 hours per sample
5. Build count matrix: 5 minutes
6. DESeq2 analysis: 5-20 minutes

**Total: 2-4 days** (mostly waiting for Trinity)

---

## What to Check After Each Step

### After FastQC:
- ✅ Mean quality score > 28?
- ✅ Adapter content < 5%?
- ✅ Sequence length uniform?
- ❌ If multiple FAILs → Run fastp

### After STAR Alignment:
- ✅ Uniquely mapped reads > 70%?
- ✅ Multi-mapped reads < 10%?
- ✅ Unmapped reads < 15%?
- ❌ If uniquely mapped < 50% → Check if you used the right reference genome!

### After Count Matrix:
- ✅ Do replicates cluster together? (PCA plot)
- ✅ Library sizes similar across samples?
- ✅ Most genes have > 10 reads?
- ❌ If one sample is very different → Possible batch effect or failed library

### After DESeq2:
- ✅ Do you have enough DEGs (typically 100s-1000s)?
- ✅ Are DEGs biologically sensible for your contrast?
- ✅ Do controls (housekeeping genes) show no change?
- ❌ If zero DEGs → Increase FDR cutoff or check if you have biological replicates

---

## Next Steps

Now that you understand the pipeline conceptually:

1. **Decide which path:** STAR or Trinity?
   - If you have a reference genome for your organism → STAR
   - If not → Trinity

2. **Run QC on your data:**
   ```bash
   fastqc raw_data/*.fastq.gz -o qc_reports/
   multiqc qc_reports/ -o qc_reports/
   # Open qc_reports/multiqc_report.html in browser
   ```

3. **If QC looks good, proceed to alignment (STAR) or assembly (Trinity)**

4. **Build count matrix**

5. **Run differential expression analysis**

6. **Interpret your results** (this is where the biology happens!)

---

## Questions to Ask Yourself

As you go through the pipeline:

1. **Do my replicates cluster together?** (PCA plot after normalization)
   - If not → Possible batch effects or sample swaps

2. **Do my DEGs make biological sense?**
   - Example: If comparing leaf vs root, are photosynthesis genes up in leaf?

3. **What's the biological story?**
   - Don't just report "500 DEGs" - what pathways? what processes?

4. **Can I validate key findings?**
   - qPCR for a few important genes
   - Check if results match published literature

---

## Common Pitfalls (Learn from Others' Mistakes!)

1. **Skipping QC** → Garbage in, garbage out
2. **Using wrong reference genome** → Most reads won't map
3. **No biological replicates** → Can't do statistics!
4. **Ignoring batch effects** → Confounds results
5. **Over-interpreting fold changes without considering base expression**
   - 10× change in a gene with 5 reads → Probably noise
   - 1.5× change in a gene with 10,000 reads → Might be real
6. **Not adjusting p-values** → False positives everywhere!
7. **Cherry-picking significant genes** → Confirmation bias

---

## Resources for Learning More

**For STAR:**
- STAR manual: https://github.com/alexdobin/STAR
- Nature Methods paper: Dobin et al. 2013

**For Trinity:**
- Trinity GitHub: https://github.com/trinityrnaseq/trinityrnaseq
- Trinity documentation: https://github.com/trinityrnaseq/trinityrnaseq/wiki

**For DESeq2:**
- DESeq2 vignette (R): https://bioconductor.org/packages/DESeq2/
- PyDESeq2 (Python): https://github.com/owkin/PyDESeq2
- Original paper: Love et al. 2014

**General RNA-seq:**
- RNA-seqlopedia: https://rnaseq.uoregon.edu/
- Harvard Chan Bioinformatics Core training: https://hbctraining.github.io/

---

## Your Specific Pipeline Commands

Based on your scripts, here's what you'd actually run:

### QC:
```bash
# Run FastQC
python -m rna_pipeline.main qc /path/to/fastq_files/ --outdir qc_reports --threads 8

# Review multiqc_report.html
# If bad quality detected, run fastp:
python -m rna_pipeline.main trim sample_R1.fq.gz sample_R2.fq.gz --threads 8
```

### STAR alignment:
```bash
# 1. Build index (once):
python -m rna_pipeline.main star-index \
    --genome genome.fa \
    --gtf annotation.gtf \
    --outdir star_index \
    --readlen 150

# 2. Align each sample:
python -m rna_pipeline.main star-align \
    --genome-index star_index \
    --reads-left sample_R1.fq.gz \
    --reads-right sample_R2.fq.gz \
    --outprefix sample1_ \
    --threads 8

# 3. Build count matrix:
python build_count_matrix.py aligned_reads/ -o count_matrices/

# 4. Differential expression:
python pydeseq2_analysis.py \
    count_matrices/gene_count_matrix.tsv \
    count_matrices/sample_metadata.tsv \
    -o results \
    --contrast-A root \
    --contrast-B leaf \
    --padj 0.05 \
    --lfc 1.0
```

---

**End of Teaching Guide**

Come back to this document whenever you're confused about a step. Each concept builds on the previous one - if something doesn't make sense, re-read the earlier sections!

**Remember:** The computer does the hard computational work, but YOU must understand what it's doing to interpret the results correctly. That's what makes you a scientist, not just a button-pusher!

Good luck with your RNA-seq project! 🧬🔬

---

*Last updated: 2026-01-29*
*Questions? Ask your PI (or me!) in lab meeting.*
