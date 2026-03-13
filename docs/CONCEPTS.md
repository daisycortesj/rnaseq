# RNA-seq Concepts Guide

Background knowledge for understanding the pipeline. Useful if you're new to
RNA-seq or need a refresher on why each step exists.

---

## Table of Contents

1. [FASTQ Files](#1-fastq-files)
2. [Quality Control](#2-quality-control)
3. [Trimming & Filtering](#3-trimming--filtering)
4. [Alignment: STAR vs Trinity](#4-alignment-star-vs-trinity)
5. [Read Counting](#5-read-counting)
6. [Differential Expression (PyDESeq2)](#6-differential-expression-pydeseq2)
7. [BLAST Annotation](#7-blast-annotation)
8. [Domain Analysis (HMMER & PROSITE)](#8-domain-analysis-hmmer--prosite)
9. [Gene Families & Cross-Species Comparison](#9-gene-families--cross-species-comparison)

---

## 1. FASTQ Files

Raw sequencing output. Each read has four lines:

```
@READ_001                    Read name
ACTGACTGACTGACTGACTG         Nucleotide sequence
+                            Separator
IIIIIIIIIIIIIIIIIII          Quality scores (Phred-encoded)
```

Quality scores map to error probability: `I` (Phred ~40) = 99.99% accuracy,
`#` (Phred ~2) = 60% accuracy.

**Paired-end reads** (R1 + R2) sequence both ends of the same RNA fragment,
improving alignment accuracy.

---

## 2. Quality Control

**FastQC** generates per-sample HTML reports covering:
- Per-base quality (should be >Q20 across the read)
- Adapter contamination (spikes at the end = adapters not trimmed)
- GC content distribution (bimodal = possible contamination)
- Sequence duplication (some expected in RNA-seq, extreme = PCR bias)
- Overrepresented sequences (adapters, rRNA)

**MultiQC** aggregates all FastQC reports into a single interactive page for
side-by-side comparison across samples.

### What to look for

| Metric | Good | Concerning |
|--------|------|------------|
| Per-base quality | >Q28 everywhere | Drops below Q20 at 3' end |
| Adapter content | <5% at all positions | Rising curve at 3' end |
| GC content | Single smooth peak | Bimodal or shifted peak |
| Duplication | <50% (RNA-seq) | >70% or one sequence dominates |

---

## 3. Trimming & Filtering

**fastp** handles adapter removal, quality trimming, and read filtering in one
pass. It automatically detects Illumina adapters.

Key operations:
- **Adapter trimming**: removes sequencing adapters from read ends
- **Quality trimming**: trims low-quality bases from the 3' end
- **Length filtering**: discards reads shorter than a threshold (e.g., 36 bp)
- **Complexity filtering**: removes low-complexity reads (e.g., polyA runs)

fastp produces a JSON/HTML report showing how many reads passed, failed, and
why — useful for confirming the trim was effective.

### When to trim

Run FastQC first. If quality is already high and adapter contamination is <1%,
trimming may not change much but is still good practice for reproducibility.

---

## 4. Alignment: STAR vs Trinity

This pipeline primarily uses **STAR** (reference-based alignment), with
Trinity available for de novo assembly when no reference genome exists.

### STAR (Spliced Transcripts Alignment to a Reference)

1. **Genome indexing** (one-time): builds a suffix array from the reference
   genome + GTF annotation so STAR can align quickly.
2. **Alignment**: maps each read to the genome, handling splice junctions
   (intron-spanning reads). Outputs BAM files.

STAR is fast and splice-aware, making it the standard for organisms with a
reference genome.

### Trinity (De Novo Assembly)

Assembles reads into transcripts *without* a reference genome:
1. Breaks reads into k-mers
2. Builds a de Bruijn graph
3. Resolves paths into contigs (assembled transcript fragments)

Use Trinity when:
- No reference genome is available
- You want to discover novel transcripts not in the annotation

This pipeline uses STAR for *D. carota* (reference genome available) and keeps
Trinity scripts for exploratory or supplemental assemblies.

---

## 5. Read Counting

**featureCounts** (from Subread) counts how many aligned reads overlap each
gene defined in the GTF annotation. STAR also outputs per-gene counts in its
`ReadsPerGene.out.tab` files.

**build_count_matrix.py** merges per-sample count files into a single matrix:
rows = genes, columns = samples. This matrix is the input for differential
expression.

### Strandedness

RNA-seq libraries can be unstranded, forward-stranded, or reverse-stranded.
The count tool must match your library prep. STAR's `ReadsPerGene.out.tab`
provides columns for all three; the build script selects the correct one.

---

## 6. Differential Expression (PyDESeq2)

PyDESeq2 is the Python implementation of the DESeq2 method. It tests whether
gene expression differs significantly between conditions (e.g., treatment vs
control).

### Three-step workflow

1. **Step 1 — Run analysis**: fits a negative binomial model per gene,
   estimates size factors (normalization), dispersion, and log2 fold changes.
   Outputs an UNFILTERED results table with all genes.

2. **Step 2 — Filter**: applies statistical cutoffs to keep significant genes.

3. **Step 3 — Plots**: generates MA plots, volcano plots, PCA, and heatmaps.

### Key columns in the results

| Column | Meaning |
|--------|---------|
| `baseMean` | Average normalized count across all samples |
| `log2FoldChange` | Effect size: positive = upregulated, negative = downregulated |
| `pvalue` | Raw p-value from the statistical test |
| `padj` | Adjusted p-value (Benjamini-Hochberg FDR correction) |

### Typical cutoffs

- `padj < 0.05` — statistically significant after multiple-testing correction
- `|log2FoldChange| > 1` — at least 2-fold change in expression

---

## 7. BLAST Annotation

BLAST searches assign functional annotations to genes by finding similar
sequences in curated databases.

### Pipeline flow

1. **Extract CDS** from the GTF + genome → nucleotide FASTA
2. **Translate** nucleotide → protein FASTA
3. **BLASTp** (protein vs protein) or **BLASTx** (nucleotide vs protein)
   against SwissProt or NR databases
4. **Combine** BLAST hits with DESeq2 results to create an annotated table

### Discovery vs Strict modes

- **Discovery** (permissive): e-value 1e-3, keeps more hits for exploratory
  analysis
- **Strict** (high-confidence): e-value 1e-10 + higher identity threshold,
  fewer false positives for publication

---

## 8. Domain Analysis (HMMER & PROSITE)

After BLAST annotation, domain analysis adds an independent layer of
functional evidence.

**HMMER** searches protein sequences against the Pfam database to identify
conserved protein domains (e.g., "p450" domain for CYP genes). Uses profile
Hidden Markov Models.

**PROSITE** (via EMBOSS `patmatmotifs`) scans for short functional motifs and
active sites. Complements HMMER by finding smaller patterns that HMMs may miss.

Together, BLAST + HMMER + PROSITE provide three independent lines of evidence
for gene function.

---

## 9. Gene Families & Cross-Species Comparison

### Gene family extraction

The pipeline extracts specific gene families (CYP450, OMT) from the annotated
results using:
- BLAST hit descriptions containing family keywords
- HMMER domain assignments (e.g., Pfam PF00067 for cytochrome P450)
- GTF keyword search (gene names/descriptions matching "cytochrome P450")

### Cross-species comparison

`compare_species.py` merges DC and DG results to identify:
- Genes with conserved expression patterns across species
- Species-specific differentially expressed genes
- Shared and unique gene family members

### Post-analysis

- **Phylogenetic trees**: MAFFT alignment → FastTree → tree visualization
  overlaid with expression data
- **Genomic clustering**: maps gene family members to chromosomal positions to
  identify tandem duplications and gene clusters
