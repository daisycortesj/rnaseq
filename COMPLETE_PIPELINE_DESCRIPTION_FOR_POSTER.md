# Complete RNA-seq Pipeline Description for Poster
*Template for describing phenylpropanoid pathway enzyme identification*

---

## OVERVIEW

This pipeline identifies candidate enzymes (including cytochrome P450s) involved in the phenylpropanoid pathway leading to myristicin and apiole biosynthesis through comparative transcriptome analysis of tissue types.

---

## STEP-BY-STEP PIPELINE

### **STEP 1: Sample Collection & RNA Sequencing**

**Input:**
- Biological samples from different tissues (root vs. leaf)
- Species: *Daucus carota* ssp. maximus and/or *Daucus glaber*
- Biological replicates per tissue type

**Process:**
- RNA extraction from fresh/frozen tissue
- Library preparation for paired-end RNA sequencing
- Illumina sequencing (typically 50-150 bp reads)

**Output:**
- FASTQ files containing raw sequencing reads (R1 and R2 paired-end)
- Typical yield: 10-50 million read pairs per sample

**Rationale:**
Tissue-specific gene expression profiles reveal where biosynthetic enzymes are active.

---

### **STEP 2: Quality Control (QC)**

**Tools:** FastQC, MultiQC, fastp

**Process:**
- **FastQC**: Assess read quality, adapter content, GC content, sequence duplication
- **MultiQC**: Aggregate QC metrics across all samples for comparison
- **fastp** (if needed): Remove adapters, trim low-quality bases, filter short reads

**Quality Metrics:**
- Target: Mean quality score >28 (Phred score)
- Adapter contamination: <5%
- Per-base quality: Maintained across read length

**Output:**
- QC reports (HTML format)
- Trimmed/cleaned FASTQ files (if needed)

**Rationale:**
Poor quality data produces unreliable results. QC ensures data integrity before analysis.

---

### **STEP 3: Read Alignment to Reference Genome**

**Tool:** STAR (Spliced Transcripts Alignment to a Reference)

**Input:**
- Reference genome (FASTA format)
- Genome annotation (GTF file with gene coordinates)
- Cleaned FASTQ files

**Process:**
1. **Build STAR index** (one-time per genome)
   - Indexes reference genome for rapid alignment
   - Incorporates splice junction information from GTF

2. **Align reads to genome** (per sample)
   - Maps each read to its genomic origin
   - Handles spliced alignments (intron-spanning reads)
   - Outputs alignment statistics

**Output:**
- BAM files (binary alignment files)
- `ReadsPerGene.out.tab` files (raw gene counts)
- Alignment statistics (% mapped, % uniquely mapped)

**Quality Metrics:**
- Target: >70% uniquely mapped reads
- Indicates alignment quality and reference genome appropriateness

**Rationale:**
STAR is splice-aware, essential for RNA-seq where introns are removed from mature transcripts.

---

### **STEP 4: Gene Quantification - Count Matrix**

**Tool:** Custom Python script (`build_count_matrix.py`)

**Input:**
- All `ReadsPerGene.out.tab` files from STAR alignment
- Sample metadata (sample names, tissue types, conditions)

**Process:**
- Aggregate read counts across all samples
- Create matrix: genes (rows) × samples (columns)
- Each cell = number of reads mapped to that gene in that sample

**Output:**
- `gene_count_matrix.tsv`: Complete count table
- `sample_metadata.tsv`: Sample information and experimental design

**Example Count Matrix:**
```
gene_id        DC1L1  DC1L2  DC1R1  DC1R2
LOC108192212   450    523    892    1045
LOC108192214   1023   1156   23     45
```
*(L = leaf, R = root; numbers are read counts)*

**Rationale:**
Raw count data is the foundation for statistical analysis of differential expression.

---

### **STEP 5: Differential Expression Analysis**

**Tool:** PyDESeq2 (Python implementation of DESeq2)

**Input:**
- Gene count matrix
- Sample metadata with experimental design

**Process:**

#### **5.1 Normalization**
- Calculate size factors to account for sequencing depth differences
- Normalizes library sizes across samples

#### **5.2 Dispersion Estimation**
- Estimate gene-wise biological variability
- Models variance-mean relationship across genes
- Shrinks dispersions toward fitted trend

#### **5.3 Statistical Testing**
- Fits negative binomial model for count data
- Tests each gene for differential expression between conditions
- Calculates:
  - **baseMean**: Average normalized expression
  - **log2FoldChange**: Expression difference (log2 scale)
  - **pvalue**: Raw statistical significance
  - **padj**: Adjusted p-value (False Discovery Rate correction)

#### **5.4 Multiple Testing Correction**
- Applies Benjamini-Hochberg correction
- Controls False Discovery Rate (FDR) when testing thousands of genes
- padj < 0.05 means <5% of significant genes are false positives

**Output:**
- `pydeseq2_results_UNFILTERED.tsv`: All genes with statistics
- `pydeseq2_results_FILTERED.tsv`: Significant DEGs only
  - Typical filters: padj < 0.05, |log2FC| > 2 (4-fold change)

**Output Interpretation:**
- **Positive log2FC**: Gene upregulated in condition A (e.g., root)
- **Negative log2FC**: Gene upregulated in condition B (e.g., leaf)
- **log2FC = +2**: 4× higher expression in root vs. leaf
- **padj < 0.05**: Statistically significant

**Rationale:**
Identifies which genes are significantly different between tissues. Genes upregulated in tissues that accumulate myristicin/apiole are prime candidates for biosynthetic enzymes.

---

### **STEP 6: Gene Annotation from Reference**

**Tool:** Custom Python script (`parse_refseq_gtf.py`)

**Input:**
- Reference genome annotation (GTF file)
- List of genes of interest (from PyDESeq2 results)

**Process:**
1. Parse GTF file to extract:
   - Gene IDs
   - Transcript IDs (XM_ accessions)
   - Gene descriptions/functions
   
2. Inner join: Keep only genes in both:
   - Your differentially expressed gene list
   - Reference annotation

3. Create comprehensive annotation table

**Output:**
- `filtered_genes_with_annotations.tsv`: Gene IDs with descriptions
  - Example: LOC108192212 | XM_017458819.2 | cytochrome P450 71A1

**Rationale:**
Links gene IDs to functional descriptions, enabling biological interpretation.

---

### **STEP 7: Sequence Extraction & Translation**

**Tools:** Custom Python scripts (`extract_cds_from_gtf.py`, `run_translate_cds.sbatch`)

**Input:**
- Reference genome (FASTA)
- Reference annotation (GTF)
- Filtered gene IDs

**Process:**

#### **7.1 Extract CDS Sequences**
Extract coding sequences (CDS) for genes of interest from the reference genome using GTF coordinates. FASTA headers include gene ID, transcript ID, and description.

```bash
sbatch scripts/extract_cds__for_blast.sbatch DC
sbatch scripts/extract_cds__for_blast.sbatch DG
```

#### **7.2 Translate CDS to Protein**
CDS nucleotide sequences are translated to amino acid sequences using Biopython, producing protein FASTA files required for BLASTp and HMMER.

```bash
sbatch scripts/run_translate_cds.sbatch DC
sbatch scripts/run_translate_cds.sbatch DG
```

**Output:**
- `all_genes_cds.fasta`: Nucleotide CDS sequences
- `all_genes_protein.fasta`: Translated protein sequences

**Rationale:**
Protein sequences enable more sensitive homology searches (BLASTp) and are required for domain identification (HMMER).

---

### **STEP 8: Functional Annotation - Homology Search & Integration**

**Tools:** BLASTp, custom scripts (`blastp_discoveryfilter.sbatch`, `run_combine_filter.sbatch`)

**Input:**
- Translated protein sequences (`all_genes_protein.fasta`)
- Protein databases (SwissProt, NCBI nr, RefSeq)

**Process:**

#### **8.1 BLASTp** (protein query vs. protein database)
Translated protein sequences are searched against curated protein databases using BLASTp in discovery mode with permissive thresholds to capture distant homologs. The script has built-in resume support -- if a job times out, resubmitting the same command continues from where it left off.

```bash
# DC against SwissProt:
sbatch scripts/blastp_discoveryfilter.sbatch DC swissprot

# DG against SwissProt:
sbatch scripts/blastp_discoveryfilter.sbatch DG swissprot

# Against NCBI nr (larger, slower):
sbatch scripts/blastp_discoveryfilter.sbatch DC nr
sbatch scripts/blastp_discoveryfilter.sbatch DG nr
```

**Discovery mode parameters:** evalue=1e-4, qcov_hsp_perc=40, max_target_seqs=25

#### **8.2 Combine BLAST + DESeq2 Results**
BLAST annotations are merged with PyDESeq2 differential expression statistics into a single table per gene, then filtered for significant candidates using adjustable thresholds.

```bash
# Combine and filter for each species:
sbatch scripts/run_combine_filter.sbatch DC swissprot discovery
sbatch scripts/run_combine_filter.sbatch DG swissprot discovery

# With stricter thresholds (optional):
sbatch scripts/run_combine_filter.sbatch DC swissprot discovery strict
```

**Filter modes:**
- `standard` -- padj < 0.05, |log2FC| > 2.0 (default)
- `strict` -- padj < 0.01, |log2FC| > 2.5 (publication-ready)
- `lenient` -- padj < 0.1, |log2FC| > 1.0 (exploratory)

**Output:**
- `{SP}_{DB}_{MODE}_annotated.tsv`: Full combined table (every gene with expression stats + BLAST annotation)
- `{SP}_{DB}_{MODE}_filtered_{FILTER}.tsv`: Significant candidates only
- `{SP}_{DB}_{MODE}_filtered_{FILTER}_CYP_only.tsv`: CYP candidate subset

**Example:**
```
gene_id         log2FC   padj      blast_description
LOC108192212    +3.4     2.1e-08   cytochrome P450 71A1 [Arabidopsis thaliana]
LOC108197651    -2.8     4.5e-05   caffeic acid O-methyltransferase [Daucus carota]
```

**Rationale:**
Merging BLAST with expression data in a single table links functional annotation to tissue-specific expression, allowing direct identification of differentially expressed enzymes.

---

### **STEP 9: Protein Domain & Motif Analysis**

**Tools:** HMMER (Pfam domains), PROSITE (motifs)

**Input:**
- Protein sequences (`all_genes_protein.fasta`)
- Pfam database (downloaded automatically)

**Process:**

#### **9.1 Download Pfam Database** (one-time)
Downloads and indexes the Pfam-A HMM database required for HMMER scans.

```bash
sbatch scripts/download_protein_databases.sbatch
```

#### **9.2 HMMER - Domain Search**
Protein sequences are scanned against the Pfam database using HMMER hmmscan, which is more sensitive than BLAST for identifying conserved protein domains such as PF00067 (CYP) and PF00891 (OMT).

```bash
# Scan all proteins for each species:
sbatch scripts/run_hmmer.sbatch DC all_genes_protein.fasta
sbatch scripts/run_hmmer.sbatch DG all_genes_protein.fasta
```

**Key domains detected:**
- Cytochrome P450 domain (PF00067) -- confirms CYP family membership
- O-methyltransferase domain (PF00891) -- confirms OMT family membership
- Methyltransferase dimerisation (PF08100)

#### **9.3 PROSITE - Motif Search** (optional)
- P450 heme-binding motif (FxxGxxxCxG)
- SAM-binding motif (methyltransferases)
- Catalytic residues

**Output:**
- `all_genes_protein_pfam_domains.txt`: Domain table (gene ID, domain name, Pfam accession, E-value, coordinates)
- `all_genes_protein_pfam_domains.out`: Human-readable domain report

**Rationale:**
HMMER provides independent evidence for enzyme family classification beyond BLAST homology. A gene with both a BLAST hit to "cytochrome P450" and a PF00067 Pfam domain has stronger evidence than BLAST alone.

---

### **STEP 10: Integration & Candidate Prioritization**

**Tools:** Custom scripts (`compare_species.py`, `extract_gene_families_combined.py`)

**Process:**

#### **10.1 Cross-Species Comparison (DC vs DG)**
Annotated results from both species are merged into a comparison table categorizing each gene by its expression pattern: shared same-direction, shared opposite-direction, DC-only, or DG-only.

```bash
sbatch scripts/run_compare_species.sbatch DC DG swissprot discovery

# With custom significance cutoffs:
sbatch scripts/run_compare_species.sbatch DC DG swissprot discovery 0.01 2.5
```

**Output:**
- `DC_vs_DG_comparison.tsv`: Full comparison table
- `DC_vs_DG_comparison_shared_same_direction.tsv`: Genes DE in both species, same direction
- `DC_vs_DG_comparison_DC_only.tsv` / `DG_only.tsv`: Species-unique DE genes

#### **10.2 Gene Family Extraction (Multi-Evidence)**
CYP and OMT genes are identified using both BLAST annotation patterns and HMMER Pfam domains; each gene receives a confidence score based on how many evidence sources agree.

**CYP BLAST patterns** (tightened to avoid false positives):
- `cytochrome P450`, `CYP` + digit, `P450` (exact word)
- Excludes: NADPH-cytochrome P450 reductase, cytochrome P450 reductase, cytochrome b5, cytochrome c
- The broad `monooxygenase` pattern was removed -- it matched many non-CYP enzymes (flavin-dependent, copper-containing, etc.)
- HMMER Pfam domain PF00067 provides independent confirmation

```bash
# BLAST + HMMER combined evidence for both species:
sbatch scripts/run_gene_family_combined.sbatch DC DG swissprot discovery full

# Single species (BLAST only):
sbatch scripts/run_gene_family_heatmap.sbatch DC swissprot discovery
```

**Evidence scoring:**
- **BLAST + HMMER** (both agree) = highest confidence
- **BLAST only** or **HMMER only** = moderate confidence
- Genes included if ANY source flags them

#### **10.3 Prioritization Criteria**

Integrate multiple lines of evidence:

1. **Differential Expression** -- Upregulated in metabolite-producing tissues (padj < 0.05, |log2FC| > 2)
2. **Gene Family Classification** -- CYP450 or OMT family confirmed by BLAST patterns
3. **Domain Confirmation** -- Pfam domain detected by HMMER (PF00067 for CYP, PF00891 for OMT)
4. **Cross-Species Conservation** -- Same direction DE in both DC and DG

**Confidence tiers:**
- **Tier 1**: All 4 criteria met (highest priority)
- **Tier 2**: 3/4 criteria met
- **Tier 3**: 2/4 criteria met

**Output:**
- `DC_DG_CYP_verified.tsv` / `DC_DG_OMT_verified.tsv`: Verified gene family members with confidence scores
- `DC_DG_CYP_OMT_combined.tsv`: All verified CYP + OMT genes
- `DC_DG_CYP_OMT_shared_same_direction.tsv`: Conserved candidates (both species, same direction)
- `DC_DG_gene_family_summary.txt`: Summary statistics

---

### **STEP 11: Visualization & Interpretation**

**Tools:** Custom scripts (`pydeseq2_generate_plots.py`, `generate_family_heatmap.py`)

**Commands:**
```bash
# Single-species plots (MA, volcano, PCA, correlation, CYP/OMT heatmaps):
sbatch scripts/run_pydeseq2_step3_plots.sbatch DC combined_annotated
sbatch scripts/run_pydeseq2_step3_plots.sbatch DG combined_annotated

# Combined two-species heatmap (DC + DG on one plot):
sbatch scripts/run_pydeseq2_step3_plots.sbatch DC combined_annotated DG
```

**Heatmap DE filter (on by default):** Only differentially expressed genes
appear in CYP/OMT heatmaps. Defaults: `PADJ_CUTOFF=0.05`, `LFC_CUTOFF=2.0`.
A gene must match a CYP/OMT BLAST pattern AND pass padj < 0.05, |log2FC| > 2
to be included. Adjust these inside the sbatch script or via CLI:

```bash
# Inside run_pydeseq2_step3_plots.sbatch:
PADJ_CUTOFF=0.05   # padj threshold for heatmap genes
LFC_CUTOFF=2.0     # |log2FC| threshold for heatmap genes
TOP_N=50            # optional: further limit to top 50 per family
```

Or pass directly to the Python script:
```bash
python3 pydeseq2_generate_plots.py results.tsv \
  --count-matrix counts.tsv --metadata meta.tsv \
  --padj-cutoff 0.05 --lfc-cutoff 2.0 --top-n 50
```

The script auto-detects the GTF annotation file from the reference directory. HMMER domain information is saved to the gene list TSV files when the Pfam scan output exists.

**Visualizations Created:**

#### **11.1 MA Plot**
Mean expression vs. log2 fold change with three-color coding and gene count annotations for plant transcriptomics.

- X-axis: Mean of normalized counts (log scale)
- Y-axis: Log2 fold change
- **Red points**: Significantly up in Root (padj < 0.05, |log2FC| > 2)
- **Blue points**: Significantly up in Leaf
- **Gray points**: Not significant
- Horizontal dashed lines at fold-change thresholds
- Gene counts for up/down annotated on the plot

#### **11.2 Enhanced Volcano Plot**
Four-category volcano plot with gene count annotations per category.

- **Gray (NS):** Neither threshold met
- **Green (|Log2FC| > 2):** Passes fold-change but not significance
- **Blue (padj < 0.05):** Passes significance but not fold-change
- **Red (both):** Passes both thresholds -- the DE candidates
- Title: "Differential Expression: Root vs Leaf"
- Bottom annotation: total genes, counts up/down per tissue
- Top 10 most significant genes labeled by ID
- Vertical and horizontal dashed lines at cutoff thresholds

#### **11.3 PCA Plot**
PCA of the top 500 most variable genes with sample labels and 95% confidence ellipses.

- PC1 vs PC2 scatter plot with variance explained (%) on each axis
- Samples labeled with short names (DC1L, DC2L, DC1R, DC2R, etc.)
- 95% confidence ellipses drawn around Root and Leaf groups
- Samples colored by tissue: orange = Root, green = Leaf
- Uses SVD on normalized, log2-transformed counts
- Reveals batch effects, outliers, and condition-level separation

#### **11.4 Sample Correlation Heatmap**
Sample-to-sample Euclidean distance matrix with short sample names (DC1L, DC1R, etc.).

- Blue gradient colormap (darker = more similar)
- Condition color bars on both rows and columns
- Dendrograms group similar samples together
- Short display names for clean presentation

#### **11.5 Single-Species Gene Family Heatmaps**
Publication-style heatmaps matching the reference figure design:

- **Gene locus IDs as row labels** on the left of the heatmap
- Genes **grouped by subfamily** with **bracket annotations** on the far left (CYP71D, CYP81, CYP86, CYP72, CYP85, CYP719, etc.)
- Subfamilies parsed automatically from BLAST hit descriptions using regex
- **Short sample names** below: Leaf first (DC1L, DC2L, DC3L), then Root (DC1R, DC2R, DC3R)
- Centered log2 normalized expression (RdBu_r scale, red = high, blue = low)
- Column color bar: green = Leaf, brown = Root
- Separate heatmaps for CYP and OMT gene families

#### **11.6 Combined Two-Species Heatmap**
Combined heatmaps with both species side by side, same design.

- **Gene locus IDs as row labels**, grouped by subfamily with bracket annotations
- **Short sample names**: DC1L, DC2L, DC1R, DC2R, DG1L, DG2L, DG1R, DG2R
- Two stacked column color bars: **Species** (top) + **Tissue** (bottom)
- Each species normalized independently before merging (own size factors)
- Separate combined heatmaps for CYP and OMT families

**Output files:**
- `ma_plot.pdf` -- MA plot (red/blue up/down, gene counts, threshold caption)
- `volcano_plot.pdf` -- Enhanced Volcano plot (4-color, boxed labels with connectors, threshold caption)
- `pca_plot.pdf` -- PCA of top 500 variable genes (sample labels, 95% confidence ellipses)
- `sample_correlation_heatmap.pdf` -- Sample distance heatmap (short names: DC1L, DC1R, etc.)
- `cyp_heatmap.pdf`, `omt_heatmap.pdf` -- Subfamily-bracketed heatmaps with gene locus labels
- `cyp_heatmap_combined.pdf`, `omt_heatmap_combined.pdf` -- Both species combined
- `cyp_gene_list.tsv`, `omt_gene_list.tsv` -- Detected gene family members

**Rationale:**
The five-panel visualization set provides a comprehensive view of the data: PCA confirms sample grouping, the correlation heatmap validates replicate consistency, the enhanced volcano and MA plots summarize genome-wide differential expression, and the gene family heatmaps reveal tissue-specific patterns for CYP and OMT candidates. The combined two-species heatmap is the key figure for this project, revealing which candidates have conserved expression in both *D. carota* and *D. glaber*.

---

## KEY RESULTS FOR POSTER

### **What This Pipeline Reveals:**

1. **Number of Differentially Expressed Genes (DEGs)**
   - Total DEGs identified (padj < 0.05)
   - Tissue-upregulated DEGs (root vs. leaf)

2. **Candidate Enzyme Identification**
   - Number of CYP450 genes identified
   - Number with expression patterns matching metabolite profiles
   - Number with homology to phenylpropanoid enzymes

3. **Pathway Reconstruction**
   - Candidate genes for each biosynthetic step:
     - Hydroxylation (CYPs)
     - Methylation (OMTs)
     - Reduction/oxidation
   - Co-expressed gene clusters

4. **Validation Evidence**
   - Multi-evidence support for top candidates
   - Tissue-specificity correlation with metabolite accumulation

---

## BIOLOGICAL INTERPRETATION

### **Why This Matters:**

**The Challenge:**
- Phenylpropanoid pathways involve sequential enzymatic steps
- Many enzymes (especially CYPs) are poorly characterized
- Genome contains thousands of potential enzyme-coding genes

**The Approach:**
- **Differential expression** narrows thousands of genes to hundreds of candidates
- **Co-expression** with known pathway genes suggests functional involvement
- **Tissue-specific expression** links genes to metabolite-producing tissues
- **Homology & domain analysis** provides functional classification

**The Impact:**
- Identifies specific CYP450s and other enzymes responsible for myristicin/apiole biosynthesis
- Provides targets for:
  - Functional validation (enzyme assays, transgenic studies)
  - Metabolic engineering
  - Understanding specialized metabolism evolution
- Contributes to knowledge of plant chemical defense mechanisms

---

## STATISTICS SUMMARY FOR POSTER

**Data Volume:**
- XX samples sequenced
- XX million reads per sample
- XX,XXX genes quantified

**Quality Metrics:**
- >70% uniquely mapped reads
- XX biological replicates per condition

**Differential Expression:**
- XX,XXX genes tested
- XXX significantly differentially expressed (padj < 0.05)
- XX upregulated in root, XX upregulated in leaf

**Candidate Enzymes:**
- XX CYP450 genes identified
- XX with tissue-specific expression
- XX with homology to phenylpropanoid enzymes
- XX Tier 1 candidates (all criteria met)

---

## SOFTWARE & RESOURCES

**Analysis Tools:**
- FastQC / MultiQC -- Quality control
- fastp -- Read trimming and adapter removal
- STAR -- Splice-aware read alignment
- PyDESeq2 -- Differential expression (Python implementation of DESeq2)
- NCBI BLAST+ -- Protein homology search (BLASTp)
- HMMER -- Profile HMM domain identification (hmmscan)
- Biopython -- CDS to protein translation
- seaborn / matplotlib -- Heatmap and plot generation (Enhanced Volcano, PCA, correlation)
- scipy -- Sample distance computation

**Databases:**
- Reference genome: *Daucus carota* ssp. maximus (GCF_001625215.2 DH1 v3.0)
- NCBI RefSeq GTF annotation
- SwissProt (curated protein database)
- NCBI nr (non-redundant protein database)
- Pfam-A (protein family HMM profiles)

**Computing:**
- Analysis performed on Virginia Tech ARC TinkerCliffs HPC cluster
- SLURM job scheduler for all compute-intensive steps
- Custom Python/Bash pipeline scripts (available on GitHub)

---

## FLOWCHART FOR POSTER

```
 Raw RNA-seq Reads (DC + DG)
              ↓
     Quality Control (FastQC/fastp)
              ↓
     Alignment (STAR → counts)
              ↓
     Count Matrix (genes × samples)
              ↓
  Differential Expression (PyDESeq2)
              ↓
  ┌───────────┼───────────┐
  ↓           ↓           ↓
GTF       CDS Extract   Count Matrix
Annotation  + Translate  (normalized)
  ↓           ↓           ↓
  ↓    ┌──────┴──────┐   ├── PCA Plot
  ↓    ↓             ↓   └── Sample Correlation
  ↓  BLASTp       HMMER
  ↓  (swissprot)  (Pfam)
  ↓    ↓             ↓
  └──→ Combine BLAST  ←──┘
       + DESeq2
          ↓
  ┌───────┴───────┐
  ↓               ↓
Species        Gene Family
Comparison     Extraction
(DC vs DG)    (BLAST+HMMER)
  ↓               ↓
  └───────┬───────┘
          ↓
  Publication Plots:
  - Enhanced Volcano (4-color, boxed labels)
  - MA Plot (red/blue up/down)
  - PCA (top 500 variable genes, 95% ellipses)
  - Sample Correlation Heatmap
  - Gene Family Heatmaps (subfamily sidebar)
  - Combined Heatmaps (DC+DG)
          ↓
  Candidate Enzyme List
   (CYPs, OMTs, ranked)
```

---

## ONE-SENTENCE PIPELINE SUMMARY

"This integrated RNA-seq approach combines differential expression analysis, homology searching, and protein domain identification to systematically discover tissue-specific biosynthetic enzymes, including cytochrome P450s, responsible for myristicin and apiole production in the phenylpropanoid pathway."

---

## FILL IN YOUR SPECIFIC DETAILS:

- [ ] Species name(s)
- [ ] Number of samples
- [ ] Tissue types compared
- [ ] Number of replicates
- [ ] Total DEGs found
- [ ] Number of CYP candidates
- [ ] Top candidate genes (with IDs)
- [ ] Correlation with metabolite data (if available)

---

**END OF TEMPLATE**

*Use this as a starting point for your poster. Adapt to your specific results and poster format. Focus on the biological story, not just the technical pipeline.*
