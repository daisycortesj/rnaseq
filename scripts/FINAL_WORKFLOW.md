# Final Updated Workflow - All Scripts Synchronized

## ✅ All Scripts Now Work Together!

I've updated all three scripts to work with the new workflow that includes:
- ✅ Automatic transcript ID (XM) lookup
- ✅ Enhanced FASTA headers with descriptions
- ✅ Consistent use of filtered gene lists

---

## 📚 Updated Scripts

### **1. `parse_refseq_gtf.py`** (Smart Parser)

**What changed:**
- ✅ Now does automatic transcript ID lookup (no more `--feature` flag needed)
- ✅ Three-pass parsing: genes → transcript lookup → merge
- ✅ One row per gene WITH transcript IDs

**Usage:**
```bash
python scripts/parse_refseq_gtf.py \
  --gtf 04_reference/dc_genomic.gtf \
  --out genes_with_XM.tsv
```

**Output:**
```
local_number	gene_id	        transcript_id	    description
1	        LOC108192212	XM_017458819.2	    MDIS1-interacting receptor...
```

---

### **2. `extract_cds_from_gtf.py`** (Enhanced Headers)

**What changed:**
- ✅ Now includes transcript_id and description in FASTA headers
- ✅ Format: `>gene_id|transcript_id description`
- ✅ Added note about using simple gene ID lists (not TSV)

**Usage:**
```bash
python scripts/extract_cds_from_gtf.py \
  04_reference/dc_genomic.gtf \
  04_reference/genome.fna \
  filtered_gene_ids.txt \
  output.fasta
```

**Output FASTA:**
```fasta
>LOC108192212|XM_017458819.2 MDIS1-interacting receptor like kinase 2
ATGGCTAGCTAGC...
```

---

### **3. `prepare_blast_input.sh`** (Master Script - Bash)

**What changed:**
- ✅ Uses smart parser (automatic XM lookup)
- ✅ Extracts filtered gene IDs from TSV (Step 2.5)
- ✅ Uses filtered gene IDs for sequence extraction
- ✅ Updated output messages

**Usage:**
```bash
bash scripts/prepare_blast_input.sh DC
```

**Does automatically:**
1. Parse GTF with smart parser (gets XM)
2. Inner join (filter for your genes)
3. Extract gene IDs from filtered TSV
4. Extract sequences using filtered IDs

---

### **4. `extract_cds__for_blast.sbatch`** (Master Script - SLURM)

**What changed:**
- ✅ Now runs complete 3-step workflow
- ✅ Uses smart parser
- ✅ Extracts filtered gene IDs
- ✅ Uses filtered IDs for sequences

**Usage:**
```bash
sbatch scripts/extract_cds__for_blast.sbatch DC
```

**Does automatically:**
Same as bash version, but submits to SLURM!

---

## 🔄 Complete Workflow

### **Manual (Step by Step):**
```bash
cd /projects/tholl_lab_1/daisy_analysis

# Step 1: Parse GTF (smart parser - auto XM lookup)
python 05_rnaseq-code/scripts/parse_refseq_gtf.py \
  --gtf 04_reference/dc_genomic.gtf \
  --out 06_analysis/blast_input_DC/genes_with_XM.tsv

# Step 2: Inner join (filter)
python 05_rnaseq-code/scripts/join_gtf_with_gene_ids.py \
  --parsed 06_analysis/blast_input_DC/genes_with_XM.tsv \
  --gene-ids 06_analysis/blast_input_DC/all_gene_ids.txt \
  --out 06_analysis/blast_input_DC/filtered_complete.tsv

# Step 2.5: Extract gene IDs from filtered TSV
tail -n +2 06_analysis/blast_input_DC/filtered_complete.tsv | cut -f2 > 06_analysis/blast_input_DC/filtered_gene_ids.txt

# Step 3: Extract sequences (using filtered IDs)
python 05_rnaseq-code/scripts/extract_cds_from_gtf.py \
  04_reference/dc_genomic.gtf \
  04_reference/GCF_001625215.2_DH1_v3.0_genomic.fna \
  06_analysis/blast_input_DC/filtered_gene_ids.txt \
  06_analysis/blast_input_DC/all_genes_cds.fasta
```

---

### **Automatic (One Command - Bash):**
```bash
cd /projects/tholl_lab_1/daisy_analysis

bash 05_rnaseq-code/scripts/prepare_blast_input.sh \
  04_reference/dc_genomic.gtf \
  04_reference/GCF_001625215.2_DH1_v3.0_genomic.fna \
  06_analysis/blast_input_DC/all_gene_ids.txt \
  06_analysis/blast_input_DC
```

---

### **Automatic (One Command - SLURM):**
```bash
cd /projects/tholl_lab_1/daisy_analysis

sbatch 05_rnaseq-code/scripts/extract_cds__for_blast.sbatch DC
```

---

## 📊 What You Get

### **File 1: `genes_with_XM.tsv`**
All genes from GTF with transcript IDs
```
local_number	gene_id	        transcript_id	    description
1	        LOC108192212	XM_017458819.2	    MDIS1-interacting...
...
39102	        LOC999999999	XM_999999999.1	    zinc finger...
```

### **File 2: `filtered_complete.tsv`** ⭐
Only YOUR genes with complete annotations
```
local_number	gene_id	        transcript_id	    description
1	        LOC108192212	XM_017458819.2	    MDIS1-interacting...
...
34420	        LOC999999999	XM_999999999.1	    protein kinase...
```

### **File 3: `filtered_gene_ids.txt`**
Simple gene ID list extracted from file 2
```
LOC108192212
LOC108192214
LOC108192216
...
```

### **File 4: `all_genes_cds.fasta`** ⭐
Sequences with enhanced headers
```fasta
>LOC108192212|XM_017458819.2 MDIS1-interacting receptor like kinase 2
ATGGCTAGCTAGC...
>LOC108192214|XM_017458821.1 heat shock protein 70
ATGCTAGCT...
```

---

## 🎯 Key Improvements

### **1. Smart Parser (Automatic XM Lookup)**
- ✅ No more manual transcript extraction
- ✅ One row per gene WITH transcript ID
- ✅ Meets PI requirements automatically

### **2. Enhanced FASTA Headers**
- ✅ Includes transcript IDs (XM numbers)
- ✅ Includes descriptions (gene functions)
- ✅ Format: `>gene_id|transcript_id description`

### **3. Consistent Filtering**
- ✅ Step 3 now uses filtered gene IDs from Step 2
- ✅ All output files have SAME genes
- ✅ Perfect alignment between TSV and FASTA

### **4. Updated Master Scripts**
- ✅ Bash script (`prepare_blast_input.sh`) - Updated
- ✅ SLURM script (`extract_cds__for_blast.sbatch`) - Updated
- ✅ Both run complete workflow automatically

---

## 📋 File Relationships

```
all_gene_ids.txt (34,408 genes - original list)
        ↓ [Step 1: Smart parser]
genes_with_XM.tsv (39,102 genes - all from GTF, with XM)
        ↓ [Step 2: Inner join]
filtered_complete.tsv (34,420 genes - only yours, with XM) ⭐
        ↓ [Step 2.5: Extract IDs]
filtered_gene_ids.txt (34,420 gene IDs - simple list)
        ↓ [Step 3: Extract sequences]
all_genes_cds.fasta (31,479 sequences - protein-coding only) ⭐
```

**All perfectly aligned!** ✅

---

## 🚀 How to Use Now

### **On HPC with SLURM (Easiest):**
```bash
cd /projects/tholl_lab_1/daisy_analysis
sbatch 05_rnaseq-code/scripts/extract_cds__for_blast.sbatch DC
```

### **On HPC without SLURM:**
```bash
cd /projects/tholl_lab_1/daisy_analysis
bash 05_rnaseq-code/scripts/prepare_blast_input.sh \
  04_reference/dc_genomic.gtf \
  04_reference/GCF_001625215.2_DH1_v3.0_genomic.fna \
  06_analysis/blast_input_DC/all_gene_ids.txt \
  06_analysis/blast_input_DC
```

### **Manual (For Learning/Testing):**
Follow the manual commands in the "Manual (Step by Step)" section above.

---

## ✅ What Changed Summary

| Script | What Changed |
|--------|--------------|
| `parse_refseq_gtf.py` | Now smart parser - auto-lookups XM numbers |
| `extract_cds_from_gtf.py` | Enhanced FASTA headers with XM + description |
| `prepare_blast_input.sh` | Extracts filtered gene IDs, uses them for Step 3 |
| `extract_cds__for_blast.sbatch` | Runs complete workflow, uses filtered IDs |

**All scripts now work together perfectly!** 🎉

---

## 🎓 Your Insights Improved the Workflow!

**Questions you asked that led to improvements:**
1. ✅ "Why manual two-step?" → Created smart parser
2. ✅ "Can we include description?" → Enhanced FASTA headers
3. ✅ "Why so many missing?" → Realized gene ID format check needed
4. ✅ "Is this using filtered list?" → Fixed workflow consistency

**Excellent critical thinking!** Your questions made this workflow much better! 🌟

---

## Part 2: Downstream Analysis (BLAST, HMMER, Gene Families)

Everything above prepares the BLAST input (CDS FASTA). The steps below run
the actual analysis pipeline: translate to protein, BLAST, combine with
DESeq2, HMMER domain scan, species comparison, and gene family extraction.

All commands assume:
```bash
cd /projects/tholl_lab_1/daisy_analysis
```

---

### **Step 4: Translate CDS to Protein**

Converts nucleotide CDS to amino acid sequences for BLASTp and HMMER.

```bash
# Single species:
sbatch 05_rnaseq-code/scripts/run_translate_cds.sbatch DC
sbatch 05_rnaseq-code/scripts/run_translate_cds.sbatch DG

# Or all species at once:
sbatch 05_rnaseq-code/scripts/run_translate_cds.sbatch all
```

**Input:**  `06_analysis/blast_input_{SP}/all_genes_cds.fasta`
**Output:** `06_analysis/blast_input_{SP}/all_genes_protein.fasta`

---

### **Step 5: Run BLASTp (Discovery Mode)**

Runs BLASTp against a protein database. Has built-in resume support -- if the
job times out, resubmit the same command and it picks up where it left off.

```bash
# DC against swissprot:
sbatch 05_rnaseq-code/scripts/blastp_discoveryfilter.sbatch DC swissprot

# DG against swissprot:
sbatch 05_rnaseq-code/scripts/blastp_discoveryfilter.sbatch DG swissprot

# Other databases (nr, refseq_protein, DCdb):
sbatch 05_rnaseq-code/scripts/blastp_discoveryfilter.sbatch DC nr
sbatch 05_rnaseq-code/scripts/blastp_discoveryfilter.sbatch DG nr
```

**Parameters (discovery mode):** evalue=1e-4, qcov_hsp_perc=40, max_target_seqs=25

**Input:**  `06_analysis/blast_input_{SP}/all_genes_protein.fasta`
**Output:** `06_analysis/blastp_{SP}/blastp_{SP}_{DB}_discovery.tsv`

---

### **Step 6: Combine BLAST + DESeq2 Results**

Merges BLAST annotations with PyDESeq2 differential expression statistics into
a single table per gene, then filters for significant candidates.

```bash
# DC + swissprot discovery:
sbatch 05_rnaseq-code/scripts/run_combine_filter.sbatch DC swissprot discovery

# DG + swissprot discovery:
sbatch 05_rnaseq-code/scripts/run_combine_filter.sbatch DG swissprot discovery

# With stricter filters (optional 4th argument):
#   standard  - padj<0.05, |log2FC|>2.0 (default)
#   strict    - padj<0.01, |log2FC|>2.5
#   lenient   - padj<0.1,  |log2FC|>1.0
#   blast_only - any DE, but must have BLAST hit
sbatch 05_rnaseq-code/scripts/run_combine_filter.sbatch DC swissprot discovery strict
```

**Input:**
- `06_analysis/pydeseq2_{SP}/` (DESeq2 step 1 results)
- `06_analysis/blastp_{SP}/blastp_{SP}_{DB}_discovery.tsv`

**Output:** `06_analysis/combined_{SP}/`
- `{SP}_{DB}_{MODE}_annotated.tsv` (full combined table)
- `{SP}_{DB}_{MODE}_filtered_{FILTER}.tsv` (significant candidates)
- `{SP}_{DB}_{MODE}_filtered_{FILTER}_CYP_only.tsv` (CYP genes)

---

### **Step 7: Run HMMER Domain Scan**

Scans protein sequences against the Pfam database to identify conserved
protein domains. More sensitive and specific than BLAST for domain detection.

Requires Pfam database (download first if not present):
```bash
sbatch 05_rnaseq-code/scripts/download_protein_databases.sbatch
```

Run HMMER:
```bash
# All proteins for DC:
sbatch 05_rnaseq-code/scripts/run_hmmer.sbatch DC all_genes_protein.fasta

# All proteins for DG:
sbatch 05_rnaseq-code/scripts/run_hmmer.sbatch DG all_genes_protein.fasta
```

**Input:**  `06_analysis/blast_input_{SP}/all_genes_protein.fasta`
**Output:** `06_analysis/hmmer_{SP}/`
- `all_genes_protein_pfam_domains.txt` (domain table)
- `all_genes_protein_pfam_domains.out` (human-readable)

---

### **Step 8: Compare Species (DC vs DG)**

Merges annotated results from both species into a comparison table showing
shared/unique differential expression. Both species must use the same
reference genome.

```bash
sbatch 05_rnaseq-code/scripts/run_compare_species.sbatch DC DG swissprot discovery

# With custom significance cutoffs:
sbatch 05_rnaseq-code/scripts/run_compare_species.sbatch DC DG swissprot discovery 0.01 2.5
```

**Input:**
- `06_analysis/combined_DC/DC_swissprot_discovery_annotated.tsv`
- `06_analysis/combined_DG/DG_swissprot_discovery_annotated.tsv`

**Output:** `06_analysis/comparison_DC_vs_DG/`
- `DC_vs_DG_{DB}_{MODE}_comparison.tsv` (full table)
- `DC_vs_DG_{DB}_{MODE}_comparison_shared_same_direction.tsv`
- `DC_vs_DG_{DB}_{MODE}_comparison_DC_only.tsv`
- `DC_vs_DG_{DB}_{MODE}_comparison_DG_only.tsv`

---

### **Step 9: Extract Gene Families (CYP + OMT)**

Searches BLAST annotations and HMMER Pfam domains to identify CYP (cytochrome
P450) and OMT (O-methyltransferase) gene families. Genes are scored by how
many evidence sources agree. HMMER is optional.

**CYP pattern specificity:** Only three precise BLAST patterns are used
(`cytochrome P450`, `CYP` + digit, `P450`). The broad `monooxygenase` pattern
was removed to avoid false positives from non-CYP monooxygenases.
`cytochrome c` is also now excluded alongside cytochrome b5 and reductases.

#### Single species:
```bash
# BLAST only:
sbatch 05_rnaseq-code/scripts/run_gene_family_heatmap.sbatch DC swissprot discovery

# BLAST + HMMER (full evidence):
sbatch 05_rnaseq-code/scripts/run_gene_family_heatmap.sbatch DC swissprot discovery full
sbatch 05_rnaseq-code/scripts/run_gene_family_heatmap.sbatch DG swissprot discovery full
```

**Output:** `06_analysis/gene_families_{SP}/`
- `{SP}_CYP_verified.tsv`, `{SP}_OMT_verified.tsv`
- `{SP}_CYP_OMT_combined.tsv`
- `{SP}_gene_family_summary.txt`
- `cyp_heatmap.pdf`, `omt_heatmap.pdf`, `cyp_omt_combined_heatmap.pdf`

#### Combined two-species:
```bash
# BLAST only:
sbatch 05_rnaseq-code/scripts/run_gene_family_combined.sbatch DC DG swissprot discovery

# BLAST + HMMER (full evidence):
sbatch 05_rnaseq-code/scripts/run_gene_family_combined.sbatch DC DG swissprot discovery full
```

**Output:** `06_analysis/gene_families_DC_DG/`
- `DC_DG_CYP_verified.tsv`, `DC_DG_OMT_verified.tsv`
- `DC_DG_CYP_OMT_combined.tsv`
- `DC_DG_CYP_OMT_shared_same_direction.tsv`
- `DC_DG_CYP_OMT_DC_only.tsv`, `DC_DG_CYP_OMT_DG_only.tsv`
- `DC_DG_gene_family_summary.txt`
- `cyp_heatmap_DC.pdf`, `cyp_heatmap_DG.pdf`
- `omt_heatmap_DC.pdf`, `omt_heatmap_DG.pdf`

---

### **Step 10: Publication-Quality Plots (MA, Volcano, PCA, Heatmaps)**

Generates a complete multi-panel figure set from the combined annotated results
and count matrix. Includes enhanced volcano (4-color, boxed labels), MA plot
(red/blue up/down), PCA (95% confidence ellipses), sample correlation heatmap,
and gene family heatmaps with CYP/OMT subfamily color sidebars.

```bash
# Single-species (all 6 plot types):
sbatch 05_rnaseq-code/scripts/run_pydeseq2_step3_plots.sbatch DC combined_annotated
sbatch 05_rnaseq-code/scripts/run_pydeseq2_step3_plots.sbatch DG combined_annotated

# Combined two-species heatmap (adds DC+DG heatmaps):
sbatch 05_rnaseq-code/scripts/run_pydeseq2_step3_plots.sbatch DC combined_annotated DG
```

**Heatmap DE filter (default ON):** Only differentially expressed genes appear
in heatmaps. The defaults are `PADJ_CUTOFF=0.05` and `LFC_CUTOFF=2.0` inside
the sbatch script. Genes must pass padj < 0.05 AND |log2FC| > 2 to be plotted,
even if their BLAST description matches a CYP/OMT pattern.

**Limit heatmap size:** Set `TOP_N=50` inside the sbatch script (or pass `--top-n 50` to the Python scripts) to further limit to the top N genes per family, ranked by adjusted p-value then |log2FC|.

The script auto-detects GTF and HMMER files from the standard directory layout:
- **GTF** (if found): used for biotype annotation in gene list TSVs
- **HMMER** (if found): used for gene family detection (combined with BLAST evidence) and domain info saved to gene list TSVs

**Gene family detection** uses both BLAST descriptions and HMMER Pfam domains.
A gene is included if EITHER source identifies it as CYP/OMT (union of both).
HMMER catches genes whose BLAST hit is vague (e.g. "hypothetical protein") but
whose protein has a clear CYP domain (PF00067) or OMT domain (PF00891).

**Plots generated:**

| Output file | Description |
|-------------|-------------|
| `ma_plot.pdf` | MA plot -- red/blue up/down coloring with gene counts, log-scale mean expression |
| `volcano_plot.pdf` | Enhanced Volcano -- 4 categories with per-category gene counts, Root vs Leaf title, top 10 labeled |
| `pca_plot.pdf` | PCA of top 500 variable genes, DC1L/DC1R sample labels, 95% confidence ellipses |
| `sample_correlation_heatmap.pdf` | Sample distance heatmap with short names (DC1L, DC1R, etc.) |
| `cyp_heatmap.pdf` | CYP heatmap -- gene locus labels, subfamily bracket annotations (CYP71D, CYP81, etc.), Leaf then Root |
| `omt_heatmap.pdf` | OMT heatmap -- same style with OMT subfamily brackets (COMT, CCoAOMT, etc.) |
| `cyp_heatmap_combined.pdf` | Combined DC+DG CYP heatmap with locus labels + subfamily brackets |
| `omt_heatmap_combined.pdf` | Combined DC+DG OMT heatmap (same layout) |
| `cyp_gene_list.tsv` | Detected CYP genes with expression stats |
| `omt_gene_list.tsv` | Detected OMT genes with expression stats |

**Input:**
- `06_analysis/combined_{SP}/{SP}_swissprot_discovery_annotated.tsv`
- `03_count_tables/{FOLDER}/gene_count_matrix.tsv`
- `03_count_tables/{FOLDER}/sample_metadata.tsv`
- `04_reference/{sp}_genomic.gtf` (auto-detected, optional)
- `06_analysis/hmmer_{SP}/all_genes_protein_pfam_domains.txt` (auto-detected, optional)

**Output:** `06_analysis/pydeseq2_{SP}_step3_plots_annotated/`

---

## Complete Pipeline at a Glance

```
all_gene_ids.txt
    |  [Step 1-3: Parse GTF, join, extract CDS]
    v
all_genes_cds.fasta
    |  [Step 4: Translate]
    v
all_genes_protein.fasta ------+--------------------+
    |                         |                    |
    |  [Step 5: BLASTp]       |  [Step 7: HMMER]  |
    v                         v                    |
blastp_{SP}_{DB}.tsv    pfam_domains.txt           |
    |                         |                    |
    |  [Step 6: Combine       |                    |
    |   BLAST + DESeq2]       |                    |
    v                         |                    |
{SP}_annotated.tsv            |                    |
    |                         |                    |
    +-------+-----------------+                    |
            |                                      |
            |  [Step 8: Compare species]           |
            v                                      |
    DC_vs_DG_comparison.tsv                        |
            |                                      |
            +--------------------------------------+
            |
            |  [Step 9: Gene families]
            v
    CYP_verified.tsv, OMT_verified.tsv
            |
            |  [Step 10: Publication plots]
            v
    ma_plot.pdf, volcano_plot.pdf, pca_plot.pdf,
    sample_correlation_heatmap.pdf,
    cyp_heatmap.pdf, omt_heatmap.pdf,
    cyp_heatmap_combined.pdf, omt_heatmap_combined.pdf
```

---

## Scripts Reference

| Script | Purpose |
|--------|---------|
| `parse_refseq_gtf.py` | Smart GTF parser with auto XM lookup |
| `extract_cds_from_gtf.py` | Extract CDS with enhanced FASTA headers |
| `prepare_blast_input.sh` | Master script: GTF to FASTA (bash) |
| `extract_cds__for_blast.sbatch` | Master script: GTF to FASTA (SLURM) |
| `run_translate_cds.sbatch` | CDS nucleotide to protein translation |
| `blastp_discoveryfilter.sbatch` | BLASTp with resume support |
| `run_combine_filter.sbatch` | Merge BLAST + DESeq2, filter candidates |
| `download_protein_databases.sbatch` | Download Pfam + PROSITE databases |
| `run_hmmer.sbatch` | HMMER domain scan against Pfam |
| `compare_species.py` | DC vs DG differential expression comparison |
| `run_compare_species.sbatch` | SLURM wrapper for species comparison |
| `extract_gene_families.py` | Single-species CYP + OMT extraction |
| `extract_gene_families_combined.py` | Two-species combined gene family extraction |
| `generate_family_heatmap.py` | Expression heatmaps for gene families |
| `run_gene_family_heatmap.sbatch` | SLURM: single-species families + heatmaps |
| `run_gene_family_combined.sbatch` | SLURM: combined two-species families + heatmaps |
| `pydeseq2_generate_plots.py` | Publication plots: MA, volcano, PCA, correlation, heatmaps |
| `run_pydeseq2_step3_plots.sbatch` | SLURM wrapper for all publication plots |
