# Directory Structure

Complete file map of the `rnaseq` repository after reorganization.
Scripts are grouped by pipeline stage in numbered directories so the
execution order is obvious.

---

## Top Level

```
rnaseq/
├── README.md                 Landing page — quick start, layout overview, key commands
├── environment.yml           Conda environment definition (all tools + Python deps)
├── pyproject.toml            Python package config — makes rna_pipeline installable
│
├── scripts/                  User-facing workflow scripts (sbatch + Python helpers)
├── rna_pipeline/             Core Python library (QC, alignment, assembly engine)
├── docs/                     All documentation
└── archive/                  Deprecated scripts and old docs (kept for reference)
```

---

## `scripts/` — Workflow Scripts by Pipeline Stage

Each subdirectory contains `.sbatch` files (SLURM jobs you submit) and `.py`
helper scripts (called automatically by the batch jobs). Run directories in
numerical order.

```
scripts/
├── 01_qc/                           Quality control
│   ├── run_fastqc.sbatch            FastQC + MultiQC on all raw samples
│   └── run_fastp.sbatch             Trim/filter reads with fastp
│
├── 02_alignment/                    Read alignment
│   ├── run_genome_index.sbatch      Build STAR genome index (one-time)
│   ├── run_star_align.sbatch        STAR align — single sample group
│   ├── run_star_align_all.sbatch    STAR align — all samples
│   └── run_samtools_sort.sbatch     Sort BAM files with samtools
│
├── 03_assembly/                     De novo assembly (Trinity)
│   ├── run_trinity.sbatch           Trinity — single sample
│   ├── run_trinity_all.sbatch       Trinity — all samples sequentially
│   ├── run_trinity_array.sbatch     Trinity — SLURM array parallelism
│   └── run_trinity_rsem.sbatch      Trinity + RSEM quantification
│
├── 04_counting/                     Read counting → count matrix
│   ├── run_featurecounts.sbatch     featureCounts from BAM files
│   ├── run_count_matrix.sbatch      Merge per-sample counts into matrix
│   └── build_count_matrix.py        Python: STAR ReadsPerGene → count matrix
│
├── 05_pydeseq2/                     Differential expression (3-step) + all-in-one
│   ├── run_step1_analysis.sbatch    Step 1: unfiltered DE analysis
│   ├── run_step2_filter.sbatch      Step 2: filter by padj / log2FC
│   ├── run_step3_plots.sbatch       Step 3: MA, volcano, PCA, heatmaps
│   ├── pydeseq2_run_analysis.py     Python: PyDESeq2 model fitting
│   ├── pydeseq2_filter_results.py   Python: statistical filtering
│   ├── pydeseq2_generate_plots.py   Python: publication-quality plots
│   └── R_pydeseq2/                  All-in-one DESeq2 + visualization script
│       ├── run_R_pydeq2.sbatch
│       ├── R_pydeq2.py
│       ├── README.md
│       └── requirements.txt
│
├── 06_blast/                        BLAST annotation
│   ├── run_blast_input.sbatch       Parse GTF + extract CDS for BLAST
│   ├── run_extract_cds.sbatch       Extract CDS sequences from GTF + genome
│   ├── run_translate_cds.sbatch     Translate CDS nucleotide → protein
│   ├── run_blastp_discovery.sbatch  BLASTp — discovery mode (permissive)
│   ├── run_blastp_strict.sbatch     BLASTp — strict mode (high-confidence)
│   ├── run_blastx_discovery.sbatch  BLASTx — discovery mode
│   ├── run_blastx_strict.sbatch     BLASTx — strict mode
│   ├── run_combine_blast_deseq.sbatch  Merge BLAST hits + DESeq2 stats
│   ├── extract_cds_from_gtf.py      Python: GTF + genome → CDS FASTA
│   ├── translate_cds_to_protein.py  Python: nucleotide FASTA → protein FASTA
│   ├── parse_refseq_gtf.py          Python: parse RefSeq GTF → gene annotation TSV
│   ├── join_gtf_with_gene_ids.py    Python: inner-join parsed GTF with gene ID list
│   ├── combine_blast_deseq.py       Python: join BLAST hits with DE results
│   └── filter_combined_results.py   Python: apply cutoffs to combined file
│
├── 07_domains/                      Protein domain & motif analysis
│   ├── run_hmmer.sbatch             HMMER Pfam domain scan
│   └── run_prosite.sbatch           PROSITE motif scan (EMBOSS)
│
├── 08_gene_families/                Gene family extraction & CYP450 analysis
│   ├── run_cyp450_database.sbatch   Build CYP master list (HMMER + keyword)
│   ├── run_cyp_extract.sbatch       Intersect CYP list with DESeq2 + extract proteins
│   ├── run_filter_genelist.sbatch   Gene list → PyDESeq2 → filter → candidates
│   ├── run_gene_families.sbatch     Extract families from annotated results (single species)
│   ├── run_gene_families_combined.sbatch  Extract families across DC + DG
│   ├── cyp_hmmer_scan.py            Python: HMMER scan for CYP Pfam domains
│   ├── cyp_gff_keyword_search.py    Python: GTF keyword search for "cytochrome P450"
│   ├── cyp_combine_results.py       Python: combine CYP-specific BLAST + expression
│   ├── cyp_intersect_pydeseq2.py    Python: CYP gene list ∩ unfiltered DESeq2
│   ├── cyp_extract_proteins.py      Python: extract protein FASTA for CYP subset
│   ├── filter_count_by_genelist.py  Python: run DE on full matrix, subset to gene list
│   ├── extract_gene_families.py     Python: extract CYP/OMT from annotated results
│   ├── extract_gene_families_combined.py  Python: families across DC + DG combined
│   └── generate_family_heatmap.py   Python: heatmap for a specific gene family
│
├── 09_comparison/                   Cross-species analysis
│   ├── run_compare_species.sbatch   Compare DC vs DG DE results
│   └── compare_species.py           Python: cross-species DE comparison
│
├── 10_post_analysis/                Phylogenetics & genomic clustering
│   ├── run_post_analysis.sbatch     Wrapper: phylo + clustering in one job
│   ├── run_phylo_heatmap.sbatch     Phylogenetic tree + expression heatmap
│   ├── run_genomic_clustering.sbatch  Genomic cluster analysis
│   ├── phylo_heatmap.py             Python: MAFFT → FastTree → tree + heatmap
│   └── genomic_cluster_analysis.py  Python: GTF coordinates → clusters + distance plots
│
├── 11_verify/                       Verification
│   ├── run_verify_genelist.sbatch   6-check verification vs previous student
│   └── verify_genelist.py           Python: gene overlap, direction, correlation checks
│
├── setup/                           One-time database downloads
│   ├── download_blastdb.sbatch      Download SwissProt / NR databases
│   ├── download_protein_databases.sbatch  Download protein reference databases
│   └── extract_daucus_carota_db.sbatch    Extract carrot-specific protein DB
```

---

## `rna_pipeline/` — Core Python Library

Installable package (`pip install -e .`). Handles upstream pipeline stages
(QC → alignment → assembly) via the CLI or programmatic API.

```
rna_pipeline/
├── __init__.py              Package marker
├── cli.py                   Argparse CLI — parses --mode (qc, index, align, etc.)
├── main.py                  Pipeline orchestration — dispatches to tool workflows
├── logging_setup.py         Console + file logging configuration
│
├── tools/                   Tool-specific command builders
│   ├── __init__.py
│   ├── qc.py               FastQC + MultiQC command builders
│   ├── star.py              STAR index + align command builders
│   ├── trinity.py           Trinity assembly command builder
│   ├── featurecounts.py     featureCounts (Subread) read counting
│   ├── blast.py             BLASTp + BLASTx search command builders
│   ├── hmmer.py             HMMER hmmscan (Pfam domain scanning)
│   ├── prosite.py           EMBOSS patmatmotifs (motif scanning)
│   └── pydeseq2.py          PyDESeq2 differential expression helpers
│
├── runners/
│   └── local.py             Subprocess runner (executes built commands)
│
└── utils/
    ├── io_utils.py           File I/O helpers
    └── sys_utils.py          System utilities
```

---

## `docs/` — Documentation

```
docs/
├── PIPELINE_WORKFLOW.md      Master runbook — Steps 0–12, full commands + expected output
├── HMMER_PROSITE_GUIDE.md    Domain analysis setup & usage
├── DIRECTORY_STRUCTURE.md    This file
├── CONCEPTS.md               RNA-seq background for beginners
└── SETUP.md                  Conda environment, HPC layout, first-time config
```

---

## `archive/` — Deprecated Files

Old scripts and documentation kept for reference. Not part of the active
pipeline.

```
archive/
├── ReadME/                   Old documentation files (superseded by docs/)
├── FINAL_WORKFLOW.md         Old workflow notes
├── build_count_matrix_old.py Legacy count matrix builder
├── pydeseq2_countmatrix.py   Old single-script DESeq2 analysis
├── pydeseq2_basic.py         Early prototype
├── old_pydeq2_analysis.py    Previous analysis script
├── compare_count_matrices.py Old matrix comparison utility
├── fetch_proteins_edirect.py Old NCBI protein fetch script
├── blast_sequences.py        Empty placeholder (never used)
└── rna_pipeline_tools_*.py   Unused tool modules from rna_pipeline
```

---

## Execution Order

Run each stage for **both DC and DG** before moving to cross-species steps.

```
Stage 01  QC              →  Are reads high quality?
Stage 02  Alignment       →  Map reads to reference genome (STAR)
Stage 03  Assembly        →  (Optional) De novo assembly if needed
Stage 04  Counting        →  Build gene × sample count matrix
Stage 05  DESeq2          →  Differential expression (3 steps)
Stage 06  BLAST           →  Annotate genes with function
Stage 07  Domains         →  HMMER + PROSITE domain verification
Stage 08  Gene Families   →  Extract CYP/OMT families
Stage 09  Comparison      →  DC vs DG cross-species analysis
Stage 10  Post-Analysis   →  Phylogenetic trees, genomic clustering
Stage 11  Verify          →  Validate against previous results
```
