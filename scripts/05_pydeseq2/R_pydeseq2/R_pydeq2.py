#!/usr/bin/env python3
"""

Deseq2 Reference Script
Steps:
  0.  Load count matrix (from featureCounts output)
  1.  Filter low-count genes (rowSums > 20)
  2.  Define experimental design (condition metadata)
  3.  Run PyDESeq2 (equivalent of DESeq2 in R)
  4.  Save normalized counts
  5.  QC plots: PCA + dispersion-like plot
  6.  Extract DE results (contrast: A vs B)
  7.  Filter significant genes (padj < 0.05, matching R subset() behavior)
  8.  Heatmap of all significant genes
  9.  Volcano plot + MA plot
  10. Filter P450 genes → P450 heatmap
  11. Filter methyltransferases → MT heatmap
  12. Phylogenetic tree + heatmap (MAFFT → FastTree → combined plot)
  13. Fetch gene annotations from NCBI (Python Entrez = R's rentrez)
  14. Genomic proximity / cluster analysis

Required Python packages:
  numpy, pandas, matplotlib, seaborn, scipy, pydeseq2, biopython

  Install all at once:
    pip install numpy pandas matplotlib seaborn scipy pydeseq2 biopython

  Or via conda:
    conda install -c conda-forge numpy pandas matplotlib seaborn scipy
    pip install pydeseq2 biopython

  External tools (for phylogenetic tree building — optional):
    mafft, FastTree  (loaded via 'module load' on HPC or installed via conda)

Usage:
  python scripts/R_pydeseq2/R_pydeq2.py \\
      --counts 03_count_tables/00_1_DC/gene_count_matrix.tsv \\
      --metadata 03_count_tables/00_1_DC/sample_metadata.tsv \\
      --p450-list 07_NRdatabase/sukman_database/P450_list_RefSeq.txt \\
      --gtf 04_reference/dc_genomic.gtf \\
      --protein-fasta 04_reference/GCF_001625215.2_DH1_v3.0_protein.faa \\
      --output-dir 06_analysis/R_pydeseq2_DC/

  Minimal (uses defaults):
  python scripts/R_pydeseq2/R_pydeq2.py

  Custom contrast (e.g. Root vs Leaf, Leaf = baseline):
  python scripts/R_pydeseq2/R_pydeq2.py --contrast-A E --contrast-B N
"""

import argparse
import subprocess
import sys
import warnings
from pathlib import Path

# ═══════════════════════════════════════════════════════════════════════════════
# DEPENDENCY CHECK — verifies all required packages before starting
# ═══════════════════════════════════════════════════════════════════════════════

REQUIRED_PACKAGES = {
    "numpy":      {"import": "numpy",        "pip": "numpy"},
    "pandas":     {"import": "pandas",       "pip": "pandas"},
    "matplotlib": {"import": "matplotlib",   "pip": "matplotlib"},
    "seaborn":    {"import": "seaborn",      "pip": "seaborn"},
    "scipy":      {"import": "scipy",        "pip": "scipy"},
    "pydeseq2":   {"import": "pydeseq2",     "pip": "pydeseq2"},
}
OPTIONAL_PACKAGES = {
    "biopython":  {"import": "Bio",          "pip": "biopython",
                   "used_for": "NCBI Entrez queries (step 13) + phylo tree reading (step 12)"},
}

def check_dependencies():
    """Verify all required and optional packages are importable."""
    missing_required = []
    missing_optional = []

    for name, info in REQUIRED_PACKAGES.items():
        try:
            __import__(info["import"])
        except ImportError:
            missing_required.append(info["pip"])

    for name, info in OPTIONAL_PACKAGES.items():
        try:
            __import__(info["import"])
        except ImportError:
            missing_optional.append((info["pip"], info["used_for"]))

    if missing_required:
        print("ERROR: Missing required Python packages:")
        for pkg in missing_required:
            print(f"  - {pkg}")
        print()
        print("Install all at once:")
        print(f"  pip install {' '.join(missing_required)}")
        print()
        print("Or on the HPC, activate your conda environment first:")
        print("  conda activate rnaseq")
        sys.exit(1)

    if missing_optional:
        print("WARNING: Missing optional packages (some steps will be skipped):")
        for pkg, reason in missing_optional:
            print(f"  - {pkg}  →  needed for: {reason}")
        print(f"  Install: pip install {' '.join(p for p, _ in missing_optional)}")
        print()

    # Check external tools (optional)
    for tool, used_for in [("mafft", "multiple sequence alignment"),
                           ("FastTree", "phylogenetic tree building")]:
        try:
            subprocess.run([tool, "--help"], capture_output=True, timeout=5)
        except (FileNotFoundError, subprocess.TimeoutExpired):
            pass  # reported later when needed

check_dependencies()

import numpy as np
import pandas as pd

warnings.filterwarnings('ignore')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist

# HPC default paths (override with CLI args)
BASE_DIR = "/projects/tholl_lab_1/daisy_analysis"
DEFAULTS = {
    "counts": f"{BASE_DIR}/03_count_tables/00_1_DC/gene_count_matrix.tsv",
    "metadata": f"{BASE_DIR}/03_count_tables/00_1_DC/sample_metadata.tsv",
    "p450_list": f"{BASE_DIR}/07_NRdatabase/sukman_database/P450_list_RefSeq.txt",
    "mt_list": f"{BASE_DIR}/07_NRdatabase/sukman_database/MT_expression_refseq.txt",
    "gtf": f"{BASE_DIR}/04_reference/dc_genomic.gtf",
    "protein_fasta": f"{BASE_DIR}/04_reference/GCF_001625215.2_DH1_v3.0_protein.faa",
    "output_dir": f"{BASE_DIR}/06_analysis/R_pydeseq2_DC",
}


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 0-1: Load and filter count matrix
# ═══════════════════════════════════════════════════════════════════════════════
# R equivalent:
#   Counts <- read.delim(file.choose(), header=TRUE, row.names=1)
#   Counts <- Counts_merged[, -c(1:5)]
#   Counts <- Counts[which(rowSums(Counts) > 20),]

def load_count_matrix(counts_path, metadata_path, min_counts=20):
    """Load count matrix and metadata, filter low-count genes."""
    print("=" * 70)
    print("STEP 0-1: Load count matrix + filter low-count genes")
    print("=" * 70)

    counts = pd.read_csv(counts_path, sep='\t', index_col=0)

    # Drop annotation columns if present (featureCounts raw output has Chr, Start, End, Strand, Length)
    annotation_cols = {'Chr', 'Start', 'End', 'Strand', 'Length'}
    drop_cols = [c for c in counts.columns if c in annotation_cols]
    if drop_cols:
        print(f"  Dropping annotation columns: {drop_cols}")
        counts = counts.drop(columns=drop_cols)

    # Ensure integer counts
    counts = counts.apply(pd.to_numeric, errors='coerce').fillna(0).astype(int)

    print(f"  Raw count matrix: {counts.shape[0]} genes × {counts.shape[1]} samples")

    # Load metadata
    metadata = pd.read_csv(metadata_path, sep='\t')
    if 'sample' in metadata.columns:
        metadata = metadata.set_index('sample')

    common = sorted(set(counts.columns) & set(metadata.index))
    counts = counts[common]
    metadata = metadata.loc[common]
    print(f"  Samples matched: {len(common)}")

    # Filter: rowSums > min_counts (student used > 20)
    gene_sums = counts.sum(axis=1)
    keep = gene_sums > min_counts
    counts_filtered = counts.loc[keep]
    print(f"  Genes before filter: {counts.shape[0]}")
    print(f"  Genes after filter (rowSums > {min_counts}): {counts_filtered.shape[0]}")
    print(f"  Genes removed: {counts.shape[0] - counts_filtered.shape[0]}")

    return counts_filtered, metadata


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 2-4: Run PyDESeq2 (equivalent of DESeq2 in R)
# ═══════════════════════════════════════════════════════════════════════════════
# R equivalent:
#   dds <- DESeqDataSetFromMatrix(countData=Counts, colData=coldata, design=~condition)
#   dds <- DESeq(dds, minReplicatesForReplace=Inf)
#   normalized_counts <- counts(dds, normalized=TRUE)
#   res <- results(dds, cooksCutoff=FALSE, independentFiltering=FALSE,
#                  contrast=c("condition","R","L"))

def run_pydeseq2(counts, metadata, contrast_factor, contrast_A, contrast_B,
                 output_dir):
    """Run PyDESeq2 differential expression analysis."""
    print("\n" + "=" * 70)
    print("STEP 2-4: Run PyDESeq2 (R equivalent: DESeq2)")
    print("=" * 70)

    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError:
        print("ERROR: PyDESeq2 not installed. Install: pip install pydeseq2")
        sys.exit(1)

    counts_t = counts.T
    meta = metadata.copy()
    if contrast_factor not in meta.columns:
        raise ValueError(f"'{contrast_factor}' not in metadata columns: {list(meta.columns)}")

    print(f"  Design: ~ {contrast_factor}")
    print(f"  Contrast: {contrast_A} vs {contrast_B}")
    print(f"  Samples: {counts_t.shape[0]}, Genes: {counts_t.shape[1]}")

    # Create DeseqDataSet — refit_cooks=False to match R's minReplicatesForReplace=Inf
    dds = DeseqDataSet(
        counts=counts_t,
        metadata=meta,
        design_factors=[contrast_factor],
        refit_cooks=False,
        n_cpus=4
    )

    # Run full DESeq2 pipeline in one call
    print("  Running DESeq2 pipeline...")
    dds.deseq2()

    # Normalized counts = raw count/ size factor
    # Size factors are per-sample scaling factors computed by dds.deseq2()
    if "size_factors" in dds.obs.columns:
        sf_values = dds.obs["size_factors"].values
        norm_counts = (counts_t.div(sf_values, axis=0)).T
        print(f"  Size factors: {dict(zip(dds.obs_names, sf_values.round(4)))}")
    else:
        print("  WARNING: size factors not found in dds.obs after deseq2().")
        print(f"  Available dds.obs columns: {list(dds.obs.columns)}")
        print("  Saving raw counts instead — check your PyDESeq2 version.")
        norm_counts = counts.astype(float)

    norm_file = output_dir / "normalized_counts_refseq.csv"
    norm_counts.to_csv(norm_file)
    print(f"  Saved: {norm_file}")

    # Get results with contrast
    contrast = (contrast_factor, contrast_A, contrast_B)
    stat_res = DeseqStats(dds, contrast=contrast, n_cpus=4)
    stat_res.summary()
    results_df = stat_res.results_df

    results_file = output_dir / "res_all.csv"
    results_df.to_csv(results_file)
    print(f"  Saved: {results_file}")
    print(f"  Total genes tested: {len(results_df)}")

    valid = results_df['padj'].notna()
    sig = results_df.loc[valid & (results_df['padj'] <= 0.05)]
    print(f"  Genes with padj ≤ 0.05: {len(sig)}")

    return dds, results_df, norm_counts


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 5: QC Plots — PCA + dispersion
# ═══════════════════════════════════════════════════════════════════════════════
# R equivalent:
#   vsdata <- vst(dds, blind=FALSE)
#   plotPCA(vsdata, intgroup="condition")
#   plotDispEsts(dds)

def _confidence_ellipse(x, y, ax, n_std=1.96, **kwargs):
    """Draw a 95% confidence ellipse around a set of points."""
    from matplotlib.patches import Ellipse
    import matplotlib.transforms as transforms
    if len(x) < 3:
        return
    cov = np.cov(x, y)
    pearson = cov[0, 1] / (np.sqrt(cov[0, 0] * cov[1, 1]) + 1e-12)
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                       **kwargs)
    scale_x = np.sqrt(cov[0, 0]) * n_std
    scale_y = np.sqrt(cov[1, 1]) * n_std
    transf = (transforms.Affine2D()
              .rotate_deg(45)
              .scale(scale_x, scale_y)
              .translate(np.mean(x), np.mean(y)))
    ellipse.set_transform(transf + ax.transData)
    ax.add_patch(ellipse)


def generate_pca_plot(norm_counts, metadata, contrast_factor, output_dir):
    """PCA of top-500 most variable genes with 95% confidence ellipses."""
    print("\n" + "=" * 70)
    print("STEP 5a: PCA plot (R equivalent: vst + plotPCA)")
    print("=" * 70)

    log_data = np.log2(norm_counts + 1)

    n_top = 500
    gene_var = log_data.var(axis=1)
    n_top_actual = min(n_top, len(gene_var))
    top_genes = gene_var.nlargest(n_top_actual).index
    subset = log_data.loc[top_genes]

    X = subset.values.T
    X_centered = X - X.mean(axis=0)

    U, S, Vt = np.linalg.svd(X_centered, full_matrices=False)
    explained = (S ** 2) / (S ** 2).sum() * 100
    pc1 = U[:, 0] * S[0]
    pc2 = U[:, 1] * S[1]

    cond_palette = {'R': '#d35400', 'L': '#27ae60',
                    'root': '#d35400', 'leaf': '#27ae60'}
    cond_labels = {'R': 'Root', 'L': 'Leaf',
                   'root': 'Root', 'leaf': 'Leaf'}

    fig, ax = plt.subplots(figsize=(8, 6))
    fig.patch.set_alpha(0)
    ax.set_facecolor('none')

    group_points = {}
    for i, sample in enumerate(subset.columns):
        cond = metadata.loc[sample, contrast_factor] if sample in metadata.index else 'unknown'
        color = cond_palette.get(cond, '#bdc3c7')
        ax.scatter(pc1[i], pc2[i], c=color, s=100, edgecolors='white',
                   linewidths=0.8, zorder=3)
        ax.annotate(sample, (pc1[i], pc2[i]), fontsize=8, fontweight='bold',
                    xytext=(6, 6), textcoords='offset points')
        group_points.setdefault(cond, ([], []))
        group_points[cond][0].append(pc1[i])
        group_points[cond][1].append(pc2[i])

    for cond, (xs, ys) in group_points.items():
        color = cond_palette.get(cond, '#bdc3c7')
        _confidence_ellipse(np.array(xs), np.array(ys), ax, n_std=1.96,
                            facecolor=color, alpha=0.15, edgecolor=color,
                            linewidth=1.5, linestyle='--')

    handles = []
    seen_labels = set()
    for cond_code, color in cond_palette.items():
        label = cond_labels.get(cond_code, cond_code)
        if label not in seen_labels:
            seen_labels.add(label)
            handles.append(plt.scatter([], [], c=color, s=60,
                                       edgecolors='white', label=label))

    ax.set_xlabel(f'PC1: {explained[0]:.1f}% variance', fontsize=12)
    ax.set_ylabel(f'PC2: {explained[1]:.1f}% variance', fontsize=12)
    ax.set_title(f'PCA  |  Top {n_top_actual} most variable genes',
                 fontsize=13, fontweight='bold')
    ax.legend(handles=handles, title='Tissue', fontsize=9,
              title_fontsize=10, loc='best', framealpha=0.9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.axhline(0, color='gray', linewidth=0.3, alpha=0.5)
    ax.axvline(0, color='gray', linewidth=0.3, alpha=0.5)

    plt.tight_layout()
    out = output_dir / "pca_plot.pdf"
    fig.savefig(out, dpi=300, bbox_inches='tight', transparent=True)
    fig.savefig(output_dir / "pca_plot.png", dpi=300, bbox_inches='tight',
                transparent=True)
    fig.savefig(output_dir / "pca_plot.svg", bbox_inches='tight',
                transparent=True)
    plt.close(fig)
    print(f"  Saved: {out}")


def generate_dispersion_plot(dds, output_dir):
    """Dispersion estimates plot (R equivalent: plotDispEsts)."""
    print("STEP 5b: Dispersion plot (R equivalent: plotDispEsts)")

    try:
        gene_means = np.array(dds.varm["_normed_means"]).flatten() if "_normed_means" in dds.varm else None
        dispersions = np.array(dds.varm["dispersions"]).flatten() if "dispersions" in dds.varm else None

        if gene_means is None or dispersions is None:
            print("  Skipping dispersion plot (dispersions not available in dds object)")
            return

        fig, ax = plt.subplots(figsize=(8, 6))
        valid = (gene_means > 0) & (dispersions > 0) & np.isfinite(gene_means) & np.isfinite(dispersions)
        ax.scatter(np.log10(gene_means[valid]), np.log10(dispersions[valid]),
                   alpha=0.3, s=5, c='black')
        ax.set_xlabel('log10(mean normalized count)')
        ax.set_ylabel('log10(dispersion)')
        ax.set_title('Dispersion Estimates')
        ax.grid(alpha=0.3)

        out = output_dir / "dispersion_plot.pdf"
        fig.savefig(out, dpi=300, bbox_inches='tight')
        fig.savefig(output_dir / "dispersion_plot.png", dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"  Saved: {out}")
    except Exception as e:
        print(f"  Skipping dispersion plot: {e}")


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 7: Filter significant genes
# ═══════════════════════════════════════════════════════════════════════════════
# R equivalent (replicating subset(res, padj < 0.05, log2FoldChange >= 2)):
#   The third arg to subset() is `select`, not a filter — so only padj < 0.05
#   is applied as the row filter.  The lfc_cutoff is used only for up/down
#   categorization and for volcano/MA plot aesthetics.

def filter_significant(results_df, padj_cutoff=0.05, lfc_cutoff=2.0, output_dir=None):
    """Filter significant DE genes (padj-only, matching R subset() behavior)."""
    print("\n" + "=" * 70)
    print("STEP 7: Filter significant genes")
    print("=" * 70)

    valid = results_df['padj'].notna() & results_df['log2FoldChange'].notna()
    sig = results_df.loc[valid & (results_df['padj'] <= padj_cutoff)]

    up = sig[sig['log2FoldChange'] >= lfc_cutoff]
    down = sig[sig['log2FoldChange'] <= -lfc_cutoff]

    print(f"  padj ≤ {padj_cutoff} (matching R subset() behavior)")
    print(f"  Significant genes: {len(sig)}")
    print(f"    Upregulated   (log2FC >= {lfc_cutoff}): {len(up)}")
    print(f"    Downregulated (log2FC <= -{lfc_cutoff}): {len(down)}")

    if output_dir:
        sig.to_csv(output_dir / "res_all_significant.csv")
        up.to_csv(output_dir / "res_Upregulated.csv")
        down.to_csv(output_dir / "res_Downregulated.csv")
        print(f"  Saved: res_all_significant.csv, res_Upregulated.csv, res_Downregulated.csv")

    return sig, up, down


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 8: Heatmap of all significant genes
# ═══════════════════════════════════════════════════════════════════════════════
# R equivalent:
#   allsig <- merge(normalized_counts, significant_res_RvsL, by=0)
#   sigCounts <- allsig[,2:7]
#   pheatmap(log2(sigCounts+1), scale='row', show_rownames=FALSE)

def generate_heatmap(norm_counts, gene_set, output_dir, filename="heatmap_all_significant",
                     title="All Significant DE Genes", show_labels=False, fontsize=6,
                     scale_rows=True):
    """
    Heatmap of log2(normalized_counts + 1).

    Parameters
    ----------
    scale_rows : bool
        True  → z-score each row (R's pheatmap scale='row'). Colors = z-score.
        False → display raw log2(count+1) values.  Colors = log2 expression.
    """
    overlap = sorted(set(gene_set) & set(norm_counts.index))
    if not overlap:
        print(f"  No genes found for heatmap '{filename}' — skipping")
        return

    data = norm_counts.loc[overlap]
    log2_data = np.log2(data + 1)

    if scale_rows:
        # Row z-score scaling (R's pheatmap scale='row')
        row_means = log2_data.mean(axis=1)
        row_stds  = log2_data.std(axis=1).replace(0, 1)
        plot_data = log2_data.sub(row_means, axis=0).div(row_stds, axis=0)
        cmap = LinearSegmentedColormap.from_list('bwr', ['#2166AC', 'white', '#B2182B'])
        finite_vals = plot_data.values[np.isfinite(plot_data.values)]
        abs_max = max(abs(finite_vals.min()), abs(finite_vals.max())) if len(finite_vals) else 1
        vmin, vmax = -abs_max, abs_max
        cbar_label = 'z-score'
    else:
        # Raw log2(count+1) values — sequential colormap (0 = no expression)
        plot_data = log2_data
        cmap = 'viridis'
        vmin, vmax = 0, plot_data.values.max()
        cbar_label = 'log2(count + 1)'

    n = len(overlap)
    n_samples = len(data.columns)
    gene_label_size = max(5, min(8, 240 // max(n, 1)))
    sample_label_size = max(10, min(14, 160 // max(n_samples, 1)))
    height = max(14, min(60, 0.35 * n + 6))
    width  = max(height * 0.8, 1.5 * n_samples + 8)

    try:
        g = sns.clustermap(
            plot_data,
            cmap=cmap,
            center=0 if scale_rows else None,
            figsize=(width, height),
            row_cluster=True,
            col_cluster=True,
            linewidths=0.2,
            linecolor='white',
            vmin=vmin, vmax=vmax,
            cbar_kws={'label': cbar_label, 'shrink': 0.4},
            yticklabels=show_labels,
            xticklabels=True,
            dendrogram_ratio=(0.15, 0.06),
            method='ward',
            colors_ratio=0.02,
            cbar_pos=(0.02, 0.92, 0.03, 0.06),
        )

        g.cax.tick_params(labelsize=8)
        g.cax.set_ylabel(cbar_label, fontsize=9)

        g.ax_heatmap.tick_params(axis='y', labelsize=gene_label_size,
                                 length=0, pad=4)
        g.ax_heatmap.tick_params(axis='x', labelsize=sample_label_size,
                                 length=0, pad=4)
        for label in g.ax_heatmap.get_xticklabels():
            label.set_rotation(45)
            label.set_ha('right')
        g.ax_heatmap.set_ylabel('')
        g.ax_heatmap.set_xlabel('')

        for ax_dend in [g.ax_row_dendrogram, g.ax_col_dendrogram]:
            for coll in ax_dend.collections:
                coll.set_linewidth(1.5)

        g.fig.suptitle(title, fontsize=14, fontweight='bold', y=1.02)

        out = output_dir / f"{filename}.png"
        g.savefig(out, dpi=300, bbox_inches='tight')
        g.savefig(output_dir / f"{filename}.pdf", dpi=300, bbox_inches='tight')

        g.fig.patch.set_alpha(0)
        g.ax_heatmap.set_facecolor('none')
        g.savefig(output_dir / f"{filename}_transparent.png", dpi=300,
                  bbox_inches='tight', transparent=True)
        plt.close('all')
        print(f"  Saved: {out} ({n} genes)")
        print(f"  Saved: {output_dir / f'{filename}_transparent.png'}")
    except Exception as e:
        print(f"  Heatmap error for '{filename}': {e}")


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 9: Volcano plot + MA plot
# ═══════════════════════════════════════════════════════════════════════════════
# R equivalent:
#   EnhancedVolcano(res.df, lab=rownames(res), x='log2FoldChange', y='pvalue')
#   plotMA(res)

def generate_volcano_plot(results_df, output_dir, padj_cutoff=0.05, lfc_cutoff=2.0,
                          contrast_A="R", contrast_B="L"):
    """Enhanced volcano plot (R equivalent: EnhancedVolcano)."""
    print("\n" + "=" * 70)
    print("STEP 9a: Volcano plot (R equivalent: EnhancedVolcano)")
    print("=" * 70)

    df = results_df.copy()
    df['neg_log10_pval'] = -np.log10(df['pvalue'].clip(lower=1e-300))

    valid = df['pvalue'].notna() & df['log2FoldChange'].notna()
    df = df.loc[valid]

    # 3-color scheme matching R's keyvals: color based purely on log2FC (±2)
    colors = np.full(len(df), 'black', dtype=object)
    colors[(df['log2FoldChange'] > 2).values] = 'red'
    colors[(df['log2FoldChange'] < -2).values] = 'blue'

    fig, ax = plt.subplots(figsize=(10, 8))

    for color, label in [('black', '-2<log2FC<2'),
                          ('red', 'log2FC>2'),
                          ('blue', 'log2FC<-2')]:
        mask = colors == color
        if mask.any():
            ax.scatter(df.loc[mask, 'log2FoldChange'], df.loc[mask, 'neg_log10_pval'],
                       c=color, s=8, alpha=0.6, label=label, edgecolors='none')

    # Horizontal line at -log10(1e-13) = 13, matching R's pCutoff = 10e-14
    ax.axhline(13, color='gray', linestyle='--', linewidth=0.8)

    ax.set_xlim(-20, 20)
    ax.set_xticks(range(-20, 22, 2))
    ax.set_xlabel(r'$\log_2$ fold change', fontsize=12)
    ax.set_ylabel(r'$-\log_{10}$ p-value', fontsize=12)
    ax.set_title(f'Volcano Plot — {contrast_A} vs {contrast_B}')
    ax.legend(fontsize=9, markerscale=2)
    ax.grid(alpha=0.2)

    n_up = (colors == 'red').sum()
    n_down = (colors == 'blue').sum()
    ax.text(0.02, 0.98, f'Up: {n_up}\nDown: {n_down}',
            transform=ax.transAxes, fontsize=10, va='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    out = output_dir / "volcano_plot.png"
    fig.savefig(out, dpi=300, bbox_inches='tight')
    fig.savefig(output_dir / "volcano_plot.pdf", dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out}")


def generate_ma_plot(results_df, output_dir, padj_cutoff=0.05, lfc_cutoff=2.0):
    """MA plot (R equivalent: plotMA)."""
    print("STEP 9b: MA plot (R equivalent: plotMA)")

    df = results_df.copy()
    valid = df['baseMean'].notna() & df['log2FoldChange'].notna()
    df = df.loc[valid]

    sig_up = (df['padj'] <= padj_cutoff) & (df['log2FoldChange'] >= lfc_cutoff)
    sig_down = (df['padj'] <= padj_cutoff) & (df['log2FoldChange'] <= -lfc_cutoff)
    ns = ~sig_up & ~sig_down

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(np.log10(df.loc[ns, 'baseMean'] + 1), df.loc[ns, 'log2FoldChange'],
               c='gray', s=5, alpha=0.3, label='NS')
    ax.scatter(np.log10(df.loc[sig_up, 'baseMean'] + 1), df.loc[sig_up, 'log2FoldChange'],
               c='red', s=10, alpha=0.6, label=f'Up ({sig_up.sum()})')
    ax.scatter(np.log10(df.loc[sig_down, 'baseMean'] + 1), df.loc[sig_down, 'log2FoldChange'],
               c='blue', s=10, alpha=0.6, label=f'Down ({sig_down.sum()})')

    ax.axhline(0, color='black', linewidth=0.8)
    ax.axhline(lfc_cutoff, color='gray', linestyle='--', linewidth=0.5)
    ax.axhline(-lfc_cutoff, color='gray', linestyle='--', linewidth=0.5)
    ax.set_xlabel(r'$\log_{10}$ mean normalized count')
    ax.set_ylabel(r'$\log_2$ fold change')
    ax.set_title('MA Plot')
    ax.legend(fontsize=9)
    ax.grid(alpha=0.2)

    out = output_dir / "ma_plot.png"
    fig.savefig(out, dpi=300, bbox_inches='tight')
    fig.savefig(output_dir / "ma_plot.pdf", dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out}")


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 10-11: Filter by gene family + heatmaps
# ═══════════════════════════════════════════════════════════════════════════════
# R equivalent:
#   P450_list_refseq <- read.table(file.choose(), header=FALSE)
#   filtered_data_P450 <- res.df[res.df$gene_id %in% P450_list_refseq$V1, ]
#   P450 <- merge(normalized_counts, filtered_data_P450, by=0)
#   pheatmap(log2(P450_counts+1), scale='row')

def load_gene_list(path):
    """Load a gene list from TXT, CSV, or TSV (auto-detects tab-separated .txt)."""
    path = Path(path)
    if not path.exists():
        return []

    with open(path) as f:
        first_line = f.readline().strip()

    if '\t' in first_line or str(path).endswith(('.tsv', '.tab', '.csv')):
        sep = '\t' if '\t' in first_line else ','
        try:
            df = pd.read_csv(path, sep=sep)
            if 'gene_id' in df.columns:
                return df['gene_id'].astype(str).tolist()
            return df.iloc[:, 0].astype(str).tolist()
        except Exception:
            pass

    with open(path) as f:
        return [line.strip() for line in f if line.strip() and not line.startswith('#')]


def filter_gene_family(results_df, gene_list, family_name, output_dir):
    """Filter DE results to a specific gene family (R equivalent: %in% filter)."""
    gene_set = set(gene_list)
    filtered = results_df[results_df.index.isin(gene_set)]

    out = output_dir / f"{family_name}_expression_refseq.csv"
    filtered.to_csv(out)
    print(f"  {family_name}: {len(filtered)} genes found in DE results → {out}")

    return filtered


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 12: Phylogenetic tree + heatmap
# ═══════════════════════════════════════════════════════════════════════════════
# R equivalent:
#   phy_tree <- read.tree(file.choose())
#   p8 <- ggtree(phy_tree) + geom_tiplab()
#   p9 <- gheatmap(p8, log2_counts, offset=0.6, width=10)

def build_tree_if_needed(protein_fasta, expressed_genes, gtf_path, output_dir, family):
    """Extract family proteins, align with MAFFT, build tree with FastTree."""
    fasta_out = output_dir / f"{family}_proteins.fasta"
    alignment_out = output_dir / f"{family}_alignment.fasta"
    tree_out = output_dir / f"{family}_tree.nwk"

    if tree_out.exists():
        print(f"  Tree already exists: {tree_out}")
        return tree_out

    if not Path(protein_fasta).exists():
        print(f"  Protein FASTA not found: {protein_fasta} — skipping tree")
        return None

    # Build gene→protein map from GTF
    gene_to_prot = {}
    gtf = Path(gtf_path)
    if gtf.exists():
        with open(gtf) as f:
            for line in f:
                if line.startswith('#') or '\tCDS\t' not in line:
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                attrs = {}
                for item in fields[8].split(';'):
                    item = item.strip()
                    if item:
                        parts = item.split(' ', 1)
                        if len(parts) == 2:
                            attrs[parts[0]] = parts[1].strip('"')
                gid = attrs.get('gene_id', '')
                pid = attrs.get('protein_id', '')
                if gid in expressed_genes and pid and gid not in gene_to_prot:
                    gene_to_prot[gid] = pid

    if not gene_to_prot:
        print(f"  No protein IDs found for {family} genes — skipping tree")
        return None

    # Extract protein sequences
    wanted_pids = set(gene_to_prot.values())
    pid_to_gene = {v: k for k, v in gene_to_prot.items()}
    written = 0
    writing = False

    with open(protein_fasta) as fin, open(fasta_out, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                acc = line[1:].split()[0]
                acc_base = acc.rsplit('.', 1)[0] if '.' in acc else acc
                if acc in wanted_pids or acc_base in {p.rsplit('.', 1)[0] for p in wanted_pids}:
                    gene_id = pid_to_gene.get(acc, pid_to_gene.get(
                        next((p for p in wanted_pids if p.startswith(acc_base)), ''), acc))
                    fout.write(f">{gene_id}\n")
                    writing = True
                    written += 1
                else:
                    writing = False
            elif writing:
                fout.write(line)

    if written < 3:
        print(f"  Only {written} sequences — need ≥3 for tree. Skipping.")
        return None
    print(f"  Extracted {written} protein sequences → {fasta_out}")

    # MAFFT alignment
    try:
        with open(alignment_out, 'w') as out:
            subprocess.run(['mafft', '--auto', '--thread', '4', str(fasta_out)],
                           stdout=out, stderr=subprocess.PIPE, check=True)
        print(f"  MAFFT alignment → {alignment_out}")
    except (FileNotFoundError, subprocess.CalledProcessError) as e:
        print(f"  MAFFT failed ({e}) — skipping tree")
        return None

    # FastTree
    for cmd_name in ['FastTree', 'fasttree']:
        try:
            with open(tree_out, 'w') as out:
                subprocess.run([cmd_name, str(alignment_out)],
                               stdout=out, stderr=subprocess.PIPE, check=True)
            print(f"  FastTree → {tree_out}")
            return tree_out
        except FileNotFoundError:
            continue
        except subprocess.CalledProcessError as e:
            print(f"  FastTree failed: {e}")
            return None

    print("  FastTree not found — skipping tree")
    return None


def generate_phylo_heatmap(tree_path, norm_counts, gene_set, output_dir,
                           family="P450", title="Phylo + Heatmap"):
    """Combined phylo tree + expression heatmap (R equivalent: ggtree + gheatmap)."""
    print(f"\n  Generating {family} phylo+heatmap...")

    overlap = sorted(set(gene_set) & set(norm_counts.index))
    if len(overlap) < 2:
        print(f"    Only {len(overlap)} genes — need ≥2. Skipping.")
        return

    log2_data = np.log2(norm_counts.loc[overlap] + 1)
    # No row-scaling: matches student's gheatmap(log2_counts) which uses absolute values
    scaled = log2_data

    # Try to load tree for leaf ordering
    leaf_order = None
    if tree_path and Path(tree_path).exists():
        try:
            from Bio import Phylo
            from io import StringIO
            tree = Phylo.read(str(tree_path), 'newick')
            leaf_order = [t.name for t in tree.get_terminals() if t.name and t.name in scaled.index]
        except Exception:
            pass

    if leaf_order and len(leaf_order) >= 2:
        remaining = [g for g in scaled.index if g not in leaf_order]
        ordered = leaf_order + remaining
        scaled = scaled.reindex([g for g in ordered if g in scaled.index])

        # Draw tree + heatmap
        n_genes = len(scaled)
        height = max(6, n_genes * 0.25)
        fig = plt.figure(figsize=(12, height))
        gs = gridspec.GridSpec(1, 3, width_ratios=[3, 8, 0.3], wspace=0.05)
        ax_tree = fig.add_subplot(gs[0])
        ax_heat = fig.add_subplot(gs[1])
        ax_cbar = fig.add_subplot(gs[2])

        # Draw tree
        try:
            from Bio import Phylo
            label_to_y = {name: i for i, name in enumerate(scaled.index)}

            def get_y(clade):
                if clade.is_terminal():
                    return label_to_y.get(clade.name, 0)
                return np.mean([get_y(c) for c in clade.clades])

            def draw_clade(clade, x_start=0):
                x_end = x_start + (clade.branch_length or 0)
                y = get_y(clade)
                ax_tree.plot([x_start, x_end], [y, y], color='black', linewidth=0.8)
                if not clade.is_terminal():
                    child_ys = [get_y(c) for c in clade.clades]
                    ax_tree.plot([x_end, x_end], [min(child_ys), max(child_ys)],
                                color='black', linewidth=0.8)
                    for child in clade.clades:
                        draw_clade(child, x_end)

            draw_clade(tree.root)
            ax_tree.set_ylim(-0.5, n_genes - 0.5)
            ax_tree.invert_xaxis()
            ax_tree.set_yticks([])
            ax_tree.spines['top'].set_visible(False)
            ax_tree.spines['right'].set_visible(False)
            ax_tree.spines['bottom'].set_visible(False)
            ax_tree.set_xticks([])
        except Exception:
            ax_tree.axis('off')

        cmap = LinearSegmentedColormap.from_list('bwr', ['#2166AC', 'white', '#B2182B'])
        im = ax_heat.imshow(scaled.values, aspect='auto', cmap=cmap,
                            interpolation='nearest')
        ax_heat.set_xticks(range(len(scaled.columns)))
        ax_heat.set_xticklabels(scaled.columns, rotation=45, ha='right', fontsize=8)
        ax_heat.set_yticks(range(n_genes))
        fs = 5 if n_genes > 40 else (7 if n_genes > 20 else 9)
        ax_heat.set_yticklabels(scaled.index, fontsize=fs)
        cb = plt.colorbar(im, cax=ax_cbar, label='log2(count+1)')
        max_val = scaled.values.max()
        cb.set_ticks(np.linspace(0, np.ceil(max_val), 5))
        fig.suptitle(title, fontsize=12, fontweight='bold', y=1.01)

        plt.tight_layout()
        out = output_dir / f"heatmap_phylogeny_{family}.png"
        fig.savefig(out, dpi=300, bbox_inches='tight')
        fig.savefig(output_dir / f"heatmap_phylogeny_{family}.pdf", dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"    Saved: {out}")
    else:
        # Fallback: clustered heatmap without tree
        generate_heatmap(norm_counts, overlap, output_dir,
                         filename=f"heatmap_phylogeny_{family}",
                         title=title, show_labels=True, fontsize=6, scale_rows=False)


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 13: Fetch gene annotations from NCBI (Python Entrez = R's rentrez)
# ═══════════════════════════════════════════════════════════════════════════════
# R equivalent:
#   library(rentrez)
#   gene_records <- entrez_summary(db="gene", id=gene_ids)
#   ... extract genomicinfo, name, otheraliases, description

def fetch_ncbi_gene_info(gene_ids, output_dir, family="P450", email="daisycortesj@vt.edu"):
    """Fetch gene info from NCBI (Python equivalent of R's rentrez)."""
    print(f"\n  Fetching NCBI gene info for {family}...")

    try:
        from Bio import Entrez
    except ImportError:
        print("    Biopython not installed — skipping NCBI fetch")
        print("    Install: pip install biopython")
        return pd.DataFrame()

    Entrez.email = email

    numeric_ids = []
    for gid in gene_ids:
        num = str(gid).replace('LOC', '')
        if num.isdigit():
            numeric_ids.append(num)

    if not numeric_ids:
        print("    No valid numeric gene IDs found")
        return pd.DataFrame()

    records = []
    batch_size = 100
    for i in range(0, len(numeric_ids), batch_size):
        batch = numeric_ids[i:i + batch_size]
        try:
            handle = Entrez.esummary(db="gene", id=",".join(batch), retmax=batch_size)
            result = Entrez.read(handle)
            handle.close()

            docs = result.get('DocumentSummarySet', {}).get('DocumentSummary', [])
            for doc in docs:
                uid = str(doc.attributes.get('uid', ''))
                genomic = doc.get('GenomicInfo', [])
                records.append({
                    'GeneID': uid,
                    'gene_id': f"LOC{uid}",
                    'Name': doc.get('Name', ''),
                    'OtherAliases': doc.get('OtherAliases', ''),
                    'Description': doc.get('Description', ''),
                    'ChrLoc': genomic[0].get('ChrLoc', '') if genomic else '',
                    'chrstart': int(genomic[0].get('ChrStart', 0)) if genomic else 0,
                    'chrstop': int(genomic[0].get('ChrStop', 0)) if genomic else 0,
                    'ChrAccVer': genomic[0].get('ChrAccVer', '') if genomic else '',
                })
        except Exception as e:
            print(f"    NCBI batch {i // batch_size + 1} failed: {e}")

    df = pd.DataFrame(records)
    if not df.empty:
        out = output_dir / f"{family}_genomic_locations_description.csv"
        df.to_csv(out, index=False)
        print(f"    Saved: {out} ({len(df)} records)")

    return df


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 14: Genomic proximity / cluster analysis
# ═══════════════════════════════════════════════════════════════════════════════
# R equivalent:
#   gene_info_df$Distance <- with(gene_info_df, abs(chrstart - chrstart[1]))
#   threshold_distance <- 50000
#   close_genes <- gene_info_df[gene_info_df$Distance < threshold_distance, ]
#   plot(gene_info_df$Name, gene_info_df$Distance, ...)

def genomic_proximity_analysis(gene_info_df, output_dir, family="P450",
                               threshold=50000, gtf_path=None, gene_ids=None):
    """Genomic proximity analysis with distance plots."""
    print(f"\n  Genomic proximity analysis for {family}...")

    # If NCBI fetch didn't work, try GTF
    if gene_info_df.empty and gtf_path and gene_ids:
        print("    Falling back to GTF for coordinates...")
        gene_ids_set = set(gene_ids)
        records = []
        with open(gtf_path) as f:
            for line in f:
                if line.startswith('#') or '\tgene\t' not in line:
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                attrs = {}
                for item in fields[8].split(';'):
                    item = item.strip()
                    if item:
                        parts = item.split(' ', 1)
                        if len(parts) == 2:
                            attrs[parts[0]] = parts[1].strip('"')
                gid = attrs.get('gene_id', '')
                if gid in gene_ids_set:
                    records.append({
                        'gene_id': gid,
                        'Name': attrs.get('gene', gid),
                        'ChrAccVer': fields[0],
                        'chrstart': int(fields[3]),
                        'chrstop': int(fields[4]),
                    })
        gene_info_df = pd.DataFrame(records)

    if gene_info_df.empty or 'chrstart' not in gene_info_df.columns:
        print("    No coordinate data — skipping proximity analysis")
        return

    # Group by chromosome
    chrom_col = 'ChrAccVer' if 'ChrAccVer' in gene_info_df.columns else 'ChrLoc'
    if chrom_col not in gene_info_df.columns:
        print("    No chromosome column found — skipping")
        return

    for chrom, group in gene_info_df.groupby(chrom_col):
        if not chrom or len(group) < 2:
            continue

        group = group.sort_values('chrstart').reset_index(drop=True)
        anchor_pos = group.loc[0, 'chrstart']
        group = group.copy()
        group['Distance'] = (group['chrstart'] - anchor_pos).abs()

        close_genes = group[group['Distance'] < threshold]
        name_col = 'Name' if 'Name' in group.columns else 'gene_id'

        fig, ax = plt.subplots(figsize=(max(10, len(group) * 0.5), 6))
        x = range(len(group))
        colors = ['blue' if d < threshold else 'gray' for d in group['Distance']]

        ax.bar(x, group['Distance'], color=colors, alpha=0.7, edgecolor='black', linewidth=0.5)
        ax.axhline(threshold, color='red', linestyle='--', linewidth=1.5,
                    label=f'Threshold ({threshold // 1000} kb)')
        ax.set_xticks(x)
        labels = group[name_col].values
        ax.set_xticklabels(labels, rotation=90, fontsize=6)
        ax.set_ylabel('Distance from anchor (bp)')
        ax.set_title(f'{family} — Distance from {labels[0]} ({chrom})')
        ax.legend()

        plt.tight_layout()
        safe = chrom.replace('.', '_')
        out = output_dir / f"{family}_distance_{safe}.png"
        fig.savefig(out, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"    Saved: {out} ({len(close_genes)} genes within {threshold // 1000} kb)")

    gene_info_df_out = output_dir / f"Genes_with_distance_{family}.csv"
    gene_info_df.to_csv(gene_info_df_out, index=False)
    print(f"    Saved: {gene_info_df_out}")


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN — run everything
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Replicate the student's R DESeq2 analysis in Python",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--counts", default=DEFAULTS["counts"],
                        help="Count matrix TSV (genes × samples)")
    parser.add_argument("--metadata", default=DEFAULTS["metadata"],
                        help="Sample metadata TSV")
    parser.add_argument("--p450-list", default=DEFAULTS["p450_list"],
                        help="P450 gene list (one ID per line)")
    parser.add_argument("--mt-list", default=DEFAULTS["mt_list"],
                        help="Methyltransferase gene list")
    parser.add_argument("--gtf", default=DEFAULTS["gtf"],
                        help="GTF annotation file")
    parser.add_argument("--protein-fasta", default=DEFAULTS["protein_fasta"],
                        help="Full protein FASTA")
    parser.add_argument("-o", "--output-dir", default=DEFAULTS["output_dir"],
                        help="Output directory")
    parser.add_argument("--contrast-factor", default="condition",
                        help="Metadata column for contrast (default: condition)")
    parser.add_argument("--contrast-A", default="R",
                        help="Numerator condition (default: R = root)")
    parser.add_argument("--contrast-B", default="L",
                        help="Denominator/baseline condition (default: L = leaf)")
    parser.add_argument("--min-counts", type=int, default=20,
                        help="Min total counts to keep gene (default: 20)")
    parser.add_argument("--padj", type=float, default=0.05,
                        help="Adjusted p-value cutoff (default: 0.05)")
    parser.add_argument("--lfc", type=float, default=2.0,
                        help="Log2 fold change cutoff for up/down categorization and plots (default: 2.0)")
    parser.add_argument("--threshold", type=int, default=50000,
                        help="Genomic cluster threshold in bp (default: 50000)")
    parser.add_argument("--email", default="daisycortesj@vt.edu",
                        help="Email for NCBI Entrez queries")
    parser.add_argument("--skip-phylo", action="store_true",
                        help="Skip phylogenetic tree building (MAFFT/FastTree)")
    parser.add_argument("--skip-ncbi", action="store_true",
                        help="Skip NCBI gene info fetching")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Validate inputs
    for name, path in [("counts", args.counts), ("metadata", args.metadata)]:
        if not Path(path).exists():
            print(f"ERROR: {name} file not found: {path}")
            print(f"Run the upstream pipeline step first.")
            sys.exit(1)

    print("╔══════════════════════════════════════════════════════════════╗")
    print("║  Student R Analysis → Python Replication                    ║")
    print("╚══════════════════════════════════════════════════════════════╝")
    print(f"  Counts:   {args.counts}")
    print(f"  Metadata: {args.metadata}")
    print(f"  Contrast: {args.contrast_A} vs {args.contrast_B}")
    print(f"  Cutoffs:  padj < {args.padj} (significance), |log2FC| >= {args.lfc} (up/down & plots)")
    print(f"  Output:   {output_dir}")
    print()

    # ── Steps 0-1: Load + filter ──
    counts, metadata = load_count_matrix(args.counts, args.metadata, args.min_counts)

    # ── Steps 2-4: PyDESeq2 ──
    dds, results_df, norm_counts = run_pydeseq2(
        counts, metadata, args.contrast_factor, args.contrast_A, args.contrast_B,
        output_dir
    )

    # ── Step 5: QC plots ──
    generate_pca_plot(norm_counts, metadata, args.contrast_factor, output_dir)
    generate_dispersion_plot(dds, output_dir)

    # ── Step 7: Filter significant ──
    sig_all, sig_up, sig_down = filter_significant(
        results_df, args.padj, args.lfc, output_dir
    )

    # ── Step 8: Heatmap of all significant genes ──
    print("\n" + "=" * 70)
    print("STEP 8: Heatmap — all significant genes")
    print("=" * 70)
    generate_heatmap(norm_counts, sig_all.index, output_dir,
                     "heatmap_all_significant_genes",
                     f"All Significant DE Genes (n={len(sig_all)})",
                     show_labels=False, scale_rows=False)

    # ── Step 9: Volcano + MA ──
    generate_volcano_plot(results_df, output_dir, args.padj, args.lfc,
                          args.contrast_A, args.contrast_B)
    generate_ma_plot(results_df, output_dir, args.padj, args.lfc)

    # ── Steps 10-11: Gene family filtering + heatmaps ──
    families = {}

    if args.p450_list and Path(args.p450_list).exists():
        print("\n" + "=" * 70)
        print("STEP 10: P450 gene family analysis")
        print("=" * 70)
        p450_ids = load_gene_list(args.p450_list)
        print(f"  P450 list: {len(p450_ids)} genes from {args.p450_list}")
        p450_de = filter_gene_family(results_df, p450_ids, "P450", output_dir)
        families["P450"] = p450_ids

        # R equivalent: filtered_data_P450 <- res.df[res.df$gene_id %in% P450_list_refseq$V1, ]
        # This is the DE results table filtered to only P450 genes (no counts)
        filtered_data_p450 = results_df.loc[results_df.index.isin(p450_ids)]
        filtered_p450_file = output_dir / "filtered_data_P450.csv"
        filtered_data_p450.to_csv(filtered_p450_file)
        print(f"  Saved: {filtered_p450_file} ({len(filtered_data_p450)} P450 genes)")

        # R equivalent: P450 <- merge(normalized_counts, filtered_data_P450, by=0)
        # Separate merged file with counts + DE stats side by side
        p450_in_counts = sorted(set(p450_ids) & set(norm_counts.index))
        if p450_in_counts:
            p450_merged = norm_counts.loc[p450_in_counts].join(
                filtered_data_p450[['baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj']],
                how='left'
            )
            merged_file = output_dir / "P450_merged_counts_results.csv"
            p450_merged.to_csv(merged_file)
            print(f"  Saved: {merged_file} ({len(p450_in_counts)} genes)")

        # R equivalent: pheatmap(log2(P450_counts+1))  — legend shows 0-10 log2 values
        generate_heatmap(norm_counts, p450_ids, output_dir,
                         "heatmap_P450", f"P450 Genes (n={len(p450_ids)})",
                         show_labels=True, fontsize=4, scale_rows=False)

        # Also generate z-scored version: pheatmap(log2(P450_counts+1), scale='row')
        generate_heatmap(norm_counts, p450_ids, output_dir,
                         "heatmap_P450_zscore", f"P450 Genes z-scored (n={len(p450_ids)})",
                         show_labels=True, fontsize=4, scale_rows=True)

        # Upregulated P450s only (from this analysis)
        p450_up = [g for g in p450_ids if g in sig_up.index]
        if p450_up:
            generate_heatmap(norm_counts, p450_up, output_dir,
                             "heatmap_P450_upregulated",
                             f"P450 Upregulated (n={len(p450_up)})",
                             show_labels=True, fontsize=5, scale_rows=False)

            generate_heatmap(norm_counts, p450_up, output_dir,
                             "heatmap_P450_upregulated_zscore",
                             f"P450 Upregulated z-scored (n={len(p450_up)})",
                             show_labels=True, fontsize=5, scale_rows=True)

        # Pre-computed P450 subsets from Geneious curation (auto-detected)
        student_dir = Path(args.p450_list).parent
        prev_subsets = [
            ("P450_DE_geneious",
             "P450_expression_refseq_logfold2.txt",
             "P450 DE Genes — Geneious"),
            ("P450_upregulated_geneious",
             "P450_upregulated.txt",
             "P450 Upregulated — Geneious"),
            ("P450_downregulated_geneious",
             "P450_downregulated_list.txt",
             "P450 Downregulated — Geneious"),
        ]
        found_any = False
        for tag, filename, title_prefix in prev_subsets:
            subset_path = student_dir / filename
            if not subset_path.exists():
                continue
            if not found_any:
                print(f"\n  Geneious P450 subsets (from {student_dir.name}/):")
                found_any = True
            raw = load_gene_list(subset_path)
            ids = [f"LOC{g}" if g.isdigit() else g for g in raw]
            overlap = [g for g in ids if g in norm_counts.index]
            if overlap:
                generate_heatmap(norm_counts, overlap, output_dir,
                                 f"heatmap_{tag}",
                                 f"{title_prefix} (n={len(overlap)})",
                                 show_labels=True, fontsize=5, scale_rows=False)
                print(f"    {filename}: {len(overlap)} genes → heatmap_{tag}.png")
            else:
                print(f"    {filename}: 0 genes overlapping count matrix")
    else:
        print(f"\n  P450 list not found: {args.p450_list} — skipping P450 steps")

    if args.mt_list and Path(args.mt_list).exists():
        print("\n" + "=" * 70)
        print("STEP 11: Methyltransferase gene family analysis")
        print("=" * 70)
        mt_ids = load_gene_list(args.mt_list)
        print(f"  MT list: {len(mt_ids)} genes from {args.mt_list}")
        mt_de = filter_gene_family(results_df, mt_ids, "MT", output_dir)
        families["OMT"] = mt_ids

        generate_heatmap(norm_counts, mt_ids, output_dir,
                         "heatmap_MT", f"Methyltransferase Genes (n={len(mt_ids)})",
                         show_labels=True, fontsize=6, scale_rows=False)
    else:
        print(f"\n  MT list not provided or not found — skipping MT steps")

    # ── Step 12: Phylogenetic tree + heatmap ──
    if not args.skip_phylo:
        print("\n" + "=" * 70)
        print("STEP 12: Phylogenetic tree + heatmap (R: ggtree + gheatmap)")
        print("=" * 70)

        for fam_name, fam_ids in families.items():
            tree_path = build_tree_if_needed(
                args.protein_fasta, set(fam_ids), args.gtf, output_dir, fam_name
            )
            generate_phylo_heatmap(tree_path, norm_counts, fam_ids, output_dir,
                                   family=fam_name,
                                   title=f"{fam_name} Phylogenetic Tree + Expression")

            up_ids = [g for g in fam_ids if g in sig_up.index]
            if len(up_ids) >= 3:
                generate_phylo_heatmap(tree_path, norm_counts, up_ids, output_dir,
                                       family=f"{fam_name}_upregulated",
                                       title=f"{fam_name} Upregulated — Phylo + Expression")

    # ── Step 13: NCBI gene info (Python's Entrez = R's rentrez) ──
    ncbi_results = {}
    if not args.skip_ncbi:
        print("\n" + "=" * 70)
        print("STEP 13: NCBI gene info (R equivalent: rentrez)")
        print("=" * 70)

        for fam_name, fam_ids in families.items():
            df = fetch_ncbi_gene_info(fam_ids, output_dir, fam_name, args.email)
            ncbi_results[fam_name] = df

        # Also fetch for downregulated genes
        down_ids = list(sig_down.index)
        if down_ids:
            df = fetch_ncbi_gene_info(down_ids, output_dir, "downregulated", args.email)
            ncbi_results["downregulated"] = df

    # ── Step 14: Genomic proximity analysis ──
    print("\n" + "=" * 70)
    print("STEP 14: Genomic proximity / cluster analysis")
    print("=" * 70)

    for fam_name, fam_ids in families.items():
        gene_df = ncbi_results.get(fam_name, pd.DataFrame())
        genomic_proximity_analysis(gene_df, output_dir, fam_name,
                                   args.threshold, args.gtf, fam_ids)

    # Downregulated genes proximity (matches R's final distance analysis section)
    down_gene_df = ncbi_results.get("downregulated", pd.DataFrame())
    down_ids = list(sig_down.index)
    if not down_gene_df.empty or down_ids:
        genomic_proximity_analysis(down_gene_df, output_dir, "downregulated",
                                   args.threshold, args.gtf, down_ids)

    # ── Summary ──
    print("\n" + "=" * 70)
    print("COMPLETE — All steps finished")
    print("=" * 70)
    print(f"\n  Output directory: {output_dir}")
    print(f"\n  Files created:")
    for f in sorted(output_dir.iterdir()):
        size = f.stat().st_size
        unit = "KB" if size > 1024 else "B"
        val = size / 1024 if size > 1024 else size
        print(f"    {f.name:50s} {val:>8.1f} {unit}")

    print(f"\n  ── Step-by-step mapping ──")
    print(f"  R step                     Python output")
    print(f"  ───────────────────────────────────────────────────")
    print(f"  normalized_counts.csv   →  normalized_counts_refseq.csv")
    print(f"  res_all.csv             →  res_all.csv")
    print(f"  res_all_significant.csv →  res_all_significant.csv")
    print(f"  res_Upregulated.csv     →  res_Upregulated.csv")
    print(f"  heatmap_all_sig.png     →  heatmap_all_significant_genes.png")
    print(f"  volcano_plot.png        →  volcano_plot.png")
    print(f"  plotMA                  →  ma_plot.png")
    print(f"  heatmap_P450.png        →  heatmap_P450.png")
    print(f"  heatmap_MT.png          →  heatmap_MT.png")
    print(f"  ggtree + gheatmap       →  heatmap_phylogeny_P450.png")
    print(f"  rentrez gene info       →  P450_genomic_locations_description.csv")
    print(f"  distance analysis       →  Genes_with_distance_P450.csv")


if __name__ == "__main__":
    main()
