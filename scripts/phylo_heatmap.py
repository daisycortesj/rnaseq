#!/usr/bin/env python3
"""
Phylogenetic Tree + Expression Heatmap — combined visualization.

Python equivalent of the previous student's R workflow using ggtree + gheatmap.
Builds a phylogenetic tree from protein sequences (MAFFT → FastTree) and
combines it with a gene expression heatmap side by side.

Three modes:
  1. Build tree from protein FASTA (runs MAFFT + FastTree)
  2. Load a pre-built Newick tree file
  3. Use hierarchical clustering of expression data (no sequences needed)

─── WHERE DO THE INPUT FILES COME FROM? ────────────────────────────────────────

  File                        How it was made                              HPC location
  ─────────────────────────── ──────────────────────────────────────────── ──────────────────────────────────────────────
  gene_count_matrix.tsv       sbatch scripts/featurecounts.sbatch DC       03_count_tables/00_1_DC/
  sample_metadata.tsv         Created by build_count_matrix.py             03_count_tables/00_1_DC/
  cyp_expressed_list.tsv      sbatch scripts/run_cyp_express_extract.sbatch 07_NRdatabase/cyp450_database/
  cyp_proteins.fasta          sbatch scripts/run_cyp_express_extract.sbatch 07_NRdatabase/cyp450_database/
  geneious_expressed_list.tsv PREV_LIST=...P450_list_RefSeq.txt sbatch     07_NRdatabase/cyp450_database/
                              scripts/run_cyp_express_extract.sbatch
  geneious_proteins.fasta     same as above                                07_NRdatabase/cyp450_database/
  omt_expressed_list.tsv      (same process but for OMT family — not yet)  07_NRdatabase/omt_database/
  omt_proteins.fasta          (same process but for OMT family — not yet)  07_NRdatabase/omt_database/

  Prerequisite chain (run in order if not done yet):
    1. featurecounts.sbatch     → gene_count_matrix.tsv + sample_metadata.tsv
    2. run_pydeseq2_step1_analysis.sbatch DC → pydeseq2_results_UNFILTERED.tsv
    3. run_cyp_express_extract.sbatch       → cyp_expressed_list.tsv + cyp_proteins.fasta

  Required Python packages:
    numpy, pandas, matplotlib, seaborn, scipy
    Optional: biopython (for reading Newick trees via Bio.Phylo)

  External tools (for tree building mode only):
    mafft     — on HPC: module load MAFFT
    FastTree  — on HPC: module load FastTree

────────────────────────────────────────────────────────────────────────────────

Usage — build tree from sequences + show expression:
  python scripts/phylo_heatmap.py \\
      --fasta cyp_proteins.fasta \\
      --expression cyp_expressed_list.tsv \\
      --counts gene_count_matrix.tsv \\
      --metadata sample_metadata.tsv \\
      --output-dir 06_analysis/phylo_heatmap_cyp/

Usage — load existing Newick tree:
  python scripts/phylo_heatmap.py \\
      --tree my_tree.nwk \\
      --expression cyp_expressed_list.tsv \\
      --counts gene_count_matrix.tsv \\
      --metadata sample_metadata.tsv \\
      --output-dir phylo_output/

Usage — expression-based clustering (no tree/sequences):
  python scripts/phylo_heatmap.py \\
      --expression cyp_expressed_list.tsv \\
      --counts gene_count_matrix.tsv \\
      --metadata sample_metadata.tsv \\
      --output-dir phylo_output/ \\
      --cluster-only
"""

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    from matplotlib.colors import LinearSegmentedColormap
    from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
    from scipy.spatial.distance import pdist
except ImportError as e:
    print(f"ERROR: Missing required package: {e}")
    sys.exit(1)


# ── Tree building ────────────────────────────────────────────────────────────

def run_mafft(fasta_path, output_path, threads=4):
    """Run MAFFT multiple sequence alignment."""
    print(f"  Running MAFFT alignment...")
    cmd = ['mafft', '--auto', '--thread', str(threads), str(fasta_path)]
    try:
        with open(output_path, 'w') as out:
            result = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE,
                                    check=True, text=True)
        n_seqs = sum(1 for line in open(output_path) if line.startswith('>'))
        print(f"  MAFFT: aligned {n_seqs} sequences → {output_path}")
        return True
    except FileNotFoundError:
        print("  ERROR: MAFFT not found. Install: conda install -c bioconda mafft")
        return False
    except subprocess.CalledProcessError as e:
        print(f"  ERROR: MAFFT failed: {e.stderr[:500]}")
        return False


def run_fasttree(alignment_path, output_path):
    """Run FastTree to build phylogenetic tree from alignment."""
    print(f"  Running FastTree...")
    cmd = ['FastTree', str(alignment_path)]
    try:
        with open(output_path, 'w') as out:
            result = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE,
                                    check=True, text=True)
        print(f"  FastTree: tree saved → {output_path}")
        return True
    except FileNotFoundError:
        try:
            cmd[0] = 'fasttree'
            with open(output_path, 'w') as out:
                subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE,
                               check=True, text=True)
            print(f"  fasttree: tree saved → {output_path}")
            return True
        except FileNotFoundError:
            print("  ERROR: FastTree not found. Install: conda install -c bioconda fasttree")
            return False
    except subprocess.CalledProcessError as e:
        print(f"  ERROR: FastTree failed: {e.stderr[:500]}")
        return False


def build_tree(fasta_path, output_dir, threads=4):
    """Build phylogenetic tree: FASTA → MAFFT alignment → FastTree."""
    alignment_path = output_dir / "msa_alignment.fasta"
    tree_path = output_dir / "phylo_tree.nwk"

    if not run_mafft(fasta_path, alignment_path, threads):
        return None, None
    if not run_fasttree(alignment_path, tree_path):
        return None, alignment_path

    return tree_path, alignment_path


# ── Newick tree parsing (minimal, no Biopython dependency) ───────────────────

def parse_newick_simple(newick_str):
    """Parse a Newick string into a tree structure for dendrogram-like display.

    Returns (leaf_order, linkage_matrix_approx) or uses Bio.Phylo if available.
    """
    try:
        from Bio import Phylo
        from io import StringIO
        tree = Phylo.read(StringIO(newick_str), 'newick')
        return tree
    except ImportError:
        pass

    leaf_names = []
    clean = newick_str.strip().rstrip(';')
    current = ''
    for ch in clean:
        if ch in '(,)':
            name = current.strip().split(':')[0].strip()
            if name:
                leaf_names.append(name)
            current = ''
        else:
            current += ch
    name = current.strip().split(':')[0].strip()
    if name:
        leaf_names.append(name)
    return leaf_names


def get_leaf_order_from_tree(tree_obj):
    """Get leaf label order from a Bio.Phylo tree (depth-first)."""
    try:
        from Bio import Phylo
        terminals = tree_obj.get_terminals()
        return [t.name for t in terminals if t.name]
    except Exception:
        if isinstance(tree_obj, list):
            return tree_obj
        return []


def draw_phylo_tree(tree_obj, ax, leaf_order):
    """Draw a phylogenetic tree on a matplotlib axes.

    Uses Bio.Phylo.draw if available, otherwise draws a simple cladogram.
    """
    try:
        from Bio import Phylo

        label_to_y = {name: i for i, name in enumerate(leaf_order)}

        def get_y(clade):
            if clade.is_terminal():
                return label_to_y.get(clade.name, 0)
            children_y = [get_y(c) for c in clade.clades]
            return np.mean(children_y)

        def get_x(clade, current_x=0):
            return current_x + (clade.branch_length or 0)

        def draw_clade(clade, x_start=0):
            x_end = x_start + (clade.branch_length or 0)
            y = get_y(clade)

            ax.plot([x_start, x_end], [y, y], color='black', linewidth=0.8)

            if not clade.is_terminal():
                child_ys = []
                for child in clade.clades:
                    child_y = get_y(child)
                    child_ys.append(child_y)
                    draw_clade(child, x_end)

                ax.plot([x_end, x_end], [min(child_ys), max(child_ys)],
                        color='black', linewidth=0.8)

        draw_clade(tree_obj.root)

        ax.set_ylim(-0.5, len(leaf_order) - 0.5)
        ax.invert_xaxis()
        ax.set_yticks(range(len(leaf_order)))
        ax.set_yticklabels([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.set_xticks([])
        return True

    except Exception as e:
        print(f"  Note: Simplified tree drawing ({e})")
        n = len(leaf_order)
        for i in range(n):
            ax.plot([0, 1], [i, i], color='black', linewidth=0.5)
        ax.set_ylim(-0.5, n - 0.5)
        ax.set_xlim(1.5, -0.5)
        ax.set_yticks(range(n))
        ax.set_yticklabels([])
        ax.axis('off')
        return True


# ── Expression data loading ──────────────────────────────────────────────────

def load_expression_data(expression_path, counts_path, metadata_path, gene_ids=None):
    """Load and prepare expression data for heatmap.

    Returns (log2_counts, sample_conditions, gene_labels) where:
    - log2_counts: DataFrame (genes × samples) of log2(normalized_count + 1)
    - sample_conditions: dict mapping sample → condition
    - gene_labels: list of gene IDs
    """
    counts = pd.read_csv(counts_path, sep='\t', index_col=0)
    metadata = pd.read_csv(metadata_path, sep='\t')

    sample_conditions = {}
    if 'sample' in metadata.columns and 'condition' in metadata.columns:
        sample_conditions = dict(zip(metadata['sample'], metadata['condition']))

    if gene_ids is not None:
        gene_ids_set = set(str(g) for g in gene_ids)
        counts = counts[counts.index.isin(gene_ids_set)]

    if expression_path:
        expr_df = pd.read_csv(expression_path, sep='\t')
        if 'gene_id' not in expr_df.columns:
            expr_df = expr_df.rename(columns={expr_df.columns[0]: 'gene_id'})
        expr_gene_ids = set(expr_df['gene_id'].astype(str))
        counts = counts[counts.index.isin(expr_gene_ids)]

    common_samples = [s for s in counts.columns if s in metadata['sample'].values]
    if not common_samples:
        common_samples = list(counts.columns)
    counts = counts[common_samples]

    counts = counts[(counts.sum(axis=1) > 0)]

    log2_counts = np.log2(counts + 1)

    return log2_counts, sample_conditions, list(counts.index)


# ── Combined plot ────────────────────────────────────────────────────────────

def generate_phylo_heatmap(tree_path_or_obj, log2_counts, output_dir,
                           sample_conditions=None, gene_family="CYP",
                           scale_rows=True, tree_is_object=False):
    """Generate the combined phylogenetic tree + expression heatmap figure.

    This is the Python equivalent of ggtree + gheatmap in R.
    """
    if tree_is_object and tree_path_or_obj is not None:
        tree_obj = tree_path_or_obj
    elif tree_path_or_obj is not None and Path(tree_path_or_obj).exists():
        newick_str = Path(tree_path_or_obj).read_text().strip()
        tree_obj = parse_newick_simple(newick_str)
    else:
        tree_obj = None

    if tree_obj is not None:
        leaf_order = get_leaf_order_from_tree(tree_obj)
        leaf_order = [l for l in leaf_order if l in log2_counts.index]
        remaining = [g for g in log2_counts.index if g not in leaf_order]
        leaf_order = leaf_order + remaining
    else:
        leaf_order = list(log2_counts.index)

    heatmap_data = log2_counts.reindex(leaf_order).dropna(how='all')
    leaf_order = list(heatmap_data.index)
    n_genes = len(leaf_order)
    n_samples = len(heatmap_data.columns)

    if scale_rows:
        row_means = heatmap_data.mean(axis=1)
        row_stds = heatmap_data.std(axis=1).replace(0, 1)
        heatmap_data = heatmap_data.sub(row_means, axis=0).div(row_stds, axis=0)

    height = max(6, n_genes * 0.25)
    tree_width_ratio = 3 if tree_obj is not None else 0.5
    heatmap_width = max(4, n_samples * 0.8)

    fig = plt.figure(figsize=(tree_width_ratio + heatmap_width + 3, height))

    if tree_obj is not None:
        gs = gridspec.GridSpec(1, 3, width_ratios=[tree_width_ratio, heatmap_width, 0.3],
                               wspace=0.05)
        ax_tree = fig.add_subplot(gs[0])
        ax_heat = fig.add_subplot(gs[1])
        ax_cbar = fig.add_subplot(gs[2])
    else:
        gs = gridspec.GridSpec(1, 2, width_ratios=[heatmap_width, 0.3], wspace=0.05)
        ax_tree = None
        ax_heat = fig.add_subplot(gs[0])
        ax_cbar = fig.add_subplot(gs[1])

    cmap = LinearSegmentedColormap.from_list('bwr', ['#2166AC', 'white', '#B2182B'])

    vmax = max(abs(heatmap_data.values.min()), abs(heatmap_data.values.max()))
    if np.isnan(vmax) or vmax == 0:
        vmax = 1

    im = ax_heat.imshow(heatmap_data.values, aspect='auto', cmap=cmap,
                         vmin=-vmax, vmax=vmax, interpolation='nearest')

    ax_heat.set_xticks(range(n_samples))
    sample_labels = list(heatmap_data.columns)
    ax_heat.set_xticklabels(sample_labels, rotation=45, ha='right', fontsize=8)
    ax_heat.set_yticks(range(n_genes))

    gene_labels = leaf_order
    fontsize = 6 if n_genes > 40 else (8 if n_genes > 20 else 10)
    ax_heat.set_yticklabels(gene_labels, fontsize=fontsize)

    if sample_conditions:
        cond_colors = {}
        palette = ['#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00']
        for i, cond in enumerate(sorted(set(sample_conditions.values()))):
            cond_colors[cond] = palette[i % len(palette)]

        for i, sample in enumerate(sample_labels):
            cond = sample_conditions.get(sample, '')
            if cond in cond_colors:
                ax_heat.plot(i, -1.2, 's', color=cond_colors[cond], markersize=8,
                             clip_on=False)

        from matplotlib.lines import Line2D
        legend_elements = [Line2D([0], [0], marker='s', color='w',
                                  markerfacecolor=c, markersize=8, label=k)
                           for k, c in cond_colors.items()]
        ax_heat.legend(handles=legend_elements, loc='upper center',
                       bbox_to_anchor=(0.5, -0.15), ncol=len(cond_colors),
                       fontsize=8, frameon=False)

    if ax_tree is not None and tree_obj is not None:
        try:
            from Bio import Phylo
            draw_phylo_tree(tree_obj, ax_tree, leaf_order)
        except ImportError:
            ax_tree.axis('off')
            ax_tree.text(0.5, 0.5, '(tree display\nrequires\nBiopython)',
                         ha='center', va='center', fontsize=8, style='italic',
                         transform=ax_tree.transAxes)

    plt.colorbar(im, cax=ax_cbar, label='Row-scaled log2(count+1)' if scale_rows
                 else 'log2(count+1)')

    title = f'{gene_family} Phylogenetic Tree + Expression Heatmap'
    if tree_obj is None:
        title = f'{gene_family} Expression Heatmap'
    fig.suptitle(title, fontsize=12, fontweight='bold', y=1.02)

    plt.tight_layout()
    out_pdf = output_dir / f"{gene_family.lower()}_phylo_heatmap.pdf"
    out_png = output_dir / f"{gene_family.lower()}_phylo_heatmap.png"
    fig.savefig(out_pdf, dpi=300, bbox_inches='tight')
    fig.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out_pdf}")
    print(f"  Saved: {out_png}")


def generate_clustered_heatmap(log2_counts, output_dir, sample_conditions=None,
                               gene_family="CYP", scale_rows=True):
    """Generate a seaborn clustermap with hierarchical clustering dendrogram.

    This is the fallback when no tree/sequences are available — clusters
    genes by expression similarity rather than sequence phylogeny.
    """
    try:
        import seaborn as sns
    except ImportError:
        print("  WARNING: seaborn not available for clustermap, skipping")
        return

    data = log2_counts.copy()
    if scale_rows:
        row_means = data.mean(axis=1)
        row_stds = data.std(axis=1).replace(0, 1)
        data = data.sub(row_means, axis=0).div(row_stds, axis=0)

    cmap = LinearSegmentedColormap.from_list('bwr', ['#2166AC', 'white', '#B2182B'])

    n_genes = len(data)
    fontsize = 5 if n_genes > 50 else (7 if n_genes > 25 else 9)

    col_colors = None
    if sample_conditions:
        palette = ['#E41A1C', '#377EB8', '#4DAF4A', '#984EA3']
        cond_list = sorted(set(sample_conditions.values()))
        cond_cmap = {c: palette[i % len(palette)] for i, c in enumerate(cond_list)}
        col_colors = [cond_cmap.get(sample_conditions.get(s, ''), 'gray')
                      for s in data.columns]

    g = sns.clustermap(data, cmap=cmap, figsize=(max(6, len(data.columns) * 0.8),
                       max(8, n_genes * 0.2)),
                       row_cluster=True, col_cluster=False,
                       yticklabels=True, xticklabels=True,
                       col_colors=col_colors,
                       dendrogram_ratio=(0.15, 0.02),
                       cbar_kws={'label': 'Row-scaled log2(count+1)' if scale_rows
                                 else 'log2(count+1)'},
                       linewidths=0)

    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=fontsize)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), fontsize=8,
                                  rotation=45, ha='right')
    g.fig.suptitle(f'{gene_family} Expression Heatmap (clustered by expression)',
                    fontsize=12, fontweight='bold', y=1.02)

    out_pdf = output_dir / f"{gene_family.lower()}_clustered_heatmap.pdf"
    out_png = output_dir / f"{gene_family.lower()}_clustered_heatmap.png"
    g.savefig(out_pdf, dpi=300, bbox_inches='tight')
    g.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close('all')
    print(f"  Saved: {out_pdf}")
    print(f"  Saved: {out_png}")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Phylogenetic tree + expression heatmap (Python ggtree equivalent)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    input_group = parser.add_argument_group("Input sources (pick one tree method)")
    input_group.add_argument("--fasta", default=None,
                             help="Protein FASTA → builds tree with MAFFT + FastTree")
    input_group.add_argument("--tree", default=None,
                             help="Pre-built Newick tree file")
    input_group.add_argument("--cluster-only", action="store_true",
                             help="Skip tree, cluster by expression similarity")

    data_group = parser.add_argument_group("Expression data")
    data_group.add_argument("--expression", default=None,
                            help="Gene family expressed list TSV (filters gene set)")
    data_group.add_argument("--counts", required=True,
                            help="Gene count matrix TSV (genes × samples)")
    data_group.add_argument("--metadata", required=True,
                            help="Sample metadata TSV (sample, condition columns)")

    out_group = parser.add_argument_group("Output")
    out_group.add_argument("--output-dir", "-o", required=True,
                           help="Output directory")
    out_group.add_argument("--family", default="CYP",
                           help="Gene family label for titles (default: CYP)")

    opt_group = parser.add_argument_group("Options")
    opt_group.add_argument("--no-scale", action="store_true",
                           help="Don't row-scale the heatmap (default: scale rows)")
    opt_group.add_argument("--threads", type=int, default=4,
                           help="CPU threads for MAFFT (default: 4)")
    opt_group.add_argument("--subset", default=None,
                           help="Comma-separated: 'significant','upregulated','downregulated' "
                                "to filter genes by direction column in expression file")

    args = parser.parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print(f"Phylo + Heatmap — {args.family}")
    print("=" * 60)

    # ── Load expression data ──
    gene_ids = None
    if args.expression:
        expr_df = pd.read_csv(args.expression, sep='\t')
        if 'gene_id' not in expr_df.columns:
            expr_df = expr_df.rename(columns={expr_df.columns[0]: 'gene_id'})

        if args.subset and 'direction' in expr_df.columns:
            keep = set()
            for s in args.subset.split(','):
                s = s.strip().lower()
                if s == 'significant':
                    keep.update(['root_up', 'leaf_up'])
                elif s == 'upregulated':
                    keep.add('root_up')
                elif s == 'downregulated':
                    keep.add('leaf_up')
                else:
                    keep.add(s)
            expr_df = expr_df[expr_df['direction'].isin(keep)]
            print(f"  Subset filter: {args.subset} → {len(expr_df)} genes")

        gene_ids = expr_df['gene_id'].astype(str).tolist()
        print(f"  Expression list: {args.expression} ({len(gene_ids)} genes)")

    log2_counts, sample_conditions, found_genes = load_expression_data(
        args.expression, args.counts, args.metadata, gene_ids
    )
    print(f"  Genes in heatmap: {len(log2_counts)}")
    print(f"  Samples: {len(log2_counts.columns)}")

    if log2_counts.empty:
        print("ERROR: No genes found in count matrix")
        sys.exit(1)

    scale = not args.no_scale

    # ── Build or load tree ──
    tree_path = None
    if args.fasta and not args.cluster_only:
        fasta_path = Path(args.fasta)
        if not fasta_path.exists():
            print(f"ERROR: FASTA not found: {fasta_path}")
            sys.exit(1)
        print(f"\n  Building phylogenetic tree from: {fasta_path}")
        tree_path, alignment_path = build_tree(fasta_path, output_dir, args.threads)
        if tree_path is None:
            print("  WARNING: Tree building failed — falling back to clustered heatmap")
            args.cluster_only = True
    elif args.tree and not args.cluster_only:
        tree_path = Path(args.tree)
        if not tree_path.exists():
            print(f"ERROR: Tree file not found: {tree_path}")
            sys.exit(1)
        print(f"  Using pre-built tree: {tree_path}")

    # ── Generate plots ──
    print(f"\n  Generating plots...")

    if not args.cluster_only and tree_path is not None:
        generate_phylo_heatmap(tree_path, log2_counts, output_dir,
                               sample_conditions, args.family, scale)

    generate_clustered_heatmap(log2_counts, output_dir, sample_conditions,
                               args.family, scale)

    print()
    print("=" * 60)
    print(f"  Complete. Output: {output_dir}")
    print("=" * 60)


if __name__ == "__main__":
    main()
