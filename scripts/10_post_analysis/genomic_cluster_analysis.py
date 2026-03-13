#!/usr/bin/env python3
"""
Genomic Cluster Analysis — find gene clusters on chromosomes.

Replaces the previous student's R workflow that used rentrez to fetch gene
coordinates from NCBI. Instead, this parses coordinates directly from the
GTF annotation file (no internet needed), calculates pairwise distances
between genes on the same chromosome, identifies clusters within a
user-defined threshold, and generates plots.

Optionally fetches gene descriptions from NCBI using Biopython's Entrez
module (equivalent of R's rentrez package).

Usage:
  python scripts/genomic_cluster_analysis.py \
      --gene-list 07_NRdatabase/cyp450_database/cyp_expressed_list.tsv \
      --gtf 04_reference/dc_genomic.gtf \
      --output-dir 06_analysis/genomic_clusters_cyp/ \
      --threshold 50000

  With NCBI descriptions (requires internet + Biopython):
  python scripts/genomic_cluster_analysis.py \
      --gene-list cyp_expressed_list.tsv \
      --gtf dc_genomic.gtf \
      --output-dir genomic_clusters/ \
      --fetch-ncbi --email you@example.com
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict

import pandas as pd
import numpy as np

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.collections import LineCollection
except ImportError:
    print("ERROR: matplotlib is required. Install with: pip install matplotlib")
    sys.exit(1)


# ── GTF parsing ──────────────────────────────────────────────────────────────

def parse_gtf_gene_coords(gtf_path, gene_ids):
    """Extract gene coordinates from a GTF file for specified gene IDs.

    Returns DataFrame with: gene_id, chromosome, start, end, strand
    """
    gene_ids = set(gene_ids)
    records = []
    seen = set()

    with open(gtf_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'gene':
                continue

            attrs = _parse_attrs(fields[8])
            gid = attrs.get('gene_id', '')
            if gid not in gene_ids or gid in seen:
                continue

            seen.add(gid)
            records.append({
                'gene_id': gid,
                'chromosome': fields[0],
                'start': int(fields[3]),
                'end': int(fields[4]),
                'strand': fields[6],
                'gene_name': attrs.get('gene', attrs.get('Name', '')),
                'description': attrs.get('description', attrs.get('product', '')),
            })

    df = pd.DataFrame(records)
    if df.empty:
        return df

    df = df.sort_values(['chromosome', 'start']).reset_index(drop=True)

    missing = gene_ids - seen
    if missing:
        print(f"  WARNING: {len(missing)} gene(s) not found in GTF")
        if len(missing) <= 10:
            for g in sorted(missing):
                print(f"    {g}")

    return df


def _parse_attrs(attr_string):
    """Parse GTF/GFF3 attribute string to dict."""
    attrs = {}
    if '=' in attr_string and '"' not in attr_string.split('=')[0]:
        for item in attr_string.split(';'):
            item = item.strip()
            if '=' in item:
                k, v = item.split('=', 1)
                attrs[k] = v
    else:
        for item in attr_string.split(';'):
            item = item.strip()
            if item:
                parts = item.split(' ', 1)
                if len(parts) == 2:
                    attrs[parts[0]] = parts[1].strip('"')
    return attrs


# ── NCBI Entrez fetching (Python equivalent of R's rentrez) ─────────────────

def fetch_ncbi_gene_info(gene_ids, email, batch_size=100):
    """Fetch gene info from NCBI using Biopython Entrez (like R's rentrez).

    For each gene, retrieves: official name, aliases, description,
    chromosome coordinates, and summary.

    gene_ids should be numeric NCBI Gene IDs (e.g. '108192285') or
    LOC-prefixed IDs (the LOC prefix is stripped automatically).
    """
    try:
        from Bio import Entrez
    except ImportError:
        print("  WARNING: Biopython not installed. Skipping NCBI fetch.")
        print("  Install with: pip install biopython")
        return pd.DataFrame()

    Entrez.email = email

    numeric_ids = []
    loc_to_numeric = {}
    for gid in gene_ids:
        gid_str = str(gid)
        if gid_str.startswith('LOC'):
            num = gid_str[3:]
            numeric_ids.append(num)
            loc_to_numeric[gid_str] = num
        else:
            numeric_ids.append(gid_str)
            loc_to_numeric[gid_str] = gid_str

    print(f"  Fetching {len(numeric_ids)} gene records from NCBI...")
    records = []

    for i in range(0, len(numeric_ids), batch_size):
        batch = numeric_ids[i:i + batch_size]
        try:
            handle = Entrez.esummary(db="gene", id=",".join(batch), retmax=batch_size)
            result = Entrez.read(handle)
            handle.close()

            doc_sums = result.get('DocumentSummarySet', {}).get('DocumentSummary', [])
            for doc in doc_sums:
                uid = str(doc.attributes.get('uid', ''))
                genomic = doc.get('GenomicInfo', [])
                chrloc = genomic[0].get('ChrLoc', '') if genomic else ''
                chrstart = genomic[0].get('ChrStart', '') if genomic else ''
                chrstop = genomic[0].get('ChrStop', '') if genomic else ''

                records.append({
                    'ncbi_gene_id': uid,
                    'gene_id': f"LOC{uid}",
                    'ncbi_name': doc.get('Name', ''),
                    'ncbi_aliases': doc.get('OtherAliases', ''),
                    'ncbi_description': doc.get('Description', ''),
                    'ncbi_chromosome': chrloc,
                    'ncbi_chrstart': chrstart,
                    'ncbi_chrstop': chrstop,
                })
        except Exception as e:
            print(f"  WARNING: NCBI fetch failed for batch {i//batch_size + 1}: {e}")

    print(f"  Retrieved {len(records)} records from NCBI")
    return pd.DataFrame(records)


# ── Cluster detection ────────────────────────────────────────────────────────

def find_clusters(coord_df, threshold=50000):
    """Identify gene clusters: groups of genes within `threshold` bp on the same chromosome.

    Uses a sliding window: genes are clustered if consecutive genes
    on the same chromosome are within `threshold` of each other.
    """
    clusters = []
    cluster_id = 0

    for chrom, group in coord_df.groupby('chromosome'):
        group = group.sort_values('start').reset_index(drop=True)
        if len(group) < 2:
            for _, row in group.iterrows():
                clusters.append({**row.to_dict(), 'cluster_id': -1, 'cluster_size': 1})
            continue

        current_cluster = [0]
        for i in range(1, len(group)):
            dist = group.loc[i, 'start'] - group.loc[current_cluster[-1], 'end']
            if dist <= threshold:
                current_cluster.append(i)
            else:
                if len(current_cluster) >= 2:
                    for idx in current_cluster:
                        clusters.append({
                            **group.loc[idx].to_dict(),
                            'cluster_id': cluster_id,
                            'cluster_size': len(current_cluster),
                        })
                    cluster_id += 1
                else:
                    clusters.append({
                        **group.loc[current_cluster[0]].to_dict(),
                        'cluster_id': -1,
                        'cluster_size': 1,
                    })
                current_cluster = [i]

        if len(current_cluster) >= 2:
            for idx in current_cluster:
                clusters.append({
                    **group.loc[idx].to_dict(),
                    'cluster_id': cluster_id,
                    'cluster_size': len(current_cluster),
                })
            cluster_id += 1
        else:
            clusters.append({
                **group.loc[current_cluster[0]].to_dict(),
                'cluster_id': -1,
                'cluster_size': 1,
            })

    result = pd.DataFrame(clusters)
    return result


def compute_distance_matrix(coord_df):
    """Compute pairwise distances between all genes on the same chromosome.

    Returns a DataFrame (gene_id x gene_id) with distances in bp.
    Genes on different chromosomes get NaN.
    """
    n = len(coord_df)
    gene_ids = coord_df['gene_id'].values
    chroms = coord_df['chromosome'].values
    midpoints = ((coord_df['start'] + coord_df['end']) / 2).values

    dist = np.full((n, n), np.nan)
    for i in range(n):
        for j in range(n):
            if chroms[i] == chroms[j]:
                dist[i, j] = abs(midpoints[i] - midpoints[j])

    return pd.DataFrame(dist, index=gene_ids, columns=gene_ids)


# ── Plotting ─────────────────────────────────────────────────────────────────

def plot_chromosome_map(coord_df, cluster_df, output_dir, gene_family="CYP",
                        threshold=50000):
    """Plot gene positions along chromosomes, highlighting clusters."""
    chroms = sorted(coord_df['chromosome'].unique(),
                    key=lambda c: (not c.replace('NC_', '').replace('.', '').isdigit(),
                                   c))

    fig, ax = plt.subplots(figsize=(14, max(6, len(chroms) * 0.5)))

    chrom_to_y = {c: i for i, c in enumerate(chroms)}

    clustered_genes = set()
    if cluster_df is not None and 'cluster_id' in cluster_df.columns:
        clustered_genes = set(cluster_df[cluster_df['cluster_id'] >= 0]['gene_id'])

    for _, row in coord_df.iterrows():
        y = chrom_to_y[row['chromosome']]
        color = 'red' if row['gene_id'] in clustered_genes else 'steelblue'
        ax.barh(y, row['end'] - row['start'], left=row['start'],
                height=0.4, color=color, edgecolor='none', alpha=0.7)

    ax.set_yticks(range(len(chroms)))
    ax.set_yticklabels(chroms, fontsize=7)
    ax.set_xlabel('Genomic position (bp)')
    ax.set_title(f'{gene_family} genes on chromosomes (red = clustered within {threshold/1000:.0f} kb)')
    ax.invert_yaxis()

    legend_elements = [
        mpatches.Patch(color='red', alpha=0.7, label=f'Clustered (≤{threshold/1000:.0f} kb)'),
        mpatches.Patch(color='steelblue', alpha=0.7, label='Singleton'),
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=8)

    plt.tight_layout()
    out = output_dir / f"{gene_family.lower()}_chromosome_map.pdf"
    fig.savefig(out, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out}")

    out_png = output_dir / f"{gene_family.lower()}_chromosome_map.png"
    fig2, ax2 = plt.subplots(figsize=(14, max(6, len(chroms) * 0.5)))
    for _, row in coord_df.iterrows():
        y = chrom_to_y[row['chromosome']]
        color = 'red' if row['gene_id'] in clustered_genes else 'steelblue'
        ax2.barh(y, row['end'] - row['start'], left=row['start'],
                 height=0.4, color=color, edgecolor='none', alpha=0.7)
    ax2.set_yticks(range(len(chroms)))
    ax2.set_yticklabels(chroms, fontsize=7)
    ax2.set_xlabel('Genomic position (bp)')
    ax2.set_title(f'{gene_family} genes on chromosomes (red = clustered within {threshold/1000:.0f} kb)')
    ax2.invert_yaxis()
    ax2.legend(handles=legend_elements, loc='upper right', fontsize=8)
    plt.tight_layout()
    fig2.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close(fig2)
    print(f"  Saved: {out_png}")


def plot_distance_from_anchor(coord_df, output_dir, gene_family="CYP",
                              threshold=50000, anchor_gene=None):
    """Plot distance of each gene from an anchor gene (like the student's
    distance-from-EIGS plot). If no anchor given, uses the first gene."""
    for chrom, group in coord_df.groupby('chromosome'):
        if len(group) < 2:
            continue

        group = group.sort_values('start').reset_index(drop=True)

        if anchor_gene and anchor_gene in group['gene_id'].values:
            anchor_idx = group[group['gene_id'] == anchor_gene].index[0]
        else:
            anchor_idx = 0

        anchor_pos = (group.loc[anchor_idx, 'start'] + group.loc[anchor_idx, 'end']) / 2
        group = group.copy()
        group['distance'] = group.apply(
            lambda r: abs((r['start'] + r['end']) / 2 - anchor_pos), axis=1
        )

        label_col = 'gene_name' if group['gene_name'].any() else 'gene_id'
        labels = group[label_col].values
        if labels[0] == '' or pd.isna(labels[0]):
            labels = group['gene_id'].values

        fig, ax = plt.subplots(figsize=(max(10, len(group) * 0.4), 6))
        x = range(len(group))
        colors = ['blue' if d <= threshold else 'gray' for d in group['distance']]

        ax.bar(x, group['distance'], color=colors, alpha=0.7, edgecolor='black', linewidth=0.5)
        ax.axhline(y=threshold, color='red', linestyle='--', linewidth=1.5,
                    label=f'Threshold ({threshold/1000:.0f} kb)')
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=90, fontsize=6, ha='center')
        ax.set_ylabel('Distance from anchor (bp)')
        anchor_label = group.loc[anchor_idx, label_col] or group.loc[anchor_idx, 'gene_id']
        ax.set_title(f'{gene_family} gene distances from {anchor_label} — {chrom}')
        ax.legend(fontsize=8)

        plt.tight_layout()
        safe_chrom = chrom.replace('.', '_')
        out = output_dir / f"{gene_family.lower()}_distance_{safe_chrom}.pdf"
        fig.savefig(out, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"  Saved: {out}")


def plot_cluster_detail(cluster_df, output_dir, gene_family="CYP"):
    """For each cluster, plot a zoomed-in view of gene positions."""
    if cluster_df is None or 'cluster_id' not in cluster_df.columns:
        return

    clustered = cluster_df[cluster_df['cluster_id'] >= 0]
    if clustered.empty:
        print("  No clusters found — skipping detail plots.")
        return

    for cid, group in clustered.groupby('cluster_id'):
        group = group.sort_values('start').reset_index(drop=True)
        chrom = group['chromosome'].iloc[0]
        label_col = 'gene_name' if group['gene_name'].any() else 'gene_id'

        fig, ax = plt.subplots(figsize=(max(8, len(group) * 1.5), 4))

        for i, (_, row) in enumerate(group.iterrows()):
            color = 'salmon' if row.get('strand', '+') == '+' else 'skyblue'
            ax.barh(0, row['end'] - row['start'], left=row['start'],
                    height=0.4, color=color, edgecolor='black', linewidth=0.8)
            label = row[label_col] if row[label_col] else row['gene_id']
            ax.text((row['start'] + row['end']) / 2, 0.35, label,
                    ha='center', va='bottom', fontsize=7, rotation=45)

        span = group['end'].max() - group['start'].min()
        ax.set_xlim(group['start'].min() - span * 0.1,
                    group['end'].max() + span * 0.1)
        ax.set_ylim(-0.5, 1.0)
        ax.set_yticks([])
        ax.set_xlabel('Genomic position (bp)')
        ax.set_title(f'{gene_family} cluster {cid} — {chrom} ({len(group)} genes, '
                     f'{span/1000:.1f} kb span)')

        legend_elements = [
            mpatches.Patch(color='salmon', label='+ strand'),
            mpatches.Patch(color='skyblue', label='- strand'),
        ]
        ax.legend(handles=legend_elements, fontsize=8)

        plt.tight_layout()
        out = output_dir / f"{gene_family.lower()}_cluster_{cid}_detail.pdf"
        fig.savefig(out, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"  Saved: {out}")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Genomic cluster analysis — find gene clusters on chromosomes",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python scripts/genomic_cluster_analysis.py \\
      --gene-list 07_NRdatabase/cyp450_database/cyp_expressed_list.tsv \\
      --gtf 04_reference/dc_genomic.gtf \\
      --output-dir 06_analysis/genomic_clusters_cyp/

  With NCBI gene descriptions:
  python scripts/genomic_cluster_analysis.py \\
      --gene-list cyp_expressed_list.tsv --gtf dc_genomic.gtf \\
      --output-dir genomic_clusters/ --fetch-ncbi --email you@vt.edu
        """
    )
    parser.add_argument("--gene-list", required=True,
                        help="Gene list TSV/CSV/TXT (must have gene_id column or one ID per line)")
    parser.add_argument("--gtf", required=True,
                        help="GTF/GFF annotation file (e.g. dc_genomic.gtf)")
    parser.add_argument("--output-dir", "-o", required=True,
                        help="Output directory for results and plots")
    parser.add_argument("--threshold", type=int, default=50000,
                        help="Distance threshold in bp for clustering (default: 50000)")
    parser.add_argument("--family", default="CYP",
                        help="Gene family label for plot titles (default: CYP)")
    parser.add_argument("--anchor-gene", default=None,
                        help="Gene ID to use as anchor for distance plot (default: first gene per chrom)")
    parser.add_argument("--fetch-ncbi", action="store_true",
                        help="Fetch gene descriptions from NCBI Entrez (requires internet + Biopython)")
    parser.add_argument("--email", default=None,
                        help="Email for NCBI Entrez queries (required with --fetch-ncbi)")
    parser.add_argument("--deseq", default=None,
                        help="PyDESeq2 results TSV to attach expression stats to output")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print(f"Genomic Cluster Analysis — {args.family}")
    print("=" * 60)

    # ── Load gene list ──
    gene_path = Path(args.gene_list)
    if not gene_path.exists():
        print(f"ERROR: Gene list not found: {gene_path}")
        sys.exit(1)

    if str(gene_path).endswith('.txt'):
        with open(gene_path) as f:
            gene_ids = [line.strip() for line in f if line.strip() and not line.startswith('#')]
    else:
        sep = '\t' if str(gene_path).endswith(('.tsv', '.tab')) else ','
        df_genes = pd.read_csv(gene_path, sep=sep)
        if 'gene_id' not in df_genes.columns:
            df_genes = df_genes.rename(columns={df_genes.columns[0]: 'gene_id'})
        gene_ids = df_genes['gene_id'].astype(str).tolist()

    print(f"  Gene list: {gene_path} ({len(gene_ids)} genes)")

    # ── Parse GTF for coordinates ──
    gtf_path = Path(args.gtf)
    if not gtf_path.exists():
        print(f"ERROR: GTF not found: {gtf_path}")
        sys.exit(1)

    print(f"  GTF: {gtf_path}")
    coord_df = parse_gtf_gene_coords(gtf_path, gene_ids)
    print(f"  Genes with coordinates: {len(coord_df)}")
    print(f"  Chromosomes: {coord_df['chromosome'].nunique()}")

    if coord_df.empty:
        print("ERROR: No gene coordinates found. Check gene_id format.")
        sys.exit(1)

    # ── Optionally fetch NCBI info ──
    if args.fetch_ncbi:
        if not args.email:
            print("ERROR: --email is required with --fetch-ncbi")
            sys.exit(1)
        ncbi_df = fetch_ncbi_gene_info(gene_ids, args.email)
        if not ncbi_df.empty:
            coord_df = coord_df.merge(ncbi_df, on='gene_id', how='left')
            ncbi_out = output_dir / f"{args.family.lower()}_ncbi_gene_info.tsv"
            ncbi_df.to_csv(ncbi_out, sep='\t', index=False)
            print(f"  Saved NCBI info: {ncbi_out}")

    # ── Optionally attach DESeq2 stats ──
    if args.deseq:
        deseq_path = Path(args.deseq)
        if deseq_path.exists():
            deseq_df = pd.read_csv(deseq_path, sep='\t', comment='#')
            if 'gene_id' not in deseq_df.columns:
                deseq_df = deseq_df.rename(columns={deseq_df.columns[0]: 'gene_id'})
            deseq_df['gene_id'] = deseq_df['gene_id'].astype(str)
            expr_cols = ['baseMean', 'log2FoldChange', 'pvalue', 'padj']
            merge_cols = ['gene_id'] + [c for c in expr_cols if c in deseq_df.columns]
            coord_df = coord_df.merge(deseq_df[merge_cols], on='gene_id', how='left')
            print(f"  Attached expression stats from: {deseq_path}")

    # ── Find clusters ──
    print(f"\n  Clustering threshold: {args.threshold/1000:.0f} kb")
    cluster_df = find_clusters(coord_df, threshold=args.threshold)

    n_clustered = len(cluster_df[cluster_df['cluster_id'] >= 0])
    n_clusters = cluster_df[cluster_df['cluster_id'] >= 0]['cluster_id'].nunique()
    n_singletons = len(cluster_df[cluster_df['cluster_id'] < 0])

    print(f"  Clusters found: {n_clusters}")
    print(f"  Genes in clusters: {n_clustered}")
    print(f"  Singleton genes: {n_singletons}")

    # ── Save results ──
    coord_out = output_dir / f"{args.family.lower()}_genomic_coordinates.tsv"
    coord_df.to_csv(coord_out, sep='\t', index=False)
    print(f"\n  Saved coordinates: {coord_out}")

    cluster_out = output_dir / f"{args.family.lower()}_clusters.tsv"
    cluster_df.to_csv(cluster_out, sep='\t', index=False)
    print(f"  Saved clusters: {cluster_out}")

    dist_df = compute_distance_matrix(coord_df)
    dist_out = output_dir / f"{args.family.lower()}_distance_matrix.tsv"
    dist_df.to_csv(dist_out, sep='\t')
    print(f"  Saved distance matrix: {dist_out}")

    # ── Print cluster summary ──
    if n_clusters > 0:
        print(f"\n  Cluster details:")
        for cid, grp in cluster_df[cluster_df['cluster_id'] >= 0].groupby('cluster_id'):
            chrom = grp['chromosome'].iloc[0]
            span = grp['end'].max() - grp['start'].min()
            genes = ', '.join(grp['gene_id'].tolist())
            print(f"    Cluster {cid}: {chrom}, {len(grp)} genes, "
                  f"{span/1000:.1f} kb span")
            print(f"      {genes}")

    # ── Generate plots ──
    print(f"\n  Generating plots...")
    plot_chromosome_map(coord_df, cluster_df, output_dir,
                        gene_family=args.family, threshold=args.threshold)
    plot_distance_from_anchor(coord_df, output_dir, gene_family=args.family,
                              threshold=args.threshold, anchor_gene=args.anchor_gene)
    plot_cluster_detail(cluster_df, output_dir, gene_family=args.family)

    print()
    print("=" * 60)
    print(f"  Complete. Output: {output_dir}")
    print("=" * 60)


if __name__ == "__main__":
    main()
