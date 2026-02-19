#!/usr/bin/env python3
"""
PyDESeq2 Step 3: Generate Plots
================================

Generates publication-quality plots from DESeq2 results:
  - MA plot (red/blue up/down, gene counts, threshold caption)
  - Enhanced Volcano plot (4-color, boxed labels with connectors)
  - PCA plot (top 500 variable genes, 95% confidence ellipses)
  - Sample correlation heatmap (short names: DC1L, DC1R, etc.)
  - CYP heatmap (subfamily sidebar: CYP71D, CYP81, etc., no gene labels)
  - OMT heatmap (subfamily sidebar: COMT, CCoAOMT, etc.)

Accepts any of these as input:
  - Unfiltered DESeq2 results (from Step 1)
  - Filtered DESeq2 results (from Step 2)
  - Combined BLAST+DESeq2 annotated file (from run_combine_filter.sbatch)
  - Gene-family-verified file (from extract_gene_families.py)

CYP and OMT genes are identified automatically:
  - If the input has a 'gene_family' column, that is used directly
  - If the input has a 'blast_description' column, regex patterns identify
    CYP and OMT genes (same patterns as extract_gene_families.py)
  - If neither column exists, gene-family heatmaps are skipped

Usage:
  # From combined annotated file (auto-detects CYP/OMT):
  python pydeseq2_generate_plots.py combined_annotated.tsv \\
      --count-matrix count_matrix.tsv --metadata metadata.tsv -o plots/

  # With HMMER domain annotations on heatmap rows:
  python pydeseq2_generate_plots.py combined_annotated.tsv \\
      --count-matrix count_matrix.tsv --metadata metadata.tsv \\
      --hmmer pfam_domains.txt -o plots/

  # Combined two-species heatmap (DC + DG on one plot):
  python pydeseq2_generate_plots.py DC_annotated.tsv \\
      --count-matrix DC_counts.tsv --metadata DC_metadata.tsv \\
      --count-matrix2 DG_counts.tsv --metadata2 DG_metadata.tsv \\
      --species1 DC --species2 DG -o plots/

  # From raw DESeq2 results (MA + volcano only, no heatmaps):
  python pydeseq2_generate_plots.py pydeseq2_results.tsv -o plots/
"""

import os
import re
import sys
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.patches import Patch
except ImportError as e:
    print(f"ERROR: Missing required package: {e}")
    print("  pip install matplotlib seaborn")
    sys.exit(1)


# ============================================================================
# GENE FAMILY PATTERNS (same as extract_gene_families.py)
# ============================================================================

GENE_FAMILIES = {
    'CYP': {
        'full_name': 'Cytochrome P450',
        'blast_patterns': [
            r'cytochrome\s+P[- ]?450',
            r'\bCYP\d',
            r'\bP450\b',
        ],
        'blast_exclude': [
            r'NADPH.cytochrome\s+P450\s+reductase',
            r'cytochrome\s+P450\s+reductase',
            r'cytochrome\s+b5',
            r'cytochrome\s+c\b',
        ],
        'hmmer_pfam': ['PF00067'],
        'hmmer_names': ['p450'],
    },
    'OMT': {
        'full_name': 'O-Methyltransferase',
        'blast_patterns': [
            r'O-methyltransferase',
            r'\bOMT\b',
            r'\bCOMT\b',
            r'caffeic\s+acid.*methyltransferase',
            r'caffeoyl.CoA\s+O-methyltransferase',
            r'\bCCoAOMT\b',
            r'catechol\s+O-methyltransferase',
            r'flavonoid.*O-methyltransferase',
            r'isoflavone.*O-methyltransferase',
            r'myricetin.*O-methyltransferase',
            r'trans-resveratrol.*O-methyltransferase',
        ],
        'blast_exclude': [
            r'DNA\s+methyltransferase',
            r'histone.*methyltransferase',
            r'rRNA.*methyltransferase',
            r'tRNA.*methyltransferase',
            r'N-methyltransferase',
            r'protein.*methyltransferase',
        ],
        'hmmer_pfam': ['PF00891', 'PF01596', 'PF08100'],
        'hmmer_names': ['Methyltransf_2', 'Methyltransf_3', 'Dimerisation'],
    },
}

FAMILY_COLORS = {'CYP': '#e74c3c', 'OMT': '#3498db'}
COND_COLORS = {'R': '#d35400', 'L': '#27ae60',
               'root': '#d35400', 'leaf': '#27ae60'}

# CYP subfamily palette -- distinct colors per subfamily
CYP_SUBFAMILY_PALETTE = {
    'CYP71':  '#e41a1c',
    'CYP71D': '#e41a1c',
    'CYP72':  '#377eb8',
    'CYP76':  '#4daf4a',
    'CYP81':  '#984ea3',
    'CYP85':  '#ff7f00',
    'CYP86':  '#a65628',
    'CYP89':  '#f781bf',
    'CYP90':  '#999999',
    'CYP94':  '#66c2a5',
    'CYP97':  '#fc8d62',
    'CYP98':  '#8da0cb',
    'CYP706': '#e78ac3',
    'CYP707': '#a6d854',
    'CYP714': '#ffd92f',
    'CYP716': '#e5c494',
    'CYP719': '#b3b3b3',
    'CYP720': '#66c2a5',
    'CYP734': '#8dd3c7',
    'CYP735': '#ffffb3',
    'CYP749': '#bebada',
    'CYP78':  '#fb8072',
    'other':  '#cccccc',
}

OMT_SUBFAMILY_PALETTE = {
    'COMT':     '#1f78b4',
    'CCoAOMT':  '#33a02c',
    'FOMT':     '#e31a1c',
    'IOMT':     '#ff7f00',
    'other':    '#cccccc',
}


def _parse_cyp_subfamily(desc):
    """Extract CYP subfamily (e.g. CYP71D, CYP81, CYP719) from a BLAST description."""
    if not desc or not isinstance(desc, str):
        return 'CYP unclassified'
    m = re.search(r'\b(CYP\d{1,3}[A-Z]?)\d*\b', desc, re.IGNORECASE)
    if m:
        return m.group(1).upper()
    m = re.search(r'P[- ]?450\s+(\d{1,3})([A-Z])?', desc)
    if m:
        num = m.group(1)
        letter = m.group(2) or ''
        return f"CYP{num}{letter}"
    if re.search(r'cytochrome\s+P[- ]?450|monooxygenase', desc, re.IGNORECASE):
        return 'CYP unclassified'
    return 'CYP unclassified'


def _parse_omt_subfamily(desc):
    """Extract OMT subfamily from a BLAST description."""
    if not desc or not isinstance(desc, str):
        return 'other'
    d = desc.lower()
    if 'ccoaomt' in d or 'caffeoyl-coa' in d or 'caffeoyl.coa' in d:
        return 'CCoAOMT'
    if re.search(r'\bcomt\b', d) or 'caffeic acid' in d:
        return 'COMT'
    if 'flavonoid' in d:
        return 'FOMT'
    if 'isoflavone' in d:
        return 'IOMT'
    return 'other'


def _build_subfamily_map(gene_ids, results_df, family_name):
    """Build gene_id -> subfamily mapping from BLAST descriptions."""
    parse_fn = _parse_cyp_subfamily if family_name == 'CYP' else _parse_omt_subfamily
    subfam_map = {}
    id_col = 'gene_id' if 'gene_id' in results_df.columns else None
    for gid in gene_ids:
        if id_col:
            match = results_df[results_df[id_col] == gid]
        else:
            match = results_df.loc[[gid]] if gid in results_df.index else pd.DataFrame()
        desc = ''
        if not match.empty:
            desc = match.iloc[0].get('blast_description', '')
        subfam_map[gid] = parse_fn(desc if pd.notna(desc) else '')
    return subfam_map


# ============================================================================
# GENE FAMILY DETECTION
# ============================================================================

def detect_gene_families(results_df, domain_map=None):
    """Identify CYP and OMT genes using BLAST descriptions AND HMMER Pfam domains.

    Returns a dict mapping gene_id -> family name.  Genes found by either
    evidence source (or both) are included.
    """
    family_map = {}

    if 'gene_family' in results_df.columns:
        print("  Detecting gene families from 'gene_family' column...")
        for idx, row in results_df.iterrows():
            gid = row.get('gene_id', idx)
            fam = row['gene_family']
            if pd.notna(fam) and fam in GENE_FAMILIES:
                family_map[gid] = fam
        print(f"    Found {len(family_map)} annotated genes")
        return family_map

    blast_map = {}
    hmmer_map = {}

    # --- BLAST-based detection ---
    if 'blast_description' in results_df.columns:
        print("  Detecting gene families from BLAST descriptions...")
        desc_col = results_df['blast_description'].fillna('')

        for fam_name, fam_def in GENE_FAMILIES.items():
            include_mask = pd.Series(False, index=results_df.index)
            for pattern in fam_def['blast_patterns']:
                include_mask |= desc_col.str.contains(pattern, case=False, regex=True)

            exclude_mask = pd.Series(False, index=results_df.index)
            for pattern in fam_def['blast_exclude']:
                exclude_mask |= desc_col.str.contains(pattern, case=False, regex=True)

            hits = results_df[include_mask & ~exclude_mask]
            for idx in hits.index:
                gid = hits.at[idx, 'gene_id'] if 'gene_id' in hits.columns else idx
                if gid not in blast_map:
                    blast_map[gid] = fam_name

            n = (include_mask & ~exclude_mask).sum()
            print(f"    BLAST  {fam_name}: {n} genes matched")
    else:
        print("  No blast_description column -- skipping BLAST-based detection")

    # --- HMMER-based detection ---
    if domain_map:
        print("  Detecting gene families from HMMER Pfam domains...")
        for gid, domains_str in domain_map.items():
            domains_lower = domains_str.lower()
            for fam_name, fam_def in GENE_FAMILIES.items():
                matched = False
                for pfam_id in fam_def.get('hmmer_pfam', []):
                    if pfam_id.lower() in domains_lower:
                        matched = True
                        break
                if not matched:
                    for name in fam_def.get('hmmer_names', []):
                        if name.lower() in domains_lower:
                            matched = True
                            break
                if matched and gid not in hmmer_map:
                    hmmer_map[gid] = fam_name

        for fam_name in GENE_FAMILIES:
            n = sum(1 for f in hmmer_map.values() if f == fam_name)
            print(f"    HMMER  {fam_name}: {n} genes matched")
    else:
        print("  No HMMER data provided -- skipping domain-based detection")

    # --- Combine: union of BLAST and HMMER ---
    all_ids = set(blast_map.keys()) | set(hmmer_map.keys())
    for gid in all_ids:
        fam = blast_map.get(gid) or hmmer_map.get(gid)
        family_map[gid] = fam

    blast_only = set(blast_map.keys()) - set(hmmer_map.keys())
    hmmer_only = set(hmmer_map.keys()) - set(blast_map.keys())
    both = set(blast_map.keys()) & set(hmmer_map.keys())
    print(f"\n    Combined: {len(family_map)} total "
          f"(BLAST only: {len(blast_only)}, HMMER only: {len(hmmer_only)}, "
          f"both: {len(both)})")

    return family_map


def parse_hmmer_domains(hmmer_file):
    """Parse HMMER domtblout file, return dict: gene_id -> domain string."""
    if hmmer_file is None or not Path(hmmer_file).exists():
        return {}

    print(f"  Reading HMMER domains: {hmmer_file}")
    domain_map = {}
    with open(hmmer_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 23:
                continue
            gid = parts[3].split('|')[0]
            domain = f"{parts[0]}({parts[1]})"
            if gid not in domain_map:
                domain_map[gid] = set()
            domain_map[gid].add(domain)

    for gid in domain_map:
        domain_map[gid] = '; '.join(sorted(domain_map[gid]))

    print(f"    Loaded domains for {len(domain_map)} genes")
    return domain_map


# ============================================================================
# GTF BIOTYPE PARSING
# ============================================================================

BIOTYPE_PALETTE = {
    'protein_coding': '#e41a1c',
    'lncRNA':         '#377eb8',
    'pseudogene':     '#4daf4a',
    'miRNA':          '#984ea3',
    'snRNA':          '#ff7f00',
    'snoRNA':         '#ffff33',
    'rRNA':           '#a65628',
    'misc_RNA':       '#f781bf',
    'ncRNA':          '#377eb8',
}


def parse_biotype_from_gtf(gtf_path):
    """Extract gene_id -> gene_biotype mapping from a GTF file.

    Scans gene-level rows for the 'gene_biotype' or 'gene_type' attribute.
    Returns an empty dict if the file doesn't exist or contains no biotypes.
    """
    if gtf_path is None or not Path(gtf_path).exists():
        return {}

    print(f"  Parsing gene biotype from GTF: {gtf_path}")
    biotype_map = {}
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] != 'gene':
                continue
            attrs = fields[8]
            gid = bt = None
            for token in attrs.split(';'):
                token = token.strip()
                if not token:
                    continue
                kv = token.split(' ', 1)
                if len(kv) != 2:
                    continue
                k, v = kv[0], kv[1].strip('"')
                if k == 'gene_id':
                    gid = v
                elif k in ('gene_biotype', 'gene_type'):
                    bt = v
            if gid and bt:
                biotype_map[gid] = bt

    print(f"    Found biotype for {len(biotype_map)} genes "
          f"({len(set(biotype_map.values()))} categories)")
    return biotype_map


def _biotype_row_colors(gene_ids, biotype_map):
    """Build a Series of colors from biotype_map for the given gene_ids."""
    if not biotype_map:
        return None
    colors = []
    for gid in gene_ids:
        bt = biotype_map.get(gid, 'other')
        colors.append(BIOTYPE_PALETTE.get(bt, '#999999'))
    return colors


# ============================================================================
# SAMPLE NAMING
# ============================================================================

def _build_sample_display_names(samples, metadata, condition_col='condition',
                                species_label=None, label_style='short'):
    """Build display names from sample metadata.

    label_style:
        'short'  -> DC1L, DC2L, DC1R  (for PCA, correlation)
        'long'   -> Leaf 1, Leaf 2, Root 1  (for single-species heatmaps)
        'species_long' -> DC Leaf 1, DG Root 2  (for combined heatmaps)
    """
    meta = metadata.copy()
    sample_col = meta.columns[0]
    if sample_col != condition_col and meta.index.name is None:
        meta = meta.set_index(sample_col)

    cond_suffix = {'R': 'R', 'L': 'L', 'root': 'R', 'leaf': 'L'}
    cond_word = {'R': 'Root', 'L': 'Leaf', 'root': 'Root', 'leaf': 'Leaf'}
    counters = {}
    rename = {}

    if species_label is None:
        first = str(samples[0]) if samples else ''
        m = re.match(r'^([A-Z]{2})', first)
        species_label = m.group(1) if m else ''

    for s in samples:
        cond = meta.loc[s, condition_col] if s in meta.index else 'X'
        suffix = cond_suffix.get(cond, str(cond)[0].upper())
        counters[suffix] = counters.get(suffix, 0) + 1
        n = counters[suffix]

        if label_style == 'long':
            word = cond_word.get(cond, cond)
            rename[s] = f"{word} {n}"
        elif label_style == 'species_long':
            word = cond_word.get(cond, cond)
            rename[s] = f"{species_label} {word} {n}"
        else:
            rename[s] = f"{species_label}{n}{suffix}"

    return rename


# ============================================================================
# NORMALIZATION
# ============================================================================

def calculate_normalized_counts(count_matrix):
    """DESeq2-style median-of-ratios normalization."""
    print("\n  Calculating normalized counts...")

    geo_means = count_matrix.apply(
        lambda row: np.exp(np.log(row[row > 0]).mean()) if (row > 0).sum() > 0 else 0,
        axis=1
    )

    size_factors = []
    for col in count_matrix.columns:
        ratios = count_matrix[col] / geo_means
        ratios = ratios[np.isfinite(ratios) & (ratios > 0)]
        size_factors.append(ratios.median())

    size_factors = pd.Series(size_factors, index=count_matrix.columns)
    print("    Size factors: " + ", ".join(f"{s}={sf:.3f}" for s, sf in size_factors.items()))

    return count_matrix.div(size_factors, axis=1)


# ============================================================================
# MA + VOLCANO PLOTS
# ============================================================================

def generate_ma_plot(results_df, output_dir, alpha=0.05, lfc_cutoff=2.0,
                     contrast_A='R', contrast_B='L'):
    """MA plot following DESeq2 plotMA conventions with directional coloring.

    Standard elements (Bioconductor DESeq2 / DESeqAnalysis plotMa):
      - baseMean (x, log-scale) vs log2FoldChange (y)
      - Directional coloring: red (up), blue (down), gray (NS)
      - LFC threshold lines
      - Gene counts per direction annotated
      - Caption with thresholds and normalization method
    """
    print("  Creating MA plot...")
    fig, ax = plt.subplots(figsize=(9, 7))

    valid = (results_df['log2FoldChange'].notna() & results_df['baseMean'].notna()
             & np.isfinite(results_df['log2FoldChange']) & np.isfinite(results_df['baseMean']))
    sig = valid & results_df['padj'].notna() & (results_df['padj'] < alpha)
    sig_up = sig & (results_df['log2FoldChange'] > lfc_cutoff)
    sig_down = sig & (results_df['log2FoldChange'] < -lfc_cutoff)
    sig_mid = sig & ~sig_up & ~sig_down
    ns = valid & ~sig

    cond_up = {'R': 'Root', 'L': 'Leaf'}.get(contrast_A, contrast_A)
    cond_down = {'R': 'Root', 'L': 'Leaf'}.get(contrast_B, contrast_B)

    ax.scatter(results_df.loc[ns, 'baseMean'],
               results_df.loc[ns, 'log2FoldChange'],
               alpha=0.3, s=1.5, c='#AAAAAA', rasterized=True, label='NS')
    if sig_mid.sum() > 0:
        ax.scatter(results_df.loc[sig_mid, 'baseMean'],
                   results_df.loc[sig_mid, 'log2FoldChange'],
                   alpha=0.5, s=3, c='#636363', rasterized=True,
                   label=f'padj < {alpha} only')
    if sig_up.sum() > 0:
        ax.scatter(results_df.loc[sig_up, 'baseMean'],
                   results_df.loc[sig_up, 'log2FoldChange'],
                   alpha=0.6, s=4, c='#d62728', rasterized=True,
                   label=f'Up in {cond_up} ({sig_up.sum():,})')
    if sig_down.sum() > 0:
        ax.scatter(results_df.loc[sig_down, 'baseMean'],
                   results_df.loc[sig_down, 'log2FoldChange'],
                   alpha=0.6, s=4, c='#1f77b4', rasterized=True,
                   label=f'Up in {cond_down} ({sig_down.sum():,})')

    ax.axhline(y=0, color='black', linestyle='-', linewidth=0.8)
    ax.axhline(y=lfc_cutoff, color='gray', linestyle='--', linewidth=0.5, alpha=0.6)
    ax.axhline(y=-lfc_cutoff, color='gray', linestyle='--', linewidth=0.5, alpha=0.6)

    n_total = valid.sum()
    n_sig = sig.sum()
    ax.set_xlabel('Mean of normalized counts', fontsize=12)
    ax.set_ylabel('Log$_2$ fold change', fontsize=12)
    ax.set_title(f'MA Plot: {cond_up} vs {cond_down}',
                 fontsize=14, fontweight='bold')
    ax.text(0.5, 1.01,
            f'{n_sig:,} DE genes of {n_total:,} total',
            transform=ax.transAxes, ha='center', va='bottom',
            fontsize=10, color='#444444')
    ax.set_xscale('log')
    ax.legend(loc='upper left', fontsize=8, framealpha=0.9, markerscale=3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.annotate(f'{sig_up.sum():,} up in {cond_up}', xy=(0.98, 0.96),
                xycoords='axes fraction', ha='right', fontsize=9,
                color='#d62728', fontweight='bold')
    ax.annotate(f'{sig_down.sum():,} up in {cond_down}', xy=(0.98, 0.04),
                xycoords='axes fraction', ha='right', fontsize=9,
                color='#1f77b4', fontweight='bold')

    fig.text(0.5, 0.01,
             f'Thresholds: padj < {alpha}, |log$_2$FC| > {lfc_cutoff}  |  '
             f'Normalization: DESeq2 median-of-ratios',
             ha='center', fontsize=8, color='#666666', style='italic')

    plt.tight_layout(rect=[0, 0.03, 1, 1])
    plt.savefig(output_dir / "ma_plot.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / "ma_plot.png", dpi=150, bbox_inches='tight')
    plt.close()
    print(f"    Saved: {output_dir / 'ma_plot.pdf'}")


def generate_volcano_plot(results_df, output_dir, padj_cutoff=0.05,
                          lfc_cutoff=2.0, top_n_labels=10,
                          contrast_A='R', contrast_B='L'):
    """Enhanced Volcano plot modeled on EnhancedVolcano (Blighe et al. 2018).

    Standard elements:
      - 4 color categories: NS, |log2FC| only, padj only, both
      - Dashed threshold lines for padj and log2FC cutoffs
      - Top N significant genes labeled with connector lines
      - Per-category gene counts in legend
      - Subtitle with DE gene summary
      - Caption with exact threshold parameters
    """
    print("  Creating enhanced volcano plot...")
    fig, ax = plt.subplots(figsize=(10, 8.5))

    valid = (results_df['log2FoldChange'].notna() & results_df['padj'].notna()
             & np.isfinite(results_df['log2FoldChange']) & np.isfinite(results_df['padj']))
    df = results_df.loc[valid].copy()
    df['neg_log10_padj'] = -np.log10(df['padj'] + 1e-300)

    pass_padj = df['padj'] < padj_cutoff
    pass_lfc = df['log2FoldChange'].abs() > lfc_cutoff
    sig_up = pass_padj & pass_lfc & (df['log2FoldChange'] > 0)
    sig_down = pass_padj & pass_lfc & (df['log2FoldChange'] < 0)

    cat_ns = ~pass_padj & ~pass_lfc
    cat_lfc = pass_lfc & ~pass_padj
    cat_padj = pass_padj & ~pass_lfc
    cat_both = pass_padj & pass_lfc

    cond_up = {'R': 'Root', 'L': 'Leaf'}.get(contrast_A, contrast_A)
    cond_down = {'R': 'Root', 'L': 'Leaf'}.get(contrast_B, contrast_B)

    for mask, color, label in [
        (cat_ns,   '#AAAAAA', 'NS'),
        (cat_lfc,  '#2ca02c', f'|Log$_2$FC| > {lfc_cutoff}'),
        (cat_padj, '#1f77b4', f'padj < {padj_cutoff}'),
        (cat_both, '#d62728', f'padj < {padj_cutoff} & |Log$_2$FC| > {lfc_cutoff}'),
    ]:
        subset = df.loc[mask]
        if len(subset) == 0:
            continue
        ax.scatter(subset['log2FoldChange'], subset['neg_log10_padj'],
                   alpha=0.5, s=2, c=color, label=f'{label} ({mask.sum():,})',
                   rasterized=True)

    ax.axhline(y=-np.log10(padj_cutoff), color='black', linestyle='--',
               linewidth=0.6, alpha=0.5)
    ax.axvline(x=lfc_cutoff, color='black', linestyle='--',
               linewidth=0.6, alpha=0.5)
    ax.axvline(x=-lfc_cutoff, color='black', linestyle='--',
               linewidth=0.6, alpha=0.5)

    if top_n_labels > 0 and cat_both.sum() > 0:
        sig_df = df.loc[cat_both].nlargest(top_n_labels, 'neg_log10_padj')
        gene_id_col = 'gene_id' if 'gene_id' in sig_df.columns else None
        for idx, row in sig_df.iterrows():
            label_text = row[gene_id_col] if gene_id_col else str(idx)
            ax.annotate(
                label_text,
                xy=(row['log2FoldChange'], row['neg_log10_padj']),
                fontsize=6, fontweight='bold', alpha=0.85,
                xytext=(8, 6), textcoords='offset points',
                bbox=dict(boxstyle='round,pad=0.2', fc='white',
                          ec='gray', alpha=0.7, lw=0.5),
                arrowprops=dict(arrowstyle='->', color='gray',
                                lw=0.6, alpha=0.6,
                                connectionstyle='arc3,rad=0.15'))

    ax.set_xlabel('Log$_2$ fold change', fontsize=12)
    ax.set_ylabel('$-$Log$_{10}$ adjusted p-value', fontsize=12)
    ax.set_title(f'Differential Expression: {cond_up} vs {cond_down}',
                 fontsize=14, fontweight='bold')

    n_total = len(df)
    ax.text(0.5, 1.01,
            f'{sig_up.sum():,} up in {cond_up}  |  '
            f'{sig_down.sum():,} up in {cond_down}  |  '
            f'{n_total:,} total genes',
            transform=ax.transAxes, ha='center', va='bottom',
            fontsize=9.5, color='#444444')

    fig.text(0.5, 0.01,
             f'Log$_2$ fold change cutoff: {lfc_cutoff}  |  '
             f'Adjusted p-value cutoff: {padj_cutoff}  |  '
             f'Test: Wald (PyDESeq2)',
             ha='center', fontsize=8, color='#666666', style='italic')

    ax.legend(loc='upper right', fontsize=7.5, framealpha=0.9, markerscale=3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout(rect=[0, 0.03, 1, 1])
    plt.savefig(output_dir / "volcano_plot.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / "volcano_plot.png", dpi=150, bbox_inches='tight')
    plt.close()
    print(f"    Saved: {output_dir / 'volcano_plot.pdf'}")


def generate_ma_volcano_plots(results_df, output_dir,
                              contrast_A='R', contrast_B='L'):
    """Generate MA and Enhanced Volcano plots."""
    print("\nGenerating MA and Volcano plots...")
    generate_ma_plot(results_df, output_dir, contrast_A=contrast_A,
                     contrast_B=contrast_B)
    generate_volcano_plot(results_df, output_dir, contrast_A=contrast_A,
                          contrast_B=contrast_B)


# ============================================================================
# PCA PLOT
# ============================================================================

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


def generate_pca_plot(count_matrix, metadata, output_dir, n_top=500,
                      condition_col='condition', species_label=None):
    """PCA of top-N most variable genes with sample labels and 95% confidence ellipses."""
    print("\nGenerating PCA plot...")

    counts = count_matrix[~count_matrix.index.str.startswith('N_')]
    norm = calculate_normalized_counts(counts)
    log_norm = np.log2(norm + 1)

    variances = log_norm.var(axis=1)
    n_top_actual = min(n_top, len(variances))
    top_genes = variances.nlargest(n_top_actual).index
    subset = log_norm.loc[top_genes]

    X = subset.values.T
    X_centered = X - X.mean(axis=0)

    U, S, Vt = np.linalg.svd(X_centered, full_matrices=False)
    explained = (S ** 2) / (S ** 2).sum() * 100
    pc1 = U[:, 0] * S[0]
    pc2 = U[:, 1] * S[1]

    meta = metadata.copy()
    sample_col = meta.columns[0]
    if sample_col != condition_col:
        meta = meta.set_index(sample_col)

    display_names = _build_sample_display_names(
        list(subset.columns), metadata, condition_col, species_label)

    fig, ax = plt.subplots(figsize=(8, 6))
    cond_palette = {'R': '#d35400', 'L': '#27ae60',
                    'root': '#d35400', 'leaf': '#27ae60'}
    cond_labels = {'R': 'Root', 'L': 'Leaf',
                   'root': 'Root', 'leaf': 'Leaf'}

    group_points = {}
    for i, sample in enumerate(subset.columns):
        cond = meta.loc[sample, condition_col] if sample in meta.index else 'unknown'
        color = cond_palette.get(cond, '#bdc3c7')
        ax.scatter(pc1[i], pc2[i], c=color, s=100, edgecolors='white',
                   linewidths=0.8, zorder=3)
        short_name = display_names.get(sample, sample)
        ax.annotate(short_name, (pc1[i], pc2[i]), fontsize=8, fontweight='bold',
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
    plt.savefig(output_dir / "pca_plot.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / "pca_plot.png", dpi=150, bbox_inches='tight')
    plt.close()
    print(f"    Saved: {output_dir / 'pca_plot.pdf'}")


# ============================================================================
# SAMPLE CORRELATION HEATMAP
# ============================================================================

def generate_sample_correlation_heatmap(count_matrix, metadata, output_dir,
                                        condition_col='condition',
                                        species_label=None):
    """Sample-to-sample distance heatmap with blue gradient and short sample names."""
    print("\nGenerating sample correlation heatmap...")

    counts = count_matrix[~count_matrix.index.str.startswith('N_')]
    norm = calculate_normalized_counts(counts)
    log_norm = np.log2(norm + 1)

    display_names = _build_sample_display_names(
        list(log_norm.columns), metadata, condition_col, species_label)

    from scipy.spatial.distance import pdist, squareform
    dist_vec = pdist(log_norm.T.values, metric='euclidean')
    dist_mat = squareform(dist_vec)
    short_labels = [display_names.get(s, s) for s in log_norm.columns]
    dist_df = pd.DataFrame(dist_mat, index=short_labels, columns=short_labels)

    meta = metadata.copy()
    sample_col = meta.columns[0]
    if sample_col != condition_col:
        meta = meta.set_index(sample_col)

    cond_palette = {'R': '#d35400', 'L': '#27ae60',
                    'root': '#d35400', 'leaf': '#27ae60'}
    col_colors = []
    for sample in log_norm.columns:
        cond = meta.loc[sample, condition_col] if sample in meta.index else ''
        col_colors.append(cond_palette.get(cond, '#bdc3c7'))
    col_colors_series = pd.Series(col_colors, index=short_labels)

    n = len(dist_df)
    fig_size = max(6, 0.6 * n + 2)

    g = sns.clustermap(
        dist_df,
        cmap='Blues_r',
        figsize=(fig_size, fig_size),
        linewidths=0.5,
        linecolor='white',
        col_colors=col_colors_series,
        row_colors=col_colors_series,
        xticklabels=True,
        yticklabels=True,
        cbar_kws={'label': 'Euclidean distance', 'shrink': 0.5},
        dendrogram_ratio=0.12,
    )

    g.ax_heatmap.set_xlabel('')
    g.ax_heatmap.set_ylabel('')
    g.fig.suptitle('Sample-to-Sample Distance', fontsize=13,
                   fontweight='bold', y=1.02)

    legend_elements = [
        Patch(facecolor='#d35400', label='Root'),
        Patch(facecolor='#27ae60', label='Leaf'),
    ]
    g.ax_heatmap.legend(handles=legend_elements, loc='upper left',
                        bbox_to_anchor=(1.12, 1.0), frameon=True, fontsize=9)

    plt.savefig(output_dir / "sample_correlation_heatmap.pdf",
                dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / "sample_correlation_heatmap.png",
                dpi=150, bbox_inches='tight')
    plt.close()
    print(f"    Saved: {output_dir / 'sample_correlation_heatmap.pdf'}")


# ============================================================================
# GENE FAMILY HEATMAP
# ============================================================================

def _draw_subfamily_brackets(g, subfam_list, n_genes):
    """Draw bracket annotations for subfamily groups to the left of the heatmap.

    Genes must be pre-sorted so that same-subfamily genes are contiguous.
    Uses the row dendrogram axes area for bracket lines and rotated labels,
    matching the reference figure style.
    """
    blocks = []
    current_sf = None
    block_start = None
    for row_i, sf in enumerate(subfam_list):
        if sf != current_sf:
            if current_sf is not None and current_sf != 'other':
                blocks.append((current_sf, block_start, row_i - 1))
            current_sf = sf
            block_start = row_i
    if current_sf is not None and current_sf != 'other':
        blocks.append((current_sf, block_start, len(subfam_list) - 1))

    if not blocks:
        return

    ax = g.ax_row_dendrogram
    ax.clear()
    hm_ylim = g.ax_heatmap.get_ylim()
    ax.set_ylim(hm_ylim)
    ax.set_xlim(0, 1)
    ax.axis('off')

    label_size = max(8, min(11, 160 // max(len(blocks), 1)))

    for sf_name, start, end in blocks:
        y_top = start + 0.05
        y_bot = end + 0.95
        y_mid = (y_top + y_bot) / 2

        bx = 0.82
        ax.plot([bx, bx], [y_top, y_bot],
                color='black', linewidth=1.5, clip_on=False)
        ax.plot([bx, 1.0], [y_top, y_top],
                color='black', linewidth=1.2, clip_on=False)
        ax.plot([bx, 1.0], [y_bot, y_bot],
                color='black', linewidth=1.2, clip_on=False)

        ax.text(0.35, y_mid, sf_name, ha='center', va='center',
                fontsize=label_size, fontweight='bold', rotation=90,
                clip_on=False)


def generate_family_heatmap(results_df, gene_ids, family_name, full_name,
                            count_matrix, metadata, domain_map, output_dir,
                            scale='center', cluster_rows=True,
                            condition_col='condition', biotype_map=None,
                            top_n=None, species_label=None):
    """Generate a heatmap matching the reference figure style:
      - Gene locus IDs as row labels on the left
      - Genes grouped by subfamily with bracket annotations on the left
      - Short sample names (DC1L, DC2L, DC1R, etc.), Leaf first then Root
      - Centered log2 RdBu_r color scale (red = high, blue = low)
    """
    print(f"\n{'='*60}")
    print(f"Generating {family_name} Heatmap ({full_name})")
    print(f"{'='*60}")

    if top_n and len(gene_ids) > top_n:
        id_col = 'gene_id' if 'gene_id' in results_df.columns else None
        if id_col:
            sub = results_df[results_df[id_col].isin(gene_ids)].copy()
        else:
            sub = results_df.loc[results_df.index.intersection(gene_ids)].copy()
        if 'padj' in sub.columns:
            sub['_abs_lfc'] = sub['log2FoldChange'].abs() if 'log2FoldChange' in sub.columns else 0
            sub = sub.sort_values(['padj', '_abs_lfc'], ascending=[True, False])
            gene_ids = (sub[id_col] if id_col else sub.index).tolist()[:top_n]
            sub.drop(columns=['_abs_lfc'], inplace=True)
        else:
            gene_ids = gene_ids[:top_n]
        print(f"  Limiting to top {top_n} genes (by padj, then |log2FC|)")

    counts = count_matrix[~count_matrix.index.str.startswith('N_')]

    present = [g for g in gene_ids if g in counts.index]
    missing = [g for g in gene_ids if g not in counts.index]
    if missing:
        print(f"  WARNING: {len(missing)} {family_name} genes not in count matrix")
    if not present:
        print(f"  No {family_name} genes in count matrix -- skipping heatmap")
        return

    print(f"  Normalizing counts for {len(present)} {family_name} genes...")
    norm_counts = calculate_normalized_counts(counts)
    heatmap_data = np.log2(norm_counts.loc[present] + 1)

    if scale == 'center':
        row_means = heatmap_data.mean(axis=1)
        heatmap_data = heatmap_data.subtract(row_means, axis=0)
        cbar_label = "Centered log$_2$(norm + 1)"
    elif scale == 'zscore':
        row_means = heatmap_data.mean(axis=1)
        row_stds = heatmap_data.std(axis=1).replace(0, 1)
        heatmap_data = heatmap_data.subtract(row_means, axis=0).div(row_stds, axis=0)
        cbar_label = "Z-score"
    else:
        cbar_label = "log$_2$(norm + 1)"

    meta = metadata.copy()
    if condition_col in meta.columns:
        sample_col = meta.columns[0] if meta.index.name is None else None
        if sample_col and sample_col != condition_col:
            meta = meta.set_index(sample_col)
        ordered = [s for s in meta.sort_values(condition_col).index if s in heatmap_data.columns]
        heatmap_data = heatmap_data[ordered]

    display_names = _build_sample_display_names(
        list(heatmap_data.columns), metadata, condition_col, species_label,
        label_style='long')
    heatmap_data.columns = [display_names.get(s, s) for s in heatmap_data.columns]

    col_colors = None
    if condition_col in meta.columns:
        col_colors_list = []
        for orig_sample in ordered:
            cond = meta.loc[orig_sample, condition_col] if orig_sample in meta.index else None
            col_colors_list.append(COND_COLORS.get(cond, '#bdc3c7'))
        col_colors = pd.Series(col_colors_list, index=heatmap_data.columns)

    subfam_map = _build_subfamily_map(present, results_df, family_name)

    gene_means = heatmap_data.mean(axis=1)
    gene_order = sorted(present,
                        key=lambda g: (subfam_map.get(g, 'other') == 'other',
                                       subfam_map.get(g, 'other'),
                                       gene_means.get(g, 0)))
    heatmap_data = heatmap_data.loc[gene_order]
    subfam_list = [subfam_map.get(g, 'other') for g in gene_order]

    n_genes = len(heatmap_data)
    n_samples = len(heatmap_data.columns)
    label_size = max(5, min(8, 250 // max(n_genes, 1)))
    fig_height = max(8, min(40, 0.18 * n_genes + 3))
    fig_width = max(8, 1.0 * n_samples + 3)

    print(f"  Plotting {n_genes} genes x {n_samples} samples...")
    has_gene_id = 'gene_id' in results_df.columns

    prefix = family_name.lower()
    try:
        g = sns.clustermap(
            heatmap_data,
            cmap='RdBu_r',
            center=0,
            figsize=(fig_width, fig_height),
            row_cluster=False,
            col_cluster=False,
            col_colors=col_colors,
            linewidths=0.3,
            linecolor='white',
            cbar_kws={'label': cbar_label, 'shrink': 0.5},
            yticklabels=True,
            xticklabels=True,
            dendrogram_ratio=(0.22, 0.02),
        )

        g.ax_heatmap.tick_params(axis='y', labelsize=label_size, length=0, pad=2)
        g.ax_heatmap.set_ylabel('')
        g.ax_heatmap.set_xlabel('')

        _draw_subfamily_brackets(g, subfam_list, n_genes)

        g.fig.suptitle(f'{full_name} ({family_name}) Expression  |  {n_genes} genes',
                       fontsize=13, fontweight='bold', y=1.02)

        legend_elements = [
            Patch(facecolor='#27ae60', label='Leaf'),
            Patch(facecolor='#d35400', label='Root'),
        ]
        g.ax_heatmap.legend(handles=legend_elements, loc='upper left',
                            bbox_to_anchor=(1.05, 1.0), frameon=True, fontsize=8,
                            title='Tissue', title_fontsize=9)

        pdf_path = output_dir / f"{prefix}_heatmap.pdf"
        g.savefig(pdf_path, dpi=300, bbox_inches='tight')
        g.savefig(output_dir / f"{prefix}_heatmap.png", dpi=150, bbox_inches='tight')
        plt.close()
        print(f"    Saved: {pdf_path}")

        matrix_path = output_dir / f"{prefix}_heatmap_matrix.tsv"
        heatmap_data.to_csv(matrix_path, sep='\t')
        print(f"    Saved: {matrix_path}")

    except Exception as e:
        print(f"  ERROR generating {family_name} heatmap: {e}")
        import traceback
        traceback.print_exc()

    gene_list_path = output_dir / f"{prefix}_gene_list.tsv"
    if has_gene_id:
        family_df = results_df[results_df['gene_id'].isin(gene_ids)]
    else:
        family_df = results_df.loc[results_df.index.intersection(gene_ids)]
    family_df.to_csv(gene_list_path, sep='\t', index=not has_gene_id)
    print(f"    Saved gene list: {gene_list_path} ({len(family_df)} genes)")


# ============================================================================
# COMBINED TWO-SPECIES HEATMAP
# ============================================================================

def _order_samples(meta, available_cols, condition_col):
    """Return sample names ordered leaf-first then root (matching reference figure)."""
    sample_col = meta.columns[0]
    root_conds = {'R', 'root'}
    leaf_conds = {'L', 'leaf'}
    root_samples, leaf_samples, other_samples = [], [], []
    for _, row in meta.iterrows():
        s = row[sample_col]
        if s not in available_cols:
            continue
        cond = row.get(condition_col, '')
        if cond in root_conds:
            root_samples.append(s)
        elif cond in leaf_conds:
            leaf_samples.append(s)
        else:
            other_samples.append(s)
    return leaf_samples + root_samples + other_samples


def _normalize_and_log(count_path, gene_ids):
    """Load count matrix, normalize, log2-transform, subset to gene_ids."""
    counts = pd.read_csv(count_path, sep='\t', index_col=0)
    counts = counts[~counts.index.str.startswith('N_')]
    present = [g for g in gene_ids if g in counts.index]
    norm = calculate_normalized_counts(counts)
    return np.log2(norm.loc[present] + 1), present


def generate_combined_family_heatmap(
        results_df, gene_ids, family_name, full_name,
        count_matrix_path1, metadata_path1,
        count_matrix_path2, metadata_path2,
        species1, species2, domain_map, output_dir,
        scale='center', cluster_rows=True, condition_col='condition',
        biotype_map=None, top_n=None):
    """Generate a single heatmap with two species side by side.

    Columns: SP1-Leaf | SP1-Root | SP2-Leaf | SP2-Root.
    Gene locus IDs as row labels, subfamily bracket annotations on left.
    Two stacked column color bars for species and tissue type.
    """

    print(f"\n{'='*60}")
    print(f"Generating Combined {family_name} Heatmap: {species1} + {species2}")
    print(f"{'='*60}")

    if not gene_ids:
        print(f"  No {family_name} genes -- skipping")
        return

    if top_n and len(gene_ids) > top_n:
        id_col = 'gene_id' if 'gene_id' in results_df.columns else None
        if id_col:
            sub = results_df[results_df[id_col].isin(gene_ids)].copy()
        else:
            sub = results_df.loc[results_df.index.intersection(gene_ids)].copy()
        if 'padj' in sub.columns:
            sub['_abs_lfc'] = sub['log2FoldChange'].abs() if 'log2FoldChange' in sub.columns else 0
            sub = sub.sort_values(['padj', '_abs_lfc'], ascending=[True, False])
            gene_ids = (sub[id_col] if id_col else sub.index).tolist()[:top_n]
            sub.drop(columns=['_abs_lfc'], inplace=True)
        else:
            gene_ids = gene_ids[:top_n]
        print(f"  Limiting to top {top_n} genes (by padj, then |log2FC|)")

    print(f"  Loading and normalizing {species1} counts...")
    log1, present1 = _normalize_and_log(count_matrix_path1, gene_ids)
    print(f"    {species1}: {len(present1)} genes found")

    print(f"  Loading and normalizing {species2} counts...")
    log2, present2 = _normalize_and_log(count_matrix_path2, gene_ids)
    print(f"    {species2}: {len(present2)} genes found")

    shared = [g for g in gene_ids if g in log1.index and g in log2.index]
    if not shared:
        print(f"  No {family_name} genes found in both count matrices -- skipping")
        return
    print(f"  Shared genes: {len(shared)}")

    log1 = log1.loc[shared]
    log2 = log2.loc[shared]

    meta1 = pd.read_csv(metadata_path1, sep='\t')
    meta2 = pd.read_csv(metadata_path2, sep='\t')

    sp1_ordered = _order_samples(meta1, set(log1.columns), condition_col)
    sp2_ordered = _order_samples(meta2, set(log2.columns), condition_col)

    heatmap_data = pd.concat([log1[sp1_ordered], log2[sp2_ordered]], axis=1)

    if scale == 'center':
        row_means = heatmap_data.mean(axis=1)
        heatmap_data = heatmap_data.subtract(row_means, axis=0)
        cbar_label = "Centered log$_2$(norm + 1)"
    elif scale == 'zscore':
        row_means = heatmap_data.mean(axis=1)
        row_stds = heatmap_data.std(axis=1).replace(0, 1)
        heatmap_data = heatmap_data.subtract(row_means, axis=0).div(row_stds, axis=0)
        cbar_label = "Z-score"
    else:
        cbar_label = "log$_2$(norm + 1)"

    sp1_display = _build_sample_display_names(
        sp1_ordered, meta1, condition_col, species1, label_style='species_long')
    sp2_display = _build_sample_display_names(
        sp2_ordered, meta2, condition_col, species2, label_style='species_long')
    all_display = {**sp1_display, **sp2_display}
    orig_columns = list(heatmap_data.columns)
    heatmap_data.columns = [all_display.get(s, s) for s in heatmap_data.columns]

    sp1_set = set(sp1_ordered)
    species_palette = {species1: '#8e44ad', species2: '#2980b9'}
    tissue_palette = {'Root': '#d35400', 'Leaf': '#27ae60'}
    root_conds = {'R', 'root'}
    leaf_conds = {'L', 'leaf'}

    sp1_sample_col = meta1.columns[0]
    sp2_sample_col = meta2.columns[0]
    meta1_idx = meta1.set_index(sp1_sample_col)
    meta2_idx = meta2.set_index(sp2_sample_col)

    species_colors = []
    tissue_colors = []
    for sample in orig_columns:
        if sample in sp1_set:
            species_colors.append(species_palette[species1])
            cond = meta1_idx.loc[sample, condition_col] if sample in meta1_idx.index else ''
        else:
            species_colors.append(species_palette[species2])
            cond = meta2_idx.loc[sample, condition_col] if sample in meta2_idx.index else ''

        if cond in root_conds:
            tissue_colors.append(tissue_palette['Root'])
        elif cond in leaf_conds:
            tissue_colors.append(tissue_palette['Leaf'])
        else:
            tissue_colors.append('#bdc3c7')

    col_colors_df = pd.DataFrame({
        'Species': species_colors,
        'Tissue': tissue_colors,
    }, index=heatmap_data.columns)

    subfam_map = _build_subfamily_map(shared, results_df, family_name)

    gene_means = heatmap_data.mean(axis=1)
    gene_order = sorted(shared,
                        key=lambda g: (subfam_map.get(g, 'other') == 'other',
                                       subfam_map.get(g, 'other'),
                                       gene_means.get(g, 0)))
    heatmap_data = heatmap_data.loc[gene_order]
    subfam_list = [subfam_map.get(g, 'other') for g in gene_order]

    n_genes = len(heatmap_data)
    n_samples = len(heatmap_data.columns)
    label_size = max(5, min(8, 250 // max(n_genes, 1)))
    fig_height = max(8, min(40, 0.18 * n_genes + 3))
    fig_width = max(10, 1.0 * n_samples + 3)

    print(f"  Plotting {n_genes} genes x {n_samples} samples...")
    prefix = family_name.lower()

    try:
        g = sns.clustermap(
            heatmap_data,
            cmap='RdBu_r',
            center=0,
            figsize=(fig_width, fig_height),
            row_cluster=False,
            col_cluster=False,
            col_colors=col_colors_df,
            linewidths=0.3,
            linecolor='white',
            cbar_kws={'label': cbar_label, 'shrink': 0.5},
            yticklabels=True,
            xticklabels=True,
            dendrogram_ratio=(0.22, 0.02),
        )

        g.ax_heatmap.tick_params(axis='y', labelsize=label_size, length=0, pad=2)
        g.ax_heatmap.set_ylabel('')
        g.ax_heatmap.set_xlabel('')

        _draw_subfamily_brackets(g, subfam_list, n_genes)

        g.fig.suptitle(
            f'{full_name} ({family_name}): {species1} vs {species2}  |  {n_genes} genes',
            fontsize=13, fontweight='bold', y=1.02,
        )

        legend_elements = [
            Patch(facecolor=species_palette[species1], label=species1),
            Patch(facecolor=species_palette[species2], label=species2),
            Patch(facecolor='none', edgecolor='none', label=''),
            Patch(facecolor='#27ae60', label='Leaf'),
            Patch(facecolor='#d35400', label='Root'),
        ]
        g.ax_heatmap.legend(handles=legend_elements, loc='upper left',
                            bbox_to_anchor=(1.05, 1.0), frameon=True, fontsize=8,
                            title='Legend', title_fontsize=9)

        pdf_path = output_dir / f"{prefix}_heatmap_combined.pdf"
        g.savefig(pdf_path, dpi=300, bbox_inches='tight')
        g.savefig(output_dir / f"{prefix}_heatmap_combined.png", dpi=150, bbox_inches='tight')
        plt.close()
        print(f"    Saved: {pdf_path}")

        matrix_path = output_dir / f"{prefix}_heatmap_combined_matrix.tsv"
        heatmap_data.to_csv(matrix_path, sep='\t')
        print(f"    Saved: {matrix_path}")

    except Exception as e:
        print(f"  ERROR generating combined {family_name} heatmap: {e}")
        import traceback
        traceback.print_exc()


# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="PyDESeq2 Step 3: Generate plots (MA, volcano, CYP/OMT heatmaps)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # From combined annotated file (auto-detects CYP/OMT):
  python pydeseq2_generate_plots.py combined_annotated.tsv \\
      --count-matrix count_matrix.tsv --metadata metadata.tsv -o plots/

  # With HMMER domain labels:
  python pydeseq2_generate_plots.py combined_annotated.tsv \\
      --count-matrix count_matrix.tsv --metadata metadata.tsv \\
      --hmmer pfam_domains.txt -o plots/

  # Raw DESeq2 results (MA + volcano only):
  python pydeseq2_generate_plots.py pydeseq2_results.tsv -o plots/
        """
    )

    parser.add_argument("results_file",
                        help="DESeq2 results TSV (unfiltered, filtered, combined annotated, or gene-family verified)")

    parser.add_argument("-o", "--output", default="plots",
                        help="Output directory (default: plots)")

    parser.add_argument("--count-matrix", default=None,
                        help="Count matrix TSV (required for heatmaps)")
    parser.add_argument("--metadata", default=None,
                        help="Metadata TSV (required for heatmaps)")
    parser.add_argument("--hmmer", default=None,
                        help="HMMER domtblout file (optional, domain info saved to gene list TSVs)")

    parser.add_argument("--contrast-factor", default="condition",
                        help="Metadata column for contrast (default: condition)")
    parser.add_argument("--contrast-A", default="R",
                        help="Numerator condition code (default: R)")
    parser.add_argument("--contrast-B", default="L",
                        help="Denominator condition code (default: L)")

    parser.add_argument("--scale", choices=['center', 'zscore'], default='center',
                        help="Heatmap scaling method (default: center)")
    parser.add_argument("--no-row-cluster", action="store_true",
                        help="Disable row clustering in heatmaps")
    parser.add_argument("--gtf", default=None,
                        help="GTF annotation file (optional, biotype info saved to gene list TSVs)")
    parser.add_argument("--top-n", type=int, default=None,
                        help="Only plot the top N genes per family, ranked by padj then |log2FC| (default: all)")
    parser.add_argument("--padj-cutoff", type=float, default=0.05,
                        help="Only include genes with padj < this value in heatmaps (default: 0.05)")
    parser.add_argument("--lfc-cutoff", type=float, default=2.0,
                        help="Only include genes with |log2FC| > this value in heatmaps (default: 2.0)")

    # Combined two-species heatmap mode
    parser.add_argument("--count-matrix2", default=None,
                        help="Second species count matrix (enables combined heatmap)")
    parser.add_argument("--metadata2", default=None,
                        help="Second species metadata TSV")
    parser.add_argument("--species1", default="SP1",
                        help="Label for species 1 (default: SP1)")
    parser.add_argument("--species2", default="SP2",
                        help="Label for species 2 (default: SP2)")

    args = parser.parse_args()

    if not os.path.exists(args.results_file):
        print(f"ERROR: Results file not found: {args.results_file}")
        return 1

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("PyDESeq2 Plot Generation (Step 3)")
    print("=" * 60)
    print(f"Input:  {args.results_file}")
    print(f"Output: {output_dir}")
    print()

    try:
        results_df = pd.read_csv(args.results_file, sep='\t', comment='#', index_col=0)
        print(f"Loaded {len(results_df)} genes")

        # Detect input type
        has_blast = 'blast_description' in results_df.columns
        has_family = 'gene_family' in results_df.columns
        print(f"  Input has blast_description: {has_blast}")
        print(f"  Input has gene_family:       {has_family}")

        sp_label = args.species1 if args.species1 != 'SP1' else None

        # --- MA + Volcano (always) ---
        if 'log2FoldChange' in results_df.columns and 'baseMean' in results_df.columns:
            generate_ma_volcano_plots(results_df, output_dir,
                                      contrast_A=args.contrast_A,
                                      contrast_B=args.contrast_B)
        else:
            print("\n  Skipping MA/Volcano (no log2FoldChange or baseMean columns)")

        # --- Plots requiring count matrix + metadata ---
        if not args.count_matrix or not args.metadata:
            print("\n  Count matrix and/or metadata not provided -- skipping heatmaps, PCA, correlation.")
            print("  For heatmaps, provide --count-matrix and --metadata")
        elif not os.path.exists(args.count_matrix):
            print(f"\n  WARNING: Count matrix not found: {args.count_matrix}")
        elif not os.path.exists(args.metadata):
            print(f"\n  WARNING: Metadata not found: {args.metadata}")
        else:
            count_matrix = pd.read_csv(args.count_matrix, sep='\t', index_col=0)
            metadata = pd.read_csv(args.metadata, sep='\t')

            # --- PCA plot ---
            try:
                generate_pca_plot(count_matrix, metadata, output_dir,
                                  condition_col=args.contrast_factor,
                                  species_label=sp_label)
            except Exception as e:
                print(f"  WARNING: PCA plot failed: {e}")

            # --- Sample correlation heatmap ---
            try:
                generate_sample_correlation_heatmap(
                    count_matrix, metadata, output_dir,
                    condition_col=args.contrast_factor,
                    species_label=sp_label)
            except Exception as e:
                print(f"  WARNING: Sample correlation heatmap failed: {e}")

            # --- Gene family heatmaps ---
            domain_map = parse_hmmer_domains(args.hmmer)
            biotype_map = parse_biotype_from_gtf(args.gtf)
            family_map = detect_gene_families(results_df, domain_map)

            if family_map and ('padj' in results_df.columns
                               and 'log2FoldChange' in results_df.columns):
                id_col = 'gene_id' if 'gene_id' in results_df.columns else None
                before_n = len(family_map)
                filtered_map = {}
                for gid, fam in family_map.items():
                    if id_col:
                        row = results_df[results_df[id_col] == gid]
                    else:
                        row = results_df.loc[[gid]] if gid in results_df.index else pd.DataFrame()
                    if row.empty:
                        continue
                    p = row.iloc[0].get('padj', np.nan)
                    lfc = row.iloc[0].get('log2FoldChange', 0)
                    if pd.notna(p) and p < args.padj_cutoff and abs(lfc) > args.lfc_cutoff:
                        filtered_map[gid] = fam
                family_map = filtered_map
                print(f"\n  DE filter: padj < {args.padj_cutoff}, |log2FC| > {args.lfc_cutoff}")
                print(f"    {before_n} pattern-matched -> {len(family_map)} differentially expressed")

            if not family_map:
                print("\n  No CYP or OMT genes detected (or none pass DE filters) -- skipping family heatmaps.")
                print("  (Input may be raw DESeq2 results without BLAST annotation)")
            else:
                for fam_name, fam_def in GENE_FAMILIES.items():
                    gene_ids = [gid for gid, fam in family_map.items() if fam == fam_name]
                    if not gene_ids:
                        print(f"\n  No {fam_name} genes found -- skipping {fam_name} heatmap")
                        continue
                    generate_family_heatmap(
                        results_df, gene_ids, fam_name, fam_def['full_name'],
                        count_matrix, metadata, domain_map, output_dir,
                        scale=args.scale,
                        cluster_rows=not args.no_row_cluster,
                        condition_col=args.contrast_factor,
                        biotype_map=biotype_map,
                        top_n=args.top_n,
                        species_label=sp_label,
                    )

                # Combined two-species heatmaps
                if (args.count_matrix2 and args.metadata2
                        and os.path.exists(args.count_matrix2)
                        and os.path.exists(args.metadata2)):
                    print(f"\n  Generating combined heatmaps ({args.species1} + {args.species2})...")
                    for fam_name, fam_def in GENE_FAMILIES.items():
                        fam_gene_ids = [gid for gid, fam in family_map.items() if fam == fam_name]
                        if not fam_gene_ids:
                            continue
                        generate_combined_family_heatmap(
                            results_df, fam_gene_ids, fam_name, fam_def['full_name'],
                            args.count_matrix, args.metadata,
                            args.count_matrix2, args.metadata2,
                            args.species1, args.species2,
                            domain_map, output_dir,
                            scale=args.scale,
                            cluster_rows=not args.no_row_cluster,
                            condition_col=args.contrast_factor,
                            biotype_map=biotype_map,
                            top_n=args.top_n,
                        )

        print("\n" + "=" * 60)
        print("Plot generation complete!")
        print("=" * 60)
        print(f"\nOutput: {output_dir}")
        return 0

    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
