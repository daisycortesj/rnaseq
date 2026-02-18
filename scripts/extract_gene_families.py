#!/usr/bin/env python3
"""
Extract and verify CYP and OMT gene families using multi-evidence approach

Searches BLAST annotations and HMMER Pfam domains to identify CYP (cytochrome
P450) and OMT (O-methyltransferase) genes. A gene is included if ANY source
flags it. Confidence is scored by how many sources agree.

HMMER is optional -- the script works with BLAST only and can be re-run
when HMMER results become available.

Usage:
  # BLAST only:
  python scripts/extract_gene_families.py \
      --blast 06_analysis/combined_DC/DC_swissprot_discovery_annotated.tsv \
      --species DC \
      --output 06_analysis/gene_families_DC

  # BLAST + HMMER:
  python scripts/extract_gene_families.py \
      --blast 06_analysis/combined_DC/DC_swissprot_discovery_annotated.tsv \
      --hmmer 06_analysis/hmmer_DC/all_genes_protein_pfam_domains.txt \
      --species DC \
      --output 06_analysis/gene_families_DC
"""

import pandas as pd
import numpy as np
import argparse
import re
import sys
from pathlib import Path
from datetime import datetime


# ============================================================================
# GENE FAMILY DEFINITIONS
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


# ============================================================================
# SEARCH FUNCTIONS
# ============================================================================

def search_blast(df, family_name, family_def):
    """Search BLAST annotations for a gene family."""
    if 'blast_description' not in df.columns:
        print(f"  WARNING: No blast_description column found")
        return set()

    desc = df['blast_description'].fillna('')

    include_mask = pd.Series(False, index=df.index)
    for pattern in family_def['blast_patterns']:
        include_mask |= desc.str.contains(pattern, case=False, regex=True)

    exclude_mask = pd.Series(False, index=df.index)
    for pattern in family_def['blast_exclude']:
        exclude_mask |= desc.str.contains(pattern, case=False, regex=True)

    hits = df[include_mask & ~exclude_mask]
    gene_ids = set(hits['gene_id'].values) if 'gene_id' in hits.columns else set(hits.index)

    print(f"  BLAST: {len(gene_ids)} {family_name} genes "
          f"(matched {include_mask.sum()}, excluded {exclude_mask.sum()})")

    return gene_ids


def parse_hmmer_domtblout(filepath):
    """Parse HMMER --domtblout format into a DataFrame."""
    rows = []
    with open(filepath) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 23:
                continue
            rows.append({
                'target_name': parts[0],
                'target_acc': parts[1],
                'query_name': parts[3],
                'query_acc': parts[4],
                'evalue': float(parts[6]),
                'score': float(parts[7]),
                'dom_evalue': float(parts[12]),
                'dom_score': float(parts[13]),
                'description': ' '.join(parts[22:]),
            })

    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows)


def search_hmmer(hmmer_file, family_name, family_def):
    """Search HMMER results for family-specific Pfam domains."""
    if hmmer_file is None or not Path(hmmer_file).exists():
        print(f"  HMMER: skipped (file not provided or not found)")
        return set(), {}

    hmmer_df = parse_hmmer_domtblout(hmmer_file)
    if hmmer_df.empty:
        print(f"  HMMER: skipped (empty results)")
        return set(), {}

    mask = pd.Series(False, index=hmmer_df.index)

    for pfam_id in family_def['hmmer_pfam']:
        mask |= hmmer_df['target_acc'].str.contains(pfam_id, case=False, na=False)

    for name in family_def['hmmer_names']:
        mask |= hmmer_df['target_name'].str.contains(name, case=False, na=False)

    hits = hmmer_df[mask]

    gene_ids = set()
    domain_map = {}
    for _, row in hits.iterrows():
        gid = row['query_name'].split('|')[0]
        gene_ids.add(gid)
        if gid not in domain_map:
            domain_map[gid] = []
        domain_map[gid].append(f"{row['target_name']}({row['target_acc']})")

    for gid in domain_map:
        domain_map[gid] = '; '.join(sorted(set(domain_map[gid])))

    print(f"  HMMER: {len(gene_ids)} {family_name} genes "
          f"({hits.shape[0]} domain hits)")

    return gene_ids, domain_map


# ============================================================================
# MAIN EXTRACTION
# ============================================================================

def extract_family(blast_df, hmmer_file, family_name, family_def):
    """Extract and verify a single gene family from all evidence sources."""
    print(f"\n{'='*60}")
    print(f"Extracting: {family_name} ({family_def['full_name']})")
    print(f"{'='*60}")

    blast_ids = search_blast(blast_df, family_name, family_def)
    hmmer_ids, hmmer_domains = search_hmmer(hmmer_file, family_name, family_def)

    all_ids = blast_ids | hmmer_ids
    print(f"\n  Union (any source): {len(all_ids)} unique genes")

    if not all_ids:
        print(f"  No {family_name} genes found from any source")
        return pd.DataFrame()

    # Build result table
    gene_id_col = 'gene_id' if 'gene_id' in blast_df.columns else blast_df.index.name
    if gene_id_col and gene_id_col in blast_df.columns:
        blast_lookup = blast_df.set_index('gene_id') if gene_id_col == 'gene_id' else blast_df
    else:
        blast_lookup = blast_df

    rows = []
    for gid in sorted(all_ids):
        row = {'gene_id': gid, 'gene_family': family_name}

        # Expression stats from BLAST combined file
        if gid in blast_lookup.index:
            gene_data = blast_lookup.loc[gid]
            if isinstance(gene_data, pd.DataFrame):
                gene_data = gene_data.iloc[0]
            for col in ['baseMean', 'log2FoldChange', 'pvalue', 'padj',
                        'blast_description', 'blast_species', 'sseqid',
                        'pident', 'qcovhsp', 'evalue', 'bitscore']:
                if col in gene_data.index:
                    row[col] = gene_data[col]

        # Evidence tracking
        sources = []
        if gid in blast_ids:
            sources.append('BLAST')
        if gid in hmmer_ids:
            sources.append('HMMER')

        row['evidence_sources'] = '+'.join(sources)
        row['evidence_count'] = len(sources)
        row['confidence'] = 'high' if len(sources) >= 2 else 'low'

        # Domain details
        row['hmmer_domains'] = hmmer_domains.get(gid, '')

        # Direction
        lfc = row.get('log2FoldChange', np.nan)
        padj = row.get('padj', np.nan)
        if pd.notna(padj) and pd.notna(lfc) and padj <= 0.05 and abs(lfc) >= 1.0:
            row['direction'] = 'root_up' if lfc > 0 else 'leaf_up'
        elif pd.notna(padj) and pd.notna(lfc):
            row['direction'] = 'ns'
        else:
            row['direction'] = 'no_data'

        rows.append(row)

    result = pd.DataFrame(rows)

    # Sort by confidence (high first) then padj
    conf_order = {'high': 0, 'medium': 1, 'low': 2}
    result['_sort'] = result['confidence'].map(conf_order)
    result = result.sort_values(['_sort', 'padj'], na_position='last').drop(columns='_sort')

    # Summary
    print(f"\n  Results:")
    for conf in ['high', 'medium', 'low']:
        n = (result['confidence'] == conf).sum()
        print(f"    {conf}: {n} genes")
    for src_combo in result['evidence_sources'].value_counts().items():
        print(f"    {src_combo[0]}: {src_combo[1]} genes")

    return result


def main():
    parser = argparse.ArgumentParser(
        description='Extract and verify CYP/OMT gene families with multi-evidence approach',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('--blast', required=True,
                        help='Combined BLAST+DESeq2 annotated file (from run_combine_filter.sbatch)')
    parser.add_argument('--hmmer', default=None,
                        help='HMMER domtblout file (optional, from run_hmmer.sbatch)')
    parser.add_argument('--species', required=True,
                        help='Species code (DC, DG, MF)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output directory')
    parser.add_argument('--families', nargs='+', default=['CYP', 'OMT'],
                        choices=list(GENE_FAMILIES.keys()),
                        help='Which families to extract (default: CYP OMT)')

    args = parser.parse_args()

    # Load BLAST combined file
    print("=" * 70)
    print(f"GENE FAMILY EXTRACTION: {args.species}")
    print("=" * 70)
    print()
    print(f"Loading BLAST+DESeq2 file: {args.blast}")

    if not Path(args.blast).exists():
        print(f"ERROR: File not found: {args.blast}")
        sys.exit(1)

    blast_df = pd.read_csv(args.blast, sep='\t', comment='#')
    print(f"  Loaded {len(blast_df)} genes")

    # Check optional inputs
    hmmer_file = args.hmmer

    if hmmer_file and not Path(hmmer_file).exists():
        print(f"  WARNING: HMMER file not found: {hmmer_file}")
        print(f"  Proceeding without HMMER results")
        hmmer_file = None

    sources_available = ['BLAST']
    if hmmer_file:
        sources_available.append('HMMER')
    print(f"\n  Evidence sources available: {', '.join(sources_available)}")

    # Create output directory
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Extract each family
    all_results = []
    for family_name in args.families:
        family_def = GENE_FAMILIES[family_name]
        result = extract_family(blast_df, hmmer_file,
                                family_name, family_def)

        if not result.empty:
            # Save individual family file
            family_file = out_dir / f"{args.species}_{family_name}_verified.tsv"
            result.to_csv(family_file, sep='\t', index=False)
            print(f"\n  Saved: {family_file} ({len(result)} genes)")
            all_results.append(result)

    # Save combined file
    if all_results:
        combined = pd.concat(all_results, ignore_index=True)
        combined_file = out_dir / f"{args.species}_CYP_OMT_combined.tsv"

        with open(combined_file, 'w') as f:
            f.write("# ========================================================================\n")
            f.write(f"# Gene Family Extraction: {args.species}\n")
            f.write("# ========================================================================\n")
            f.write(f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"# BLAST input: {args.blast}\n")
            f.write(f"# HMMER input: {hmmer_file or 'not provided'}\n")
            f.write(f"# Evidence sources: {', '.join(sources_available)}\n")
            f.write("#\n")
            f.write("# CONFIDENCE LEVELS:\n")
            f.write("#   high = 2 evidence sources agree (BLAST+HMMER)\n")
            f.write("#   low  = 1 evidence source only\n")
            f.write("#\n")
            f.write("# NOTE: Hit species (blast_species) shows the closest characterized\n")
            f.write("# homolog, NOT contamination. Arabidopsis hits mean your carrot gene\n")
            f.write("# is most similar to that Arabidopsis protein in the database.\n")
            f.write("# ========================================================================\n")

        combined.to_csv(combined_file, sep='\t', index=False, mode='a')
        print(f"\n  Saved combined: {combined_file} ({len(combined)} genes)")

    # Print final summary
    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print()

    if not all_results:
        print("  No gene family members found.")
        return 0

    combined = pd.concat(all_results, ignore_index=True)

    for fam in args.families:
        fam_df = combined[combined['gene_family'] == fam]
        if fam_df.empty:
            print(f"  {fam}: 0 genes found")
            continue

        print(f"  {fam} ({GENE_FAMILIES[fam]['full_name']}):")
        print(f"    Total genes:     {len(fam_df)}")

        for conf in ['high', 'medium', 'low']:
            n = (fam_df['confidence'] == conf).sum()
            if n > 0:
                print(f"    {conf} confidence: {n}")

        sig = fam_df[fam_df['direction'].isin(['root_up', 'leaf_up'])]
        if not sig.empty:
            root_up = (sig['direction'] == 'root_up').sum()
            leaf_up = (sig['direction'] == 'leaf_up').sum()
            print(f"    Significant DE:  {len(sig)} (root_up={root_up}, leaf_up={leaf_up})")

        print()

    # Save summary
    summary_file = out_dir / f"{args.species}_gene_family_summary.txt"
    with open(summary_file, 'w') as f:
        f.write(f"Gene Family Extraction Summary: {args.species}\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Evidence sources: {', '.join(sources_available)}\n")
        f.write("\n")

        for fam in args.families:
            fam_df = combined[combined['gene_family'] == fam]
            f.write(f"{fam} ({GENE_FAMILIES[fam]['full_name']}):\n")
            f.write(f"  Total: {len(fam_df)}\n")
            for conf in ['high', 'medium', 'low']:
                n = (fam_df['confidence'] == conf).sum()
                f.write(f"  {conf}: {n}\n")
            for src, count in fam_df['evidence_sources'].value_counts().items():
                f.write(f"  {src}: {count}\n")
            sig = fam_df[fam_df['direction'].isin(['root_up', 'leaf_up'])]
            root_up = (sig['direction'] == 'root_up').sum()
            leaf_up = (sig['direction'] == 'leaf_up').sum()
            f.write(f"  DE root_up: {root_up}\n")
            f.write(f"  DE leaf_up: {leaf_up}\n")
            f.write(f"  Not significant: {len(fam_df) - len(sig)}\n")
            f.write("\n")

    print(f"  Summary saved: {summary_file}")
    print()
    print("=" * 70)
    print("EXTRACTION COMPLETE!")
    print("=" * 70)

    return 0


if __name__ == '__main__':
    sys.exit(main())
