#!/usr/bin/env python3
"""
Extract CYP and OMT gene families for heatmap generation

CYP genes: either
  (a) from a pre-built CYP master list (cyp_master_list.csv from GTF+HMMER
      pipeline), or
  (b) by blastp against a curated P450 reference database
OMT genes: identified by keyword matching on BLAST annotation descriptions

Expression stats (log2FC, padj, etc.) are pulled from the DESeq2-annotated
file (pydeseq2 step1 output or combined BLAST+DESeq2 TSV) for all identified
genes.

Usage:
  # CYP from master list (recommended — uses your 396 P450s, no BLAST):
  python scripts/extract_gene_families.py \
      --blast 06_analysis/pydeseq2_DC_step1_unfiltered/pydeseq2_results_UNFILTERED.tsv \
      --cyp-list 07_NRdatabase/cyp450_database/cyp_master_list.csv \
      --species DC \
      --output 06_analysis/gene_families_DC \
      --families CYP

  # CYP from curated BLAST (legacy):
  python scripts/extract_gene_families.py \
      --blast 06_analysis/combined_DC/DC_swissprot_discovery_annotated.tsv \
      --cyp-db 07_NRdatabase/sukman_database/P450_sequences_8-17-23_Alignment_2.fasta \
      --proteins 04_proteins/DC_proteins.fasta \
      --species DC \
      --output 06_analysis/gene_families_DC

  # OMT only (no CYP):
  python scripts/extract_gene_families.py \
      --blast 06_analysis/combined_DC/DC_swissprot_discovery_annotated.tsv \
      --species DC \
      --output 06_analysis/gene_families_DC \
      --families OMT
"""

import pandas as pd
import numpy as np
import argparse
import re
import sys
import subprocess
import tempfile
from pathlib import Path
from datetime import datetime


# ============================================================================
# FAMILY DEFINITIONS
# ============================================================================

FAMILY_NAMES = {
    'CYP': 'Cytochrome P450',
    'OMT': 'O-Methyltransferase',
}

OMT_BLAST_PATTERNS = [
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
]

OMT_BLAST_EXCLUDE = [
    r'DNA\s+methyltransferase',
    r'histone.*methyltransferase',
    r'rRNA.*methyltransferase',
    r'tRNA.*methyltransferase',
    r'N-methyltransferase',
    r'protein.*methyltransferase',
]

# Legacy config for extract_gene_families_combined.py only.
# This script does NOT use GENE_FAMILIES for CYP: CYP comes from
# --cyp-list (master list) or --cyp-db/--proteins (BLAST).
GENE_FAMILIES = {
    'CYP': {
        'full_name': 'Cytochrome P450',
        'blast_patterns': [],
        'blast_exclude': [],
        'hmmer_pfam': ['PF00067'],
        'hmmer_names': ['p450'],
    },
    'OMT': {
        'full_name': 'O-Methyltransferase',
        'blast_patterns': OMT_BLAST_PATTERNS,
        'blast_exclude': OMT_BLAST_EXCLUDE,
        'hmmer_pfam': ['PF00891', 'PF01596', 'PF08100'],
        'hmmer_names': ['Methyltransf_2', 'Methyltransf_3', 'Dimerisation'],
    },
}


# ============================================================================
# CYP MASTER LIST (from GTF + HMMER pipeline)
# ============================================================================

def load_cyp_master_list(csv_path, evidence_keep=None):
    """Load CYP gene IDs from cyp_master_list.csv (GTF + HMMER pipeline).

    By default keeps only evidence in ('Both', 'HMMER_only') so GTF_only
    (false positives like YUCCA, P450 reductase) are excluded.

    Returns (gene_ids set, hit_info dict). hit_info[gene_id] has keys like
    description, evidence, hmmer_score, gene_name, etc. for build_result_table.
    """
    if evidence_keep is None:
        evidence_keep = {'Both', 'HMMER_only'}
    path = Path(csv_path)
    if not path.exists():
        print(f"  ERROR: CYP master list not found: {path}")
        return set(), {}
    df = pd.read_csv(path)
    if 'gene_id' not in df.columns or 'evidence' not in df.columns:
        print(f"  ERROR: CYP master list must have columns gene_id, evidence")
        return set(), {}
    df = df[df['evidence'].isin(evidence_keep)]
    gene_ids = set(df['gene_id'].dropna().astype(str))
    hit_info = {}
    for _, row in df.iterrows():
        gid = str(row['gene_id'])
        hit_info[gid] = {}
        for col in ['evidence', 'gene_name', 'description', 'chromosome',
                    'start', 'stop', 'strand', 'hmmer_score', 'hmmer_evalue',
                    'protein_id', 'matched_keyword']:
            if col in row.index and pd.notna(row[col]):
                hit_info[gid][col] = row[col]
    print(f"  Loaded {len(gene_ids)} CYP genes from master list "
          f"(evidence in {evidence_keep})")
    return gene_ids, hit_info


# ============================================================================
# CURATED DATABASE SEARCH (CYP)
# ============================================================================

def strip_alignment_gaps(fasta_path):
    """Remove gap characters from an aligned FASTA for use as a BLAST DB.

    Writes a *_clean.fasta next to the original and reuses it if already
    up-to-date.
    """
    src = Path(fasta_path)
    clean_path = src.parent / (src.stem + '_clean.fasta')
    if clean_path.exists() and clean_path.stat().st_mtime >= src.stat().st_mtime:
        return str(clean_path)

    print(f"  Stripping alignment gaps -> {clean_path.name}")
    with open(src) as fin, open(clean_path, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                fout.write(line)
            else:
                fout.write(line.replace('-', '').replace('.', ''))
    return str(clean_path)


def search_curated_db(proteins_fasta, curated_fasta, evalue=1e-10,
                      pident_min=30.0):
    """Identify CYP genes by blastp against a curated P450 reference FASTA.

    Returns (gene_ids, hit_info) where hit_info maps gene_id to a dict of
    best-hit details (curated_hit, curated_pident, curated_evalue,
    curated_bitscore).
    """
    clean_db = strip_alignment_gaps(curated_fasta)

    # Build BLAST DB if needed
    db_extensions = ['.phr', '.pin', '.psq', '.pdb']
    if not any(Path(clean_db + ext).exists() for ext in db_extensions):
        print(f"  Building BLAST database from curated FASTA...")
        cmd = ['makeblastdb', '-in', clean_db, '-dbtype', 'prot']
        try:
            result = subprocess.run(cmd, capture_output=True, text=True,
                                    timeout=60)
        except FileNotFoundError:
            print("  ERROR: makeblastdb not found. Is BLAST+ in your PATH?")
            sys.exit(1)
        if result.returncode != 0:
            print(f"  ERROR: makeblastdb failed:\n{result.stderr.strip()}")
            sys.exit(1)
        print(f"  BLAST database ready")

    # Run blastp
    print(f"  Running blastp: {Path(proteins_fasta).name} vs curated P450 DB")
    with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv',
                                     delete=False) as tmp:
        tmp_out = tmp.name

    cmd = [
        'blastp',
        '-query', str(proteins_fasta),
        '-db', clean_db,
        '-out', tmp_out,
        '-outfmt', '6 qseqid sseqid pident length mismatch gapopen '
                   'qstart qend sstart send evalue bitscore',
        '-evalue', str(evalue),
        '-max_target_seqs', '5',
        '-num_threads', '2',
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True,
                                timeout=600)
    except FileNotFoundError:
        print("  ERROR: blastp not found. Is BLAST+ in your PATH?")
        Path(tmp_out).unlink(missing_ok=True)
        sys.exit(1)

    if result.returncode != 0:
        print(f"  ERROR: blastp failed:\n{result.stderr.strip()}")
        Path(tmp_out).unlink(missing_ok=True)
        sys.exit(1)

    # Parse hits — keep best hit per gene
    gene_ids = set()
    hit_info = {}

    with open(tmp_out) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 12:
                continue
            qseqid = parts[0]
            sseqid = parts[1]
            pident = float(parts[2])
            evalue_val = float(parts[10])
            bitscore = float(parts[11])

            if pident < pident_min:
                continue

            gid = qseqid.split('|')[0].split()[0]
            gene_ids.add(gid)

            if gid not in hit_info:
                hit_info[gid] = {
                    'curated_hit': sseqid,
                    'curated_pident': pident,
                    'curated_evalue': evalue_val,
                    'curated_bitscore': bitscore,
                }

    Path(tmp_out).unlink(missing_ok=True)

    print(f"  Result: {len(gene_ids)} CYP genes "
          f"(evalue < {evalue}, pident > {pident_min}%)")

    return gene_ids, hit_info


# ============================================================================
# REGEX SEARCH (OMT)
# ============================================================================

def search_blast_regex(df, family_name, patterns, exclude_patterns):
    """Search BLAST annotation descriptions for gene family keywords."""
    if 'blast_description' not in df.columns:
        print(f"  WARNING: No blast_description column in annotated file")
        return set()

    desc = df['blast_description'].fillna('')

    include_mask = pd.Series(False, index=df.index)
    for pattern in patterns:
        include_mask |= desc.str.contains(pattern, case=False, regex=True)

    exclude_mask = pd.Series(False, index=df.index)
    for pattern in exclude_patterns:
        exclude_mask |= desc.str.contains(pattern, case=False, regex=True)

    hits = df[include_mask & ~exclude_mask]
    gene_ids = (set(hits['gene_id'].values) if 'gene_id' in hits.columns
                else set(hits.index))

    print(f"  Result: {len(gene_ids)} {family_name} genes "
          f"(matched {include_mask.sum()}, excluded {exclude_mask.sum()})")

    return gene_ids


# Legacy wrappers kept so extract_gene_families_combined.py still imports OK
def search_blast(df, family_name, family_def):
    """Legacy wrapper — delegates to search_blast_regex."""
    return search_blast_regex(df, family_name,
                              family_def['blast_patterns'],
                              family_def['blast_exclude'])


def search_hmmer(hmmer_file, family_name, family_def):
    """Legacy stub — returns empty results (HMMER no longer used)."""
    print(f"  HMMER: disabled (using curated database for CYP)")
    return set(), {}


# ============================================================================
# BUILD RESULT TABLE
# ============================================================================

def build_result_table(gene_ids, deseq_df, family_name, hit_info=None):
    """Build verified gene table with expression stats and hit details."""
    if not gene_ids:
        return pd.DataFrame()

    if 'gene_id' in deseq_df.columns:
        lookup = deseq_df.set_index('gene_id')
    else:
        lookup = deseq_df

    rows = []
    for gid in sorted(gene_ids):
        row = {'gene_id': gid, 'gene_family': family_name}

        if gid in lookup.index:
            gene_data = lookup.loc[gid]
            if isinstance(gene_data, pd.DataFrame):
                gene_data = gene_data.iloc[0]
            for col in ['baseMean', 'log2FoldChange', 'pvalue', 'padj',
                        'blast_description', 'blast_species', 'sseqid',
                        'pident', 'qcovhsp', 'evalue', 'bitscore']:
                if col in gene_data.index:
                    row[col] = gene_data[col]

        if hit_info and gid in hit_info:
            row.update(hit_info[gid])

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
    result = result.sort_values('padj', na_position='last')

    return result


# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Extract CYP (curated DB) and OMT (regex) gene families',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument('--blast', required=True,
                        help='DESeq2 results TSV (or BLAST+DESeq2 combined) '
                             'with gene_id, log2FoldChange, padj, etc.')
    parser.add_argument('--cyp-list', default=None,
                        help='CYP master list CSV (from GTF+HMMER pipeline). '
                             'When set, CYP genes come from this list instead of BLAST.')
    parser.add_argument('--cyp-db', default=None,
                        help='Curated P450 FASTA (required for CYP only if not using --cyp-list)')
    parser.add_argument('--proteins', default=None,
                        help='Species protein FASTA for blastp (required for CYP only if not using --cyp-list)')
    parser.add_argument('--species', required=True,
                        help='Species code (DC, DG, MF)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output directory')
    parser.add_argument('--families', nargs='+', default=['CYP', 'OMT'],
                        choices=['CYP', 'OMT'],
                        help='Which families to extract (default: CYP OMT)')
    parser.add_argument('--evalue', type=float, default=1e-10,
                        help='E-value cutoff for curated DB blastp (default: 1e-10)')
    parser.add_argument('--pident', type=float, default=30.0,
                        help='Min percent identity for curated DB hits (default: 30)')

    args = parser.parse_args()

    # --- Validate inputs ---
    if not Path(args.blast).exists():
        print(f"ERROR: DESeq2/file not found: {args.blast}")
        sys.exit(1)

    if 'CYP' in args.families:
        if args.cyp_list:
            if not Path(args.cyp_list).exists():
                print(f"ERROR: CYP list not found: {args.cyp_list}")
                sys.exit(1)
        else:
            if not args.cyp_db or not args.proteins:
                print("ERROR: For CYP use either --cyp-list OR (--cyp-db and --proteins)")
                sys.exit(1)
            for label, path in [('Curated CYP DB', args.cyp_db),
                                ('Proteins FASTA', args.proteins)]:
                if not Path(path).exists():
                    print(f"ERROR: {label} not found: {path}")
                    sys.exit(1)

    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- Load DESeq2 annotated file ---
    print("=" * 70)
    print(f"GENE FAMILY EXTRACTION: {args.species}")
    print("=" * 70)
    print()
    print(f"DESeq2 / annotated file: {args.blast}")

    deseq_df = pd.read_csv(args.blast, sep='\t', comment='#')
    print(f"  {len(deseq_df)} genes loaded")

    if 'CYP' in args.families:
        if args.cyp_list:
            print(f"CYP source: master list  {args.cyp_list}")
        else:
            print(f"CYP curated database:   {args.cyp_db}")
            print(f"Protein sequences:      {args.proteins}")
    print()

    all_results = []
    cyp_source_note = None  # for combined header / summary

    # ==========  CYP: master list or curated database  ==========
    if 'CYP' in args.families:
        if args.cyp_list:
            print("=" * 60)
            print("CYP (Cytochrome P450) — from master list (GTF + HMMER)")
            print("=" * 60)
            cyp_ids, cyp_hits = load_cyp_master_list(args.cyp_list)
            cyp_source_note = f"master list ({args.cyp_list})"
        else:
            print("=" * 60)
            print("CYP (Cytochrome P450) — curated database search")
            print("=" * 60)
            cyp_ids, cyp_hits = search_curated_db(
                args.proteins, args.cyp_db,
                evalue=args.evalue, pident_min=args.pident,
            )
            cyp_source_note = f"curated database ({args.cyp_db})"

        cyp_result = build_result_table(
            cyp_ids, deseq_df, 'CYP', hit_info=cyp_hits,
        )

        if not cyp_result.empty:
            cyp_file = out_dir / f"{args.species}_CYP_verified.tsv"
            cyp_result.to_csv(cyp_file, sep='\t', index=False)
            print(f"\n  Saved: {cyp_file} ({len(cyp_result)} genes)")
            all_results.append(cyp_result)
        else:
            print("\n  No CYP genes found")

    # ==========  OMT: regex on BLAST descriptions  ==========
    if 'OMT' in args.families:
        print()
        print("=" * 60)
        print("OMT (O-Methyltransferase) — regex search")
        print("=" * 60)

        omt_ids = search_blast_regex(
            deseq_df, 'OMT',
            OMT_BLAST_PATTERNS, OMT_BLAST_EXCLUDE,
        )

        omt_result = build_result_table(omt_ids, deseq_df, 'OMT')

        if not omt_result.empty:
            omt_file = out_dir / f"{args.species}_OMT_verified.tsv"
            omt_result.to_csv(omt_file, sep='\t', index=False)
            print(f"\n  Saved: {omt_file} ({len(omt_result)} genes)")
            all_results.append(omt_result)
        else:
            print("\n  No OMT genes found")

    # ==========  Combined output  ==========
    if all_results:
        combined = pd.concat(all_results, ignore_index=True)
        combined_file = out_dir / f"{args.species}_CYP_OMT_combined.tsv"

        with open(combined_file, 'w') as f:
            f.write("# ======================================"
                    "==========================================\n")
            f.write(f"# Gene Family Extraction: {args.species}\n")
            f.write("# ======================================"
                    "==========================================\n")
            f.write(f"# Generated: "
                    f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"# DESeq2 annotated file: {args.blast}\n")
            if 'CYP' in args.families and cyp_source_note:
                f.write(f"# CYP source: {cyp_source_note}\n")
            f.write("#\n")
            f.write("# CYP: from master list (GTF+HMMER) or curated BLAST\n")
            f.write("# OMT: identified by keyword matching on BLAST "
                    "descriptions\n")
            f.write("# ======================================"
                    "==========================================\n")

        combined.to_csv(combined_file, sep='\t', index=False, mode='a')
        print(f"\n  Saved combined: {combined_file} ({len(combined)} genes)")

    # ==========  Summary  ==========
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
        full_name = FAMILY_NAMES.get(fam, fam)
        fam_df = combined[combined['gene_family'] == fam]
        if fam_df.empty:
            print(f"  {fam}: 0 genes found")
            continue

        print(f"  {fam} ({full_name}):")
        print(f"    Total genes: {len(fam_df)}")

        sig = fam_df[fam_df['direction'].isin(['root_up', 'leaf_up'])]
        if not sig.empty:
            root_up = (sig['direction'] == 'root_up').sum()
            leaf_up = (sig['direction'] == 'leaf_up').sum()
            print(f"    Significant DE: {len(sig)} "
                  f"(root_up={root_up}, leaf_up={leaf_up})")

        ns = (fam_df['direction'] == 'ns').sum()
        no_data = (fam_df['direction'] == 'no_data').sum()
        if ns > 0:
            print(f"    Not significant: {ns}")
        if no_data > 0:
            print(f"    No DE data: {no_data}")
        print()

    summary_file = out_dir / f"{args.species}_gene_family_summary.txt"
    with open(summary_file, 'w') as f:
        f.write(f"Gene Family Extraction Summary: {args.species}\n")
        f.write(f"Generated: "
                f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        if 'CYP' in args.families and cyp_source_note:
            f.write(f"CYP method: {cyp_source_note}\n")
        f.write(f"OMT method: regex on BLAST descriptions\n")
        f.write("\n")

        for fam in args.families:
            full_name = FAMILY_NAMES.get(fam, fam)
            fam_df = combined[combined['gene_family'] == fam]
            f.write(f"{fam} ({full_name}):\n")
            f.write(f"  Total: {len(fam_df)}\n")
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
