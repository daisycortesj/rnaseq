#!/usr/bin/env python3
"""
===============================================================================
STEP 3: COMBINE GFF KEYWORD + HMMER RESULTS INTO CYP450 MASTER LIST
===============================================================================

Merge the two CYP450 candidate lists from Step 1 (GFF keyword search) and
Step 2 (HMMER domain scan) into a single deduplicated master list.

Each gene is flagged with its evidence source:
  - "Both"       → found by GFF keywords AND confirmed by HMMER (highest confidence)
  - "HMMER_only" → has P450 domain but no keyword match (may have unusual annotation)
  - "GTF_only"   → keyword match but no P450 domain hit (may be pseudogene, fragment,
                    or false positive from keyword matching)

Input files:
  - cyp_gtf_candidates.csv   (from Step 1)
  - cyp_hmmer_confirmed.csv  (from Step 2)

Output:
  - cyp_master_list.csv

IMPORTANT NOTE ON GENE ID MATCHING:
  The GFF file uses gene IDs (e.g., LOC108192212), while the protein FASTA
  from NCBI uses protein accessions (e.g., XP_017313308.1). These are
  different identifiers for the same gene.

  This script tries multiple matching strategies:
    1. Direct gene_id match (works if both files use same IDs)
    2. If protein IDs (XP_*) are in HMMER results, it attempts to build
       a mapping using a GFF/GTF file (if available) or a feature table

  If your IDs don't match, you may need to provide a mapping file
  (see --id-map option).

USAGE:
  python scripts/cyp_combine_results.py

  # Override defaults:
  python scripts/cyp_combine_results.py \
      --gtf-csv cyp_gtf_candidates.csv \
      --hmmer-csv cyp_hmmer_confirmed.csv \
      --output cyp_master_list.csv

  # With an ID mapping file (protein_id → gene_id):
  python scripts/cyp_combine_results.py \
      --id-map protein_to_gene_map.tsv

Author: Daisy Cortes
===============================================================================
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

# ============================================================================
# USER CONFIGURATION
# ============================================================================

GTF_CSV = "07_NRdatabase/cyp450_database/cyp_gtf_candidates.csv"
HMMER_CSV = "07_NRdatabase/cyp450_database/cyp_hmmer_confirmed.csv"
OUTPUT_FILE = "07_NRdatabase/cyp450_database/cyp_master_list.csv"

# GFF file to build protein→gene ID mapping (protein XP_* → gene LOC*)
GFF_FILE = "04_reference/dc_genomic.gtf"


# ============================================================================
# ID MAPPING HELPERS
# ============================================================================

def _parse_attrs_auto(attr_string):
    """Auto-detect GFF3 vs GTF attribute format and parse to dict."""
    attrs = {}
    if '=' in attr_string and '"' not in attr_string.split('=')[0]:
        # GFF3 format: key=value;key=value
        for item in attr_string.split(';'):
            item = item.strip()
            if '=' in item:
                k, v = item.split('=', 1)
                attrs[k] = v.replace('%20', ' ').replace('%2C', ',')
    else:
        # GTF format: key "value"; key "value";
        for item in attr_string.split(';'):
            item = item.strip()
            if item:
                parts = item.split(' ', 1)
                if len(parts) == 2:
                    attrs[parts[0]] = parts[1].strip('"')
    return attrs


def build_protein_to_gene_map_from_gff(gff_path):
    """
    Build a protein_id → gene_id mapping by parsing CDS features.

    Handles both GFF3 and GTF formats:
      - GTF CDS lines have gene_id and protein_id on the same line
      - GFF3 requires chaining: CDS → Parent mRNA → Parent gene
    """
    gff_path = Path(gff_path)
    if not gff_path.exists():
        print(f"  WARNING: GFF/GTF file not found: {gff_path}")
        return {}

    print(f"  Building protein→gene mapping from: {gff_path.name}")
    print(f"  (This may take a minute for large files...)")

    protein_to_gene = {}
    # For GFF3 two-pass approach
    rna_to_gene = {}
    protein_to_rna = {}

    detected_format = None
    lines_read = 0

    with open(gff_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            feature = fields[2]
            attrs = _parse_attrs_auto(fields[8])

            # Detect format on first parsed line
            if detected_format is None:
                if 'gene_id' in attrs and ('transcript_id' in attrs
                                           or 'protein_id' in attrs):
                    detected_format = 'GTF'
                elif 'ID' in attrs or 'Parent' in attrs:
                    detected_format = 'GFF3'
                else:
                    detected_format = 'GTF'
                print(f"  Detected format: {detected_format}")

            if detected_format == 'GTF':
                # GTF: gene_id and protein_id are on the same CDS line
                if feature == 'CDS':
                    gene_id = attrs.get('gene_id', '')
                    protein_id = attrs.get('protein_id', '')
                    if gene_id and protein_id:
                        if protein_id not in protein_to_gene:
                            protein_to_gene[protein_id] = gene_id

            else:
                # GFF3: need to chain CDS → mRNA → gene
                if feature in ('mRNA', 'transcript'):
                    rna_id = attrs.get('ID', '').replace('rna-', '')
                    parent = attrs.get('Parent', '').replace('gene-', '')
                    if rna_id and parent:
                        rna_to_gene[rna_id] = parent
                        raw_rna = attrs.get('ID', '')
                        if raw_rna != rna_id:
                            rna_to_gene[raw_rna] = parent

                elif feature == 'CDS':
                    parent = attrs.get('Parent', '')
                    protein_id = attrs.get('protein_id', '')
                    if not protein_id:
                        dbxref = attrs.get('Dbxref', '')
                        for ref in dbxref.split(','):
                            if 'enbank:' in ref:
                                protein_id = ref.split(':')[1]
                                break
                    if protein_id and parent:
                        clean_parent = parent.replace('rna-', '').split(',')[0]
                        protein_to_rna[protein_id] = clean_parent

            lines_read += 1
            if lines_read % 1_000_000 == 0:
                print(f"  ...processed {lines_read:,} lines")

    # For GFF3, chain the two maps
    if detected_format == 'GFF3':
        for prot_id, rna_id in protein_to_rna.items():
            gene_id = rna_to_gene.get(
                rna_id, rna_to_gene.get('rna-' + rna_id, ''),
            )
            if gene_id:
                protein_to_gene[prot_id] = gene_id

    print(f"  Mapped {len(protein_to_gene)} protein IDs to gene IDs")
    if protein_to_gene:
        sample = list(protein_to_gene.items())[:3]
        for prot, gene in sample:
            print(f"    e.g. {prot} → {gene}")

    return protein_to_gene


def load_id_map(map_path):
    """Load a two-column TSV mapping file (protein_id → gene_id)."""
    df = pd.read_csv(map_path, sep='\t', header=None, names=['protein_id', 'gene_id'])
    return dict(zip(df['protein_id'], df['gene_id']))


# ============================================================================
# LOAD INPUT FILES
# ============================================================================

def load_gtf_csv(path):
    """Load Step 1 GFF keyword results."""
    path = Path(path)
    if not path.exists():
        print(f"WARNING: GTF candidates file not found: {path}")
        print("  (Did you run Step 1 first?)")
        return pd.DataFrame()

    df = pd.read_csv(path)
    print(f"  Loaded GTF candidates: {len(df)} genes from {path}")

    if 'gene_id' not in df.columns:
        print(f"  ERROR: 'gene_id' column not found. Available columns: "
              f"{list(df.columns)}")
        return pd.DataFrame()

    return df


def load_hmmer_csv(path):
    """Load Step 2 HMMER results."""
    path = Path(path)
    if not path.exists():
        print(f"WARNING: HMMER results file not found: {path}")
        print("  (Did you run Step 2 first?)")
        return pd.DataFrame()

    df = pd.read_csv(path)
    print(f"  Loaded HMMER confirmed: {len(df)} genes from {path}")

    if 'gene_id' not in df.columns and 'protein_id' not in df.columns:
        print(f"  ERROR: Neither 'gene_id' nor 'protein_id' column found. "
              f"Available columns: {list(df.columns)}")
        return pd.DataFrame()

    return df


# ============================================================================
# MERGE AND FLAG
# ============================================================================

def reconcile_ids(hmmer_df, protein_to_gene):
    """
    Ensure HMMER results have gene-level IDs that match the GFF results.

    If HMMER gene_ids are protein accessions (XP_*), translate them
    using the mapping.
    """
    if hmmer_df.empty or not protein_to_gene:
        return hmmer_df

    # Check if HMMER gene_ids look like protein accessions
    sample_ids = hmmer_df['gene_id'].head(10).tolist()
    looks_like_protein = any(
        str(gid).startswith(('XP_', 'NP_', 'YP_'))
        for gid in sample_ids
    )

    if not looks_like_protein:
        print("  HMMER gene_ids already appear to be gene-level IDs")
        return hmmer_df

    print("  HMMER gene_ids appear to be protein accessions — translating...")

    # Use the protein_id column if available, otherwise gene_id
    id_col = 'protein_id' if 'protein_id' in hmmer_df.columns else 'gene_id'

    mapped = 0
    unmapped = 0
    new_gene_ids = []

    for _, row in hmmer_df.iterrows():
        prot_id = row[id_col]
        gene_id = protein_to_gene.get(prot_id, '')
        if gene_id:
            new_gene_ids.append(gene_id)
            mapped += 1
        else:
            # Keep the protein ID as fallback
            new_gene_ids.append(prot_id)
            unmapped += 1

    hmmer_df = hmmer_df.copy()
    hmmer_df['gene_id'] = new_gene_ids
    print(f"  Mapped: {mapped}, Unmapped (kept protein ID): {unmapped}")

    return hmmer_df


def combine_results(gtf_df, hmmer_df):
    """
    Combine GFF keyword and HMMER results into one master list.

    Returns a DataFrame with an 'evidence' column indicating how each
    gene was identified.
    """
    gtf_ids = set()
    hmmer_ids = set()

    if not gtf_df.empty and 'gene_id' in gtf_df.columns:
        gtf_ids = set(gtf_df['gene_id'].unique())

    if not hmmer_df.empty and 'gene_id' in hmmer_df.columns:
        hmmer_ids = set(hmmer_df['gene_id'].unique())

    both = gtf_ids & hmmer_ids
    gtf_only = gtf_ids - hmmer_ids
    hmmer_only = hmmer_ids - gtf_ids

    print(f"\n  Overlap analysis:")
    print(f"    GTF keyword hits:     {len(gtf_ids)}")
    print(f"    HMMER domain hits:    {len(hmmer_ids)}")
    print(f"    Found by Both:        {len(both)}")
    print(f"    GTF only:             {len(gtf_only)}")
    print(f"    HMMER only:           {len(hmmer_only)}")

    # Build the master DataFrame
    rows = []

    # Add genes found by both methods
    for gid in sorted(both):
        row = {'gene_id': gid, 'evidence': 'Both'}
        # Pull annotation info from GFF results
        if not gtf_df.empty:
            gtf_row = gtf_df[gtf_df['gene_id'] == gid].iloc[0]
            for col in ['gene_name', 'chromosome', 'start', 'stop', 'strand',
                        'description', 'matched_keyword']:
                if col in gtf_row.index:
                    row[col] = gtf_row[col]
        # Pull HMMER stats
        if not hmmer_df.empty:
            hm_row = hmmer_df[hmmer_df['gene_id'] == gid].iloc[0]
            for col in ['protein_id', 'hmmer_evalue', 'hmmer_score',
                        'hmmer_bias', 'domain_name', 'domain_accession']:
                if col in hm_row.index:
                    row[col] = hm_row[col]
        rows.append(row)

    # Add GTF-only genes
    for gid in sorted(gtf_only):
        row = {'gene_id': gid, 'evidence': 'GTF_only'}
        if not gtf_df.empty:
            gtf_row = gtf_df[gtf_df['gene_id'] == gid].iloc[0]
            for col in ['gene_name', 'chromosome', 'start', 'stop', 'strand',
                        'description', 'matched_keyword']:
                if col in gtf_row.index:
                    row[col] = gtf_row[col]
        rows.append(row)

    # Add HMMER-only genes
    for gid in sorted(hmmer_only):
        row = {'gene_id': gid, 'evidence': 'HMMER_only'}
        if not hmmer_df.empty:
            hm_row = hmmer_df[hmmer_df['gene_id'] == gid].iloc[0]
            for col in ['protein_id', 'hmmer_evalue', 'hmmer_score',
                        'hmmer_bias', 'domain_name', 'domain_accession']:
                if col in hm_row.index:
                    row[col] = hm_row[col]
        rows.append(row)

    master = pd.DataFrame(rows)

    if master.empty:
        print("\n  WARNING: Combined list is empty!")
        return master

    # Sort: Both first, then by chromosome/position if available
    evidence_order = {'Both': 0, 'HMMER_only': 1, 'GTF_only': 2}
    master['_sort'] = master['evidence'].map(evidence_order)

    sort_cols = ['_sort']
    if 'chromosome' in master.columns:
        sort_cols.append('chromosome')
    if 'start' in master.columns:
        sort_cols.append('start')
    sort_cols.append('gene_id')

    master = master.sort_values(sort_cols).reset_index(drop=True)
    master = master.drop(columns=['_sort'])

    return master


# ============================================================================
# OUTPUT AND SUMMARY
# ============================================================================

def save_master_list(master, output_path):
    """Save the combined master list to CSV."""
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    master.to_csv(output_path, index=False)
    print(f"\nSaved: {output_path}")
    print(f"  Total genes: {len(master)}")
    print(f"  Columns: {', '.join(master.columns)}")


def print_summary(master):
    """Print final summary of the CYP450 master list."""
    if master.empty:
        print("\nNo CYP450 genes found by either method.")
        return

    print()
    print("=" * 70)
    print("STEP 3 SUMMARY — CYP450 MASTER LIST")
    print("=" * 70)

    total = len(master)
    both_n = (master['evidence'] == 'Both').sum()
    gtf_n  = (master['evidence'] == 'GTF_only').sum()
    hmm_n  = (master['evidence'] == 'HMMER_only').sum()

    print(f"""
  TOTAL CYP450 GENES:      {total}
  ─────────────────────────────────────
  Both (highest confidence): {both_n:4d}  ({100*both_n/total:.1f}%)
  HMMER only:                {hmm_n:4d}  ({100*hmm_n/total:.1f}%)
  GTF keyword only:          {gtf_n:4d}  ({100*gtf_n/total:.1f}%)
""")

    if 'chromosome' in master.columns:
        n_chr = master['chromosome'].dropna().nunique()
        print(f"  Distributed across {n_chr} chromosomes/scaffolds")

    if 'hmmer_score' in master.columns:
        scored = master.dropna(subset=['hmmer_score'])
        if not scored.empty:
            print(f"\n  HMMER score range (for confirmed hits):")
            print(f"    Min:  {scored['hmmer_score'].min():.1f}")
            print(f"    Max:  {scored['hmmer_score'].max():.1f}")
            print(f"    Mean: {scored['hmmer_score'].mean():.1f}")

    # Show some example genes from each category
    for evidence_type in ['Both', 'HMMER_only', 'GTF_only']:
        subset = master[master['evidence'] == evidence_type]
        if subset.empty:
            continue
        print(f"\n  Example {evidence_type} genes (first 3):")
        for _, row in subset.head(3).iterrows():
            gid  = row['gene_id']
            desc = str(row.get('description', '')).strip()[:55]
            score = row.get('hmmer_score', '')
            score_str = f"  score={score:.1f}" if pd.notna(score) and score != '' else ''
            print(f"    {gid:25s}  {desc}{score_str}")

    print()
    print("=" * 70)
    print("INTERPRETATION GUIDE:")
    print("=" * 70)
    print("""
  'Both' genes are your highest-confidence CYP450s — they have both
  a functional annotation AND a confirmed P450 protein domain.

  'HMMER_only' genes have the P450 domain but weren't caught by keyword
  search. These may have incomplete or generic annotations in the GFF.
  They are still real P450s (the domain doesn't lie).

  'GTF_only' genes matched keywords but lack the P450 domain. Possible
  reasons:
    - Pseudogenes or gene fragments (truncated, no full domain)
    - False positives from keyword matching (e.g., "P450 reductase"
      is NOT a P450 — it's a partner enzyme)
    - The protein wasn't in the FASTA (non-coding annotation)

  RECOMMENDATION: For your CYP450 database, use ALL 'Both' genes,
  ALL 'HMMER_only' genes, and manually review 'GTF_only' genes.
""")
    print("=" * 70)


# ============================================================================
# CLI
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Step 3: Combine GFF keyword + HMMER results into '
                    'CYP450 master list',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        '--gtf-csv', default=GTF_CSV,
        help=f'Step 1 output (default: {GTF_CSV})',
    )
    parser.add_argument(
        '--hmmer-csv', default=HMMER_CSV,
        help=f'Step 2 output (default: {HMMER_CSV})',
    )
    parser.add_argument(
        '--output', '-o', default=OUTPUT_FILE,
        help=f'Output CSV (default: {OUTPUT_FILE})',
    )
    parser.add_argument(
        '--gff', default=GFF_FILE,
        help='GFF file for building protein→gene ID mapping '
             f'(default: {GFF_FILE})',
    )
    parser.add_argument(
        '--id-map', default=None,
        help='Optional TSV mapping file: protein_id<TAB>gene_id',
    )
    args = parser.parse_args()

    print("=" * 70)
    print("CYP450 DATABASE — STEP 3: COMBINE RESULTS")
    print("=" * 70)
    print(f"  GTF candidates:   {args.gtf_csv}")
    print(f"  HMMER confirmed:  {args.hmmer_csv}")
    print(f"  Output:           {args.output}")
    print()

    # Load input files
    print("Loading input files...")
    gtf_df = load_gtf_csv(args.gtf_csv)
    hmmer_df = load_hmmer_csv(args.hmmer_csv)

    if gtf_df.empty and hmmer_df.empty:
        print("\nERROR: Both input files are empty or missing. "
              "Run Steps 1 and 2 first.")
        sys.exit(1)

    # Build or load protein→gene ID mapping
    protein_to_gene = {}
    if args.id_map:
        print(f"\nLoading ID mapping from: {args.id_map}")
        protein_to_gene = load_id_map(args.id_map)
    elif not hmmer_df.empty:
        # Check if HMMER IDs are protein accessions that need translation
        sample = hmmer_df['gene_id'].head(5).tolist()
        needs_mapping = any(
            str(gid).startswith(('XP_', 'NP_', 'YP_'))
            for gid in sample
        )
        if needs_mapping:
            print(f"\nHMMER results use protein accessions — "
                  f"building mapping from GFF...")
            protein_to_gene = build_protein_to_gene_map_from_gff(args.gff)

    # Reconcile IDs
    if not hmmer_df.empty and protein_to_gene:
        print("\nReconciling gene IDs...")
        hmmer_df = reconcile_ids(hmmer_df, protein_to_gene)

    # Combine
    print("\nCombining results...")
    master = combine_results(gtf_df, hmmer_df)

    # Save and summarize
    save_master_list(master, args.output)
    print_summary(master)

    return 0


if __name__ == '__main__':
    sys.exit(main())
