#!/usr/bin/env python3
"""
===============================================================================
STEP 1: GENE-FAMILY GFF/GTF KEYWORD SEARCH
===============================================================================

# CHANGED: 2026-03-31 — Renamed from cyp_gff_keyword_search.py → gene_gff_keyword_search.py.
#   Now uses --family arg (cyp or omt) to select keyword set.

Search a GFF3 or GTF annotation file for gene-family candidate genes
using keyword matching on description, product, and gene name fields.

CYP keywords (case-insensitive):
  "cytochrome P450", "CYP", "P450", "monooxygenase"
OMT keywords (case-insensitive):
  "O-methyltransferase", "methyltransferase", "COMT", "CCoAOMT", "FOMT",
  "IOMT", "caffeic acid", "caffeoyl-CoA"

Output:  <family>_gtf_candidates.csv
Columns: gene_id, gene_name, chromosome, start, stop, strand, description,
         matched_keyword

USAGE:
  python scripts/08_gene_families/gene_gff_keyword_search.py --family cyp
  python scripts/08_gene_families/gene_gff_keyword_search.py --family omt

  # Override defaults via command line:
  python scripts/08_gene_families/gene_gff_keyword_search.py \
      --gff /path/to/genomic.gff \
      --output /path/to/omt_gtf_candidates.csv \
      --family omt

NOTE:
  Handles both GFF3 (key=value) and GTF (key "value";) attribute formats.
  "CYP" matching uses word-boundary logic to avoid false positives (e.g.,
  it will NOT match "encryption" or "CYProsin").

Author: Daisy Cortes
===============================================================================
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

# ============================================================================
# USER CONFIGURATION — change these to match your file locations
# ============================================================================

GFF_FILE = "04_reference/dc_genomic.gtf"
OUTPUT_FILE = "07_NRdatabase/cyp450_database/cyp_gtf_candidates.csv"

# CHANGED: 2026-03-31 — Added --family support. Use --family omt for OMT keywords.
# Keywords per family (case-insensitive)
# Each tuple is (pattern, label) — label is stored in the matched_keyword column
FAMILY_KEYWORDS = {
    'cyp': [
        (r'cytochrome\s+P450',   'cytochrome P450'),
        (r'\bCYP\d',             'CYP'),          # CYP followed by a digit (CYP71, CYP76, ...)
        (r'\bCYP\b',             'CYP'),          # Standalone CYP
        (r'\bP450\b',            'P450'),
        (r'monooxygenase',       'monooxygenase'),
    ],
    'omt': [
        (r'O-methyltransferase',          'O-methyltransferase'),
        (r'methyltransferase',            'methyltransferase'),
        (r'caffeic\s+acid.*methyltransf', 'COMT'),
        (r'caffeoyl.CoA.*methyltransf',   'CCoAOMT'),
        (r'\bCOMT\b',                     'COMT'),
        (r'\bCCoAOMT\b',                  'CCoAOMT'),
    ],
}

# Default keywords (backward-compatible)
KEYWORDS = FAMILY_KEYWORDS['cyp']


# ============================================================================
# GFF / GTF ATTRIBUTE PARSERS
# ============================================================================

def parse_gff3_attributes(attr_string):
    """Parse GFF3 attribute string (key=value;key=value)."""
    attrs = {}
    for item in attr_string.strip().split(';'):
        item = item.strip()
        if '=' in item:
            key, value = item.split('=', 1)
            attrs[key] = value.replace('%20', ' ').replace('%2C', ',')
    return attrs


def parse_gtf_attributes(attr_string):
    """Parse GTF attribute string (key "value"; key "value";)."""
    attrs = {}
    for item in attr_string.strip().split(';'):
        item = item.strip()
        if item:
            parts = item.split(' ', 1)
            if len(parts) == 2:
                key, val = parts
                attrs[key] = val.strip('"')
    return attrs


def parse_attributes(attr_string):
    """Auto-detect GFF3 vs GTF format and parse accordingly."""
    if '=' in attr_string and '"' not in attr_string.split('=')[0]:
        return parse_gff3_attributes(attr_string)
    return parse_gtf_attributes(attr_string)


# ============================================================================
# KEYWORD MATCHING
# ============================================================================

def match_keywords(text, keywords=None):
    """Check if text matches any keyword. Returns matched label or None."""
    if not text:
        return None
    if keywords is None:
        keywords = KEYWORDS
    for pattern, label in keywords:
        if re.search(pattern, text, re.IGNORECASE):
            return label
    return None


# ============================================================================
# MAIN GFF PARSER
# ============================================================================

def search_gff_for_cyp(gff_path, keywords=None, family_label="CYP"):
    """
    Parse a GFF/GTF file and extract genes matching keyword patterns.

    Strategy:
      - Read every feature line (gene, mRNA, CDS, transcript, etc.)
      - Check description/product/Name/gene fields for keyword hits
      - Collect gene-level info; deduplicate by gene_id keeping the
        most informative description
    """
    if keywords is None:
        keywords = KEYWORDS
    gff_path = Path(gff_path)
    if not gff_path.exists():
        print(f"ERROR: GFF file not found: {gff_path}")
        sys.exit(1)

    print(f"Parsing GFF/GTF file: {gff_path}")
    print(f"  File size: {gff_path.stat().st_size / 1e6:.1f} MB")

    # gene_id -> dict of best info found so far
    gene_hits = {}
    lines_read = 0
    lines_matched = 0

    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            line = line.strip()
            if not line:
                continue

            lines_read += 1
            fields = line.split('\t')
            if len(fields) < 9:
                continue

            chrom     = fields[0]
            feature   = fields[2]
            start     = fields[3]
            stop      = fields[4]
            strand    = fields[6]
            attr_str  = fields[8]

            attrs = parse_attributes(attr_str)

            # Build a searchable text block from all descriptive fields
            searchable_fields = []
            for key in ['description', 'product', 'Name', 'gene',
                        'gene_biotype', 'note', 'Note']:
                val = attrs.get(key, '')
                if val:
                    searchable_fields.append(val)

            searchable_text = ' | '.join(searchable_fields)
            matched = match_keywords(searchable_text, keywords)

            if not matched:
                continue

            lines_matched += 1

            # Extract gene_id (GFF3 may use gene-LOC..., GTF uses gene_id)
            gene_id = (
                attrs.get('gene_id')
                or attrs.get('gene')
                or attrs.get('Name')
                or attrs.get('ID', '').replace('gene-', '')
            )

            if not gene_id:
                # For CDS/mRNA features in GFF3, gene_id is in the Parent field
                parent = attrs.get('Parent', '')
                if parent:
                    gene_id = parent.replace('gene-', '').split(',')[0]

            if not gene_id:
                continue

            # Clean the gene_id (GFF3 IDs may have prefixes like "rna-", "cds-")
            clean_id = gene_id
            for prefix in ['gene-', 'rna-', 'cds-']:
                if clean_id.startswith(prefix):
                    clean_id = clean_id[len(prefix):]

            gene_name = attrs.get('gene', attrs.get('Name', ''))
            description = attrs.get('description',
                                    attrs.get('product', ''))

            # Store gene-level info (prefer gene features over CDS/mRNA
            # for coordinates, but CDS/mRNA often has richer descriptions)
            if clean_id not in gene_hits:
                gene_hits[clean_id] = {
                    'gene_id':         clean_id,
                    'gene_name':       gene_name,
                    'chromosome':      chrom,
                    'start':           int(start),
                    'stop':            int(stop),
                    'strand':          strand,
                    'description':     description,
                    'matched_keyword': matched,
                }
            else:
                existing = gene_hits[clean_id]
                # Update description if the new one is longer/more informative
                if len(description) > len(existing['description']):
                    existing['description'] = description
                # Update coordinates to gene-level if this is a gene feature
                if feature == 'gene':
                    existing['chromosome'] = chrom
                    existing['start'] = int(start)
                    existing['stop'] = int(stop)
                    existing['strand'] = strand

            if lines_read % 500_000 == 0:
                print(f"  ...processed {lines_read:,} lines, "
                      f"{len(gene_hits)} unique {family_label} genes so far")

    print(f"  Total lines parsed: {lines_read:,}")
    print(f"  Lines with keyword hits: {lines_matched:,}")
    print(f"  Unique {family_label} candidate genes: {len(gene_hits)}")

    return gene_hits


# ============================================================================
# OUTPUT
# ============================================================================

def save_results(gene_hits, output_path):
    """Save CYP candidates to CSV."""
    if not gene_hits:
        print("\nWARNING: No CYP candidates found! Check your GFF file and keywords.")
        pd.DataFrame(columns=[
            'gene_id', 'gene_name', 'chromosome', 'start', 'stop',
            'strand', 'description', 'matched_keyword',
        ]).to_csv(output_path, index=False)
        return pd.DataFrame()

    df = pd.DataFrame(gene_hits.values())

    # Sort by chromosome then start position
    df = df.sort_values(['chromosome', 'start']).reset_index(drop=True)

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)

    print(f"\nSaved: {output_path}")
    print(f"  Rows: {len(df)}")
    print(f"  Columns: {', '.join(df.columns)}")

    return df


# ============================================================================
# SUMMARY
# ============================================================================

def print_summary(df):
    """Print a summary of keyword hits."""
    if df.empty:
        return

    print()
    print("=" * 60)
    print("STEP 1 SUMMARY — GFF KEYWORD SEARCH")
    print("=" * 60)
    print(f"  Total unique CYP candidate genes: {len(df)}")

    if 'matched_keyword' in df.columns:
        print("\n  Hits by keyword:")
        for kw, count in df['matched_keyword'].value_counts().items():
            print(f"    {kw:25s} {count:4d} genes")

    if 'chromosome' in df.columns:
        n_chr = df['chromosome'].nunique()
        print(f"\n  Distributed across {n_chr} chromosomes/scaffolds")

    if 'strand' in df.columns:
        plus  = (df['strand'] == '+').sum()
        minus = (df['strand'] == '-').sum()
        print(f"  Strand: {plus} (+) / {minus} (-)")

    print()
    print("  First 5 candidates:")
    preview_cols = ['gene_id', 'gene_name', 'description']
    cols = [c for c in preview_cols if c in df.columns]
    for _, row in df.head(5).iterrows():
        name = row.get('gene_name', '')
        desc = row.get('description', '')[:60]
        print(f"    {row['gene_id']:25s}  {name:15s}  {desc}")

    print("=" * 60)


# ============================================================================
# CLI
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Step 1: Search GFF/GTF for gene-family candidates by keyword',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        '--gff', default=GFF_FILE,
        help=f'Path to GFF3 or GTF annotation file (default: {GFF_FILE})',
    )
    parser.add_argument(
        '--output', '-o', default=OUTPUT_FILE,
        help=f'Output CSV path (default: {OUTPUT_FILE})',
    )
    parser.add_argument(
        '--family', default='cyp', choices=list(FAMILY_KEYWORDS.keys()),
        help='Gene family to search for (default: cyp). Options: cyp, omt',
    )
    args = parser.parse_args()

    family = args.family.lower()
    family_upper = family.upper()
    keywords = FAMILY_KEYWORDS[family]

    print("=" * 60)
    print(f"{family_upper} DATABASE — STEP 1: GFF KEYWORD SEARCH")
    print("=" * 60)
    print(f"  Family:   {family_upper}")
    print(f"  Keywords: {', '.join(label for _, label in keywords)}")
    print(f"  GFF file: {args.gff}")
    print(f"  Output:   {args.output}")
    print()

    gene_hits = search_gff_for_cyp(args.gff, keywords=keywords, family_label=family_upper)
    df = save_results(gene_hits, args.output)
    print_summary(df)

    return 0


if __name__ == '__main__':
    sys.exit(main())
