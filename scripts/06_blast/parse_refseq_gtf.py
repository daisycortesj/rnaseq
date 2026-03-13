#!/usr/bin/env python3
"""
===============================================================================
parse_refseq_gtf.py - GTF PARSER WITH AUTO TRANSCRIPT LOOKUP
===============================================================================

Parse RefSeq/Gnomon GTF files and extract gene-level data WITH transcript IDs.

WHAT THIS DOES:
This script extracts gene information from GTF files and automatically looks up
transcript IDs (XM numbers) from transcript/CDS rows. Your PI needs these XM
numbers for NCBI lookups, reproducibility, and publications.

HOW IT WORKS:
1. Pass 1: Extract gene-level features (one row per gene)
2. Pass 2: Look up transcript IDs from transcript/CDS rows
3. Pass 3: Merge together - genes WITH transcript IDs!

OUTPUT FORMAT:
    local_number   - Unique index (1, 2, 3, ...)
    gene_id        - Gene identifier (e.g., LOC135151205)
    transcript_id  - XM number (automatically looked up!)
    description    - Gene/product description

USAGE:
    python parse_refseq_gtf.py --gtf annotation.gtf --out genes_parsed.tsv

EXAMPLES:
    python parse_refseq_gtf.py --gtf 04_reference/dc_genomic.gtf --out genes.tsv

Author: Daisy Cortes
===============================================================================
"""

import argparse
import sys


def parse_gtf_attributes(attr_string):
    """Parse GTF attribute string into a dictionary."""
    attributes = {}
    parts = attr_string.split(';')
    
    for part in parts:
        part = part.strip()
        if not part:
            continue
        
        tokens = part.split(' ', 1)
        if len(tokens) == 2:
            key = tokens[0].strip()
            value = tokens[1].strip()
            if value.startswith('"') and value.endswith('"'):
                value = value[1:-1]
            attributes[key] = value
    
    return attributes


def parse_gtf_line(line):
    """Parse a single GTF line into its components."""
    fields = line.strip().split('\t')
    
    if len(fields) != 9:
        return None
    
    record = {
        'seqname': fields[0],
        'source': fields[1],
        'feature': fields[2],
        'start': fields[3],
        'end': fields[4],
        'score': fields[5],
        'strand': fields[6],
        'frame': fields[7],
    }
    
    attributes = parse_gtf_attributes(fields[8])
    record.update(attributes)
    
    return record


def parse_gtf_smart(gtf_path):
    """
    Smart GTF parser that extracts genes AND looks up their transcript IDs.
    
    Does TWO passes:
    1. Extract all gene-level features
    2. Build gene_id → transcript_id mapping from transcript/CDS rows
    3. Merge them together
    
    Returns:
        - gene_records: List of gene dictionaries with transcript IDs added
        - total_lines: Total data lines read
        - genes_found: Number of genes found
        - transcripts_mapped: Number of genes with transcript IDs
    """
    
    print(f"Reading GTF file: {gtf_path}")
    print("Mode: Smart parsing (genes + transcript lookup)")
    print()
    
    # PASS 1: Extract gene features
    print("Pass 1: Extracting gene-level features...")
    gene_records = []
    total_lines = 0
    
    with open(gtf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            if not line.strip():
                continue
            
            total_lines += 1
            record = parse_gtf_line(line)
            
            if record is None:
                continue
            
            # Extract gene features
            if record['feature'] == 'gene':
                gene_records.append(record)
    
    print(f"  Found {len(gene_records)} genes")
    
    # PASS 2: Build gene_id → transcript_id mapping
    print("Pass 2: Looking up transcript IDs from transcript/CDS rows...")
    gene_to_transcript = {}  # gene_id -> (transcript_id, description)
    
    with open(gtf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            if not line.strip():
                continue
            
            record = parse_gtf_line(line)
            if record is None:
                continue
            
            # Look at transcript and CDS rows (they have transcript_ids)
            if record['feature'] in ['transcript', 'CDS']:
                gene_id = record.get('gene_id', '')
                transcript_id = record.get('transcript_id', '')
                
                # Get description from various possible fields
                description = (record.get('description') or 
                              record.get('product') or 
                              record.get('gene', ''))
                
                if gene_id and transcript_id:
                    # Store first transcript_id found for each gene
                    # (genes can have multiple transcripts, we take the first)
                    if gene_id not in gene_to_transcript:
                        gene_to_transcript[gene_id] = (transcript_id, description)
    
    print(f"  Found transcript IDs for {len(gene_to_transcript)} genes")
    
    # PASS 3: Merge transcript info into gene records
    print("Pass 3: Merging transcript IDs with gene records...")
    transcripts_mapped = 0
    
    for gene in gene_records:
        gene_id = gene.get('gene_id', '')
        
        if gene_id in gene_to_transcript:
            transcript_id, trans_description = gene_to_transcript[gene_id]
            
            # Add transcript_id (overwrite empty one from gene row)
            gene['transcript_id'] = transcript_id
            transcripts_mapped += 1
            
            # If gene row has no description but transcript does, use transcript description
            if not gene.get('description') and trans_description:
                gene['description'] = trans_description
    
    print(f"  Mapped transcript IDs to {transcripts_mapped} genes")
    print(f"  Genes without transcript IDs: {len(gene_records) - transcripts_mapped}")
    
    return gene_records, total_lines, len(gene_records), transcripts_mapped


def extract_output_fields(record, local_number):
    """Extract fields for output in correct order."""
    output_fields = [
        'gene_id',
        'transcript_id',
        'description',
    ]
    
    values = [str(local_number)]
    
    for field in output_fields:
        value = record.get(field, '')
        values.append(value)
    
    return values


def write_tsv_output(records, output_path):
    """Write parsed records to TSV file."""
    headers = [
        'local_number',
        'gene_id',
        'transcript_id',
        'description'
    ]
    
    with open(output_path, 'w') as f:
        f.write('\t'.join(headers) + '\n')
        
        for index, record in enumerate(records, start=1):
            values = extract_output_fields(record, index)
            f.write('\t'.join(values) + '\n')


def count_unique_genes(records):
    """Count unique gene_id values."""
    gene_ids = set()
    
    for record in records:
        gene_id = record.get('gene_id', '')
        if gene_id:
            gene_ids.add(gene_id)
    
    return len(gene_ids)


def main():
    parser = argparse.ArgumentParser(
        description='Smart GTF parser - extracts genes WITH transcript IDs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Extract genes with automatic transcript ID lookup
  python parse_refseq_gtf_smart.py --gtf annotation.gtf --out genes_with_XM.tsv
  
  # This does everything in ONE step:
  # - Extracts gene features (one row per gene)
  # - Automatically finds transcript IDs from transcript/CDS rows
  # - Includes descriptions
  # - Result: Complete gene table with XM numbers!
        """
    )
    
    parser.add_argument(
        '--gtf',
        required=True,
        help='Path to input GTF file'
    )
    
    parser.add_argument(
        '--out',
        required=True,
        help='Path to output TSV file'
    )
    
    args = parser.parse_args()
    
    # Parse the GTF file (smart mode)
    try:
        records, total_lines, genes_found, transcripts_mapped = parse_gtf_smart(args.gtf)
    except FileNotFoundError:
        print(f"Error: GTF file not found: {args.gtf}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading GTF file: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Count unique genes
    unique_genes = count_unique_genes(records)
    
    # Write output
    try:
        write_tsv_output(records, args.out)
    except Exception as e:
        print(f"Error writing output file: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Print summary
    print("\n" + "="*60)
    print("SMART PARSING COMPLETE")
    print("="*60)
    print(f"Total non-comment rows read:     {total_lines}")
    print(f"Genes extracted:                 {genes_found}")
    print(f"Genes with transcript IDs:       {transcripts_mapped}")
    print(f"Genes without transcript IDs:    {genes_found - transcripts_mapped}")
    print(f"Unique gene_id values:           {unique_genes}")
    print(f"\nOutput columns:")
    print(f"  1. local_number (unique index)")
    print(f"  2. gene_id (gene identifier)")
    print(f"  3. transcript_id (XM number) ← AUTOMATICALLY FOUND!")
    print(f"  4. description (gene function)")
    print(f"\nOutput written to:               {args.out}")
    print("="*60)


if __name__ == '__main__':
    main()
