#!/usr/bin/env python3
"""
Extract CDS (coding sequence) for genes from GTF + genome FASTA.

Use this when you have genome + GTF but NO protein FASTA.
Creates a nucleotide FASTA file for BLAST (blastx vs NR).

Usage:
  # With species code (DC, DG, MF) - uses default paths:
  python scripts/extract_cds_from_gtf.py DC
  python scripts/extract_cds_from_gtf.py DG
  python scripts/extract_cds_from_gtf.py MF

  # With explicit paths (override defaults):
  python scripts/extract_cds_from_gtf.py <gtf> <genome.fna> <gene_ids.txt> <output.fasta>

  # Custom base directory:
  BASE_DIR=/path/to/project python scripts/extract_cds_from_gtf.py DC
"""

import argparse
import os
import sys
from pathlib import Path

# Species → (GTF, genome, gene_ids_subdir, output_subdir)
# BASE_DIR is read from env or defaults to current project structure
SPECIES_CONFIG = {
    "DC": {
        "gtf": "04_reference/dc_genomic.gtf",
        "genome": "04_reference/GCF_001625215.2_DH1_v3.0_genomic.fna",
        "gene_ids": "06_analysis/pydeseq2_DC_step1_unfiltered/all_gene_ids.txt",
        "output": "06_analysis/blast_input_DC/all_genes_cds.fasta",
    },
    "DG": {
        "gtf": "04_reference/dc_genomic.gtf",
        "genome": "04_reference/GCF_001625215.2_DH1_v3.0_genomic.fna",
        "gene_ids": "06_analysis/pydeseq2_DG_step1_unfiltered/all_gene_ids.txt",
        "output": "06_analysis/blast_input_DG/all_genes_cds.fasta",
    },
    "MF": {
        # Nutmeg: default genome; GTF may not exist - use explicit paths if needed
        "gtf": "04_reference/mf_genomic.gtf",
        "genome": "04_reference/MYU_GWHGDHP00000000.1_genomic.fna",
        "gene_ids": "06_analysis/pydeseq2_MF_step1_unfiltered/all_gene_ids.txt",
        "output": "06_analysis/blast_input_MF/all_genes_cds.fasta",
    },
}


def parse_gtf_attributes(attr_string):
    """Parse GTF attribute string into dict (e.g. gene_id "LOC123"; ...)"""
    attrs = {}
    for item in attr_string.strip().split(';'):
        item = item.strip()
        if item:
            parts = item.split(' ', 1)
            if len(parts) == 2:
                key, val = parts
                attrs[key] = val.strip('"')
    return attrs


def extract_cds_sequences(gtf_file, genome_fasta, gene_ids_file, output_fasta):
    """Extract CDS sequences for specified genes from GTF + genome."""
    
    # Load gene IDs
    print(f"Loading gene IDs from {gene_ids_file}...")
    with open(gene_ids_file) as f:
        target_genes = set(line.strip() for line in f if line.strip())
    print(f"  Target genes: {len(target_genes)}")
    
    # Load genome
    print(f"Loading genome from {genome_fasta}...")
    genome = {}
    current_chr = None
    current_seq = []
    
    with open(genome_fasta) as f:
        for line in f:
            if line.startswith('>'):
                if current_chr and current_seq:
                    genome[current_chr] = ''.join(current_seq)
                current_chr = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.strip())
    if current_chr and current_seq:
        genome[current_chr] = ''.join(current_seq)
    
    print(f"  Chromosomes loaded: {len(genome)}")
    
    # Parse GTF and extract CDS coordinates
    print(f"Parsing GTF from {gtf_file}...")
    gene_cds = {}  # gene_id -> [(chr, start, end, strand)]
    
    with open(gtf_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            feature_type = fields[2]
            if feature_type != 'CDS':
                continue
            chrom = fields[0]
            start = int(fields[3]) - 1  # GTF is 1-based
            end = int(fields[4])
            strand = fields[6]
            attrs = parse_gtf_attributes(fields[8])
            gene_id = attrs.get('gene_id') or attrs.get('gene')
            if not gene_id or gene_id not in target_genes:
                continue
            if gene_id not in gene_cds:
                gene_cds[gene_id] = []
            gene_cds[gene_id].append((chrom, start, end, strand))
    
    print(f"  Genes with CDS found: {len(gene_cds)}")
    
    # Extract sequences and write FASTA
    print(f"Extracting sequences...")
    Path(output_fasta).parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_fasta, 'w') as out:
        for gene_id in sorted(gene_cds.keys()):
            cds_regions = gene_cds[gene_id]
            cds_regions.sort(key=lambda x: (x[0], x[1]))
            strand = cds_regions[0][3]
            seq_parts = []
            for chrom, start, end, _ in cds_regions:
                if chrom not in genome:
                    print(f"  WARNING: Chromosome {chrom} not in genome for {gene_id}")
                    continue
                seq_parts.append(genome[chrom][start:end])
            if not seq_parts:
                continue
            seq = ''.join(seq_parts)
            if strand == '-':
                complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
                seq = ''.join(complement.get(b, 'N') for b in reversed(seq))
            out.write(f">{gene_id}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + '\n')
    
    print(f"✓ Sequences written to {output_fasta}")
    missing = len(target_genes) - len(gene_cds)
    if missing > 0:
        print(f"  Note: {missing} genes had no CDS in GTF (may be non-coding or different ID format)")


def main():
    parser = argparse.ArgumentParser(
        description="Extract CDS sequences for genes from GTF + genome (for BLAST)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "species_or_gtf",
        nargs="?",
        help="Species code (DC, DG, MF) OR path to GTF file (if using 4-arg mode)",
    )
    parser.add_argument("genome", nargs="?", help="Genome FASTA (optional; required in 4-arg mode)")
    parser.add_argument("gene_ids", nargs="?", help="Gene IDs file (optional; required in 4-arg mode)")
    parser.add_argument("output", nargs="?", help="Output FASTA (optional; required in 4-arg mode)")
    parser.add_argument(
        "--base-dir",
        default=os.environ.get("BASE_DIR", "."),
        help="Base directory for relative paths (default: env BASE_DIR or .)",
    )
    args = parser.parse_args()

    base_dir = Path(args.base_dir).resolve()

    # Mode 1: Species code (e.g. extract_cds_from_gtf.py DC)
    if args.species_or_gtf and args.species_or_gtf.upper() in SPECIES_CONFIG and not args.genome:
        species = args.species_or_gtf.upper()
        cfg = SPECIES_CONFIG[species]
        gtf_file = base_dir / cfg["gtf"]
        genome_fasta = base_dir / cfg["genome"]
        gene_ids_file = base_dir / cfg["gene_ids"]
        output_fasta = base_dir / cfg["output"]
        print(f"Species: {species}")
        print(f"  GTF:      {gtf_file}")
        print(f"  Genome:   {genome_fasta}")
        print(f"  Gene IDs: {gene_ids_file}")
        print(f"  Output:   {output_fasta}")
        print()

    # Mode 2: Explicit paths (4 args)
    elif args.species_or_gtf and args.genome and args.gene_ids and args.output:
        gtf_file = Path(args.species_or_gtf)
        genome_fasta = Path(args.genome)
        gene_ids_file = Path(args.gene_ids)
        output_fasta = Path(args.output)

    else:
        parser.print_help()
        print("\nExamples:")
        print("  python scripts/extract_cds_from_gtf.py DC")
        print("  python scripts/extract_cds_from_gtf.py path/to/gtf path/to/genome.fna path/to/gene_ids.txt path/to/output.fasta")
        sys.exit(1)

    for f in [gtf_file, genome_fasta, gene_ids_file]:
        if not Path(f).exists():
            print(f"ERROR: File not found: {f}")
            sys.exit(1)

    extract_cds_sequences(str(gtf_file), str(genome_fasta), str(gene_ids_file), str(output_fasta))


if __name__ == "__main__":
    main()
