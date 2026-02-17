#!/usr/bin/env python3
"""
Translate CDS nucleotide sequences to protein sequences.

This script takes CDS FASTA files (from extract_cds_from_gtf.py) and translates
them to protein sequences using the standard genetic code.

Usage:
    python translate_cds_to_protein.py DC
    python translate_cds_to_protein.py DG
    python translate_cds_to_protein.py MF

Author: Daisy Cortes
Date: 2024
"""

import sys
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

def translate_cds_to_protein(input_fasta, output_fasta):
    """
    Translate CDS sequences to proteins.
    
    Args:
        input_fasta: Path to input CDS FASTA file
        output_fasta: Path to output protein FASTA file
    """
    print("=" * 80)
    print("CDS → Protein Translation")
    print("=" * 80)
    print(f"Input:  {input_fasta}")
    print(f"Output: {output_fasta}")
    print()
    
    # Check if input file exists
    if not Path(input_fasta).exists():
        print(f"ERROR: Input file not found: {input_fasta}")
        print("Make sure to run: python scripts/extract_cds_from_gtf.py <SPECIES>")
        sys.exit(1)
    
    # Create output directory if needed
    Path(output_fasta).parent.mkdir(parents=True, exist_ok=True)
    
    # Translate sequences
    translated_count = 0
    warning_count = 0
    
    with open(output_fasta, 'w') as out_handle:
        for record in SeqIO.parse(input_fasta, "fasta"):
            try:
                # Translate CDS to protein
                # to_stop=True stops at first stop codon
                # cds=True validates it's a proper CDS (starts with ATG, ends with stop, length divisible by 3)
                protein_seq = record.seq.translate(to_stop=True, cds=False)
                
                # Check for internal stop codons (indicates potential issues)
                if '*' in str(protein_seq):
                    warning_count += 1
                    if warning_count <= 5:  # Only show first 5 warnings
                        print(f"WARNING: Internal stop codon in {record.id}")
                
                # Create new protein record
                protein_record = record
                protein_record.seq = protein_seq
                protein_record.description = f"{record.description} [translated]"
                
                # Write to output
                SeqIO.write(protein_record, out_handle, "fasta")
                translated_count += 1
                
            except Exception as e:
                print(f"ERROR translating {record.id}: {e}")
                continue
    
    print()
    print("=" * 80)
    print("Translation Summary")
    print("=" * 80)
    print(f"✓ Successfully translated: {translated_count} sequences")
    if warning_count > 0:
        print(f"⚠ Warnings (internal stops): {warning_count} sequences")
    print(f"✓ Output saved to: {output_fasta}")
    print()
    
    return translated_count

def main():
    """Main function."""
    if len(sys.argv) != 2:
        print("Usage: python translate_cds_to_protein.py <SPECIES>")
        print("  SPECIES: DC (Daucus carota), DG (Daucus glochidiatus), or MF (Melaleuca)")
        print()
        print("Examples:")
        print("  python scripts/translate_cds_to_protein.py DC")
        print("  python scripts/translate_cds_to_protein.py DG")
        print("  python scripts/translate_cds_to_protein.py MF")
        sys.exit(1)
    
    species = sys.argv[1].upper()
    
    # Validate species
    if species not in ['DC', 'DG', 'MF']:
        print(f"ERROR: Invalid species '{species}'")
        print("Must be one of: DC, DG, MF")
        sys.exit(1)
    
    # Set up paths
    base_dir = Path("/projects/tholl_lab_1/daisy_analysis/06_analysis")
    input_dir = base_dir / f"blast_input_{species}"
    output_dir = base_dir / f"blast_input_{species}"
    
    input_fasta = input_dir / "all_genes_cds.fasta"
    output_fasta = output_dir / "all_genes_protein.fasta"
    
    # Translate
    translate_cds_to_protein(input_fasta, output_fasta)
    
    print(f"Next step: Run BLASTp using {output_fasta}")
    print(f"  sbatch scripts/blastp_strictfilter.sbatch {species} nr")
    print()

if __name__ == "__main__":
    main()
