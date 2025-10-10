from . import _helpers
import os
from typing import Tuple  # Import Tuple for Python 3.8 compatibility

def sjdb_overhang_from_readlen(readlen: int) -> int:
    # STAR best practice: overhang = readlen - 1
    return max(1, int(readlen) - 1)

def get_genome_size(fasta: str) -> int:
    """Estimate genome size from FASTA file in base pairs."""
    total_size = 0
    try:
        with open(fasta, 'r') as f:
            for line in f:
                if not line.startswith('>'):
                    total_size += len(line.strip())
    except (IOError, OSError):
        # If we can't read the file, assume it's small
        return 1000000  # 1MB default
    return total_size

def get_star_parameters(genome_size: int) -> Tuple[int, int]:  # Fixed: Tuple instead of tuple
    """
    Get appropriate STAR parameters based on genome size.
    
    Returns:
        Tuple: (genomeSAindexNbases, genomeChrBinNbits)
    """
    if genome_size < 100_000_000:  # < 100MB - small genome
        return (3, 8)
    elif genome_size < 1_000_000_000:  # < 1GB - medium genome  
        return (10, 12)
    else:  # >= 1GB - large genome
        return (13, 14)

def build_star_index_cmd(fasta: str, gtf: str, outdir: str, threads: int, readlen: int):
    sj = sjdb_overhang_from_readlen(readlen)
    genome_size = get_genome_size(fasta)
    sa_index_nbases, chr_bin_nbits = get_star_parameters(genome_size)
    
    # Log the parameters being used for debugging
    print(f"[STAR] Genome size: {genome_size:,} bp")
    print(f"[STAR] Using --genomeSAindexNbases {sa_index_nbases} --genomeChrBinNbits {chr_bin_nbits}")
    
    return [
        "STAR",
        "--runThreadN", str(threads),
        "--runMode", "genomeGenerate",
        "--genomeDir", outdir,
        "--genomeFastaFiles", fasta,
        "--sjdbGTFfile", gtf,
        "--sjdbOverhang", str(sj),
        "--genomeSAindexNbases", str(sa_index_nbases),
        "--genomeChrBinNbits", str(chr_bin_nbits),
    ]
