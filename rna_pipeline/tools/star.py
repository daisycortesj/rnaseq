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
    """Build STAR command for genome index generation."""
    sj = sjdb_overhang_from_readlen(readlen)
    genome_size = get_genome_size(fasta)
    sa_index_nbases, chr_bin_nbits = get_star_parameters(genome_size)
    
    # Log the parameters being used for debugging
    print(f"[STAR] Genome size: {genome_size:,} bp")
    print(f"[STAR] Using --genomeSAindexNbases {sa_index_nbases} --genomeChrBinNbits {chr_bin_nbits}")
    
    # Base command
    cmd = [
        "STAR",
        "--runThreadN", str(threads),
        "--runMode", "genomeGenerate",
        "--genomeDir", outdir,
        "--genomeFastaFiles", fasta,
    ]
    
    # Add GTF and splice junction parameters only if GTF provided
    if gtf:
        print(f"[STAR] Using GTF annotation: {gtf}")
        cmd.extend([
            "--sjdbGTFfile", gtf,
            "--sjdbOverhang", str(sj),
        ])
    else:
        print("[STAR] WARNING: No GTF provided - building index without splice junction annotations")
    
    # Add genome size parameters
    cmd.extend([
        "--genomeSAindexNbases", str(sa_index_nbases),
        "--genomeChrBinNbits", str(chr_bin_nbits),
    ])
    
    return cmd


def build_star_align_cmd(genome_index: str, reads_left: str, reads_right: str, 
                         outprefix: str, threads: int, quant_mode: bool = True):
    """
    Build STAR command for read alignment.
    
    Parameters
    ----------
    genome_index : str
        Path to STAR genome index directory
    reads_left : str
        Path to R1/forward reads FASTQ file
    reads_right : str
        Path to R2/reverse reads FASTQ file (empty string for single-end)
    outprefix : str
        Output file prefix (e.g., "sample1_" will create sample1_Aligned.out.bam)
    threads : int
        Number of threads to use
    quant_mode : bool, default True
        Whether to quantify gene counts (requires GTF in index)
    
    Returns
    -------
    list
        STAR alignment command as list of strings
    """
    print(f"[STAR] Aligning reads to genome index: {genome_index}")
    print(f"[STAR] Reads: {reads_left}" + (f" + {reads_right}" if reads_right else " (single-end)"))
    
    # Base alignment command
    cmd = [
        "STAR",
        "--runThreadN", str(threads),
        "--runMode", "alignReads",
        "--genomeDir", genome_index,
    ]
    
    # Add read files
    if reads_right:  # Paired-end
        cmd.extend(["--readFilesIn", reads_left, reads_right])
    else:  # Single-end
        cmd.extend(["--readFilesIn", reads_left])
    
    # Handle compressed files
    if reads_left.endswith('.gz'):
        cmd.extend(["--readFilesCommand", "zcat"])
    
    # Output settings
    cmd.extend([
        "--outFileNamePrefix", outprefix,
        "--outSAMtype", "BAM", "SortedByCoordinate",
        "--outSAMunmapped", "Within",  # Include unmapped reads in BAM
        "--outSAMattributes", "Standard",  # Standard SAM attributes
    ])
    
    # Gene quantification (if GTF was used in index)
    if quant_mode:
        print("[STAR] Gene quantification enabled (requires GTF in index)")
        cmd.extend(["--quantMode", "GeneCounts"])
    else:
        print("[STAR] Gene quantification disabled (no GTF in index)")
    
    return cmd
