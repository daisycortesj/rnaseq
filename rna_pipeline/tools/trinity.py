
def build_trinity_cmd(left: str, right: str, outdir: str, threads: int, mem_gb: int, resume: bool = False, no_normalize_reads: bool = False):
    """
    Build Trinity de novo assembly command.
    
    Args:
        left: Path to left/R1 reads (FASTQ)
        right: Path to right/R2 reads (FASTQ)
        outdir: Output directory for assembly
        threads: Number of CPU threads
        mem_gb: Maximum memory in GB
        resume: If True, resume a previous incomplete assembly (Trinity auto-detects checkpoints)
        no_normalize_reads: If True, skip read normalization (faster, less memory)
    
    Returns:
        List of command arguments for Trinity
    
    Note:
        Trinity v2.15+ automatically resumes from checkpoints by detecting
        existing .ok files and completed command lists. No special flag needed.
        
        This command builds the transcriptome assembly (Trinity.fasta) only.
        For abundance estimation/quantification, run RSEM separately using
        run_trinity_rsem_all.sbatch after assembly completes.
    """
    cmd = [
        "Trinity",
        "--seqType", "fq",
        "--left", left,
        "--right", right,
        "--CPU", str(threads),
        "--max_memory", f"{mem_gb}G",
        "--output", outdir,
    ]
    
    # Optional: Skip read normalization (can help avoid some dependency issues)
    if no_normalize_reads:
        cmd.append("--no_normalize_reads")
    
    # Note: Trinity auto-resumes by detecting checkpoint files (.ok files,
    # recursive_trinity.cmds.completed, etc). No --resume flag needed.
    # The 'resume' parameter is kept for logging/documentation purposes.
    
    return cmd
