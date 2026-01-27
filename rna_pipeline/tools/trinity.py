
def build_trinity_cmd(left: str, right: str, outdir: str, threads: int, mem_gb: int, resume: bool = False):
    """
    Build Trinity de novo assembly command.
    
    Args:
        left: Path to left/R1 reads (FASTQ)
        right: Path to right/R2 reads (FASTQ)
        outdir: Output directory for assembly
        threads: Number of CPU threads
        mem_gb: Maximum memory in GB
        resume: If True, resume a previous incomplete assembly
    
    Returns:
        List of command arguments for Trinity
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
    
    # Resume from checkpoint if requested (for timed-out jobs)
    if resume:
        cmd.append("--resume")
    
    return cmd
