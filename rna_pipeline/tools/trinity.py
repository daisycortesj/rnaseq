
def build_trinity_cmd(left: str, right: str, outdir: str, threads: int, mem_gb: int, resume: bool = False):
    """
    Build Trinity de novo assembly command.
    
    Args:
        left: Path to left/R1 reads (FASTQ)
        right: Path to right/R2 reads (FASTQ)
        outdir: Output directory for assembly
        threads: Number of CPU threads
        mem_gb: Maximum memory in GB
        resume: If True, resume a previous incomplete assembly (Trinity auto-detects checkpoints)
    
    Returns:
        List of command arguments for Trinity
    
    Note:
        Trinity v2.15+ automatically resumes from checkpoints by detecting
        existing .ok files and completed command lists. No special flag needed.
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
    
    # Note: Trinity auto-resumes by detecting checkpoint files (.ok files,
    # recursive_trinity.cmds.completed, etc). No --resume flag needed.
    # The 'resume' parameter is kept for logging/documentation purposes.
    
    return cmd
