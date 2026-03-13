"""
featureCounts wrapper — count aligned reads per gene.

featureCounts is part of the Subread package. It reads BAM files and a GTF
annotation, then counts how many read pairs overlap each gene's exons.
"""

from typing import List
import subprocess


def check_featurecounts() -> bool:
    """Return True if featureCounts is available on PATH."""
    try:
        subprocess.run(
            ["featureCounts", "-v"],
            capture_output=True, check=False
        )
        return True
    except FileNotFoundError:
        return False


def build_featurecounts_cmd(
    gtf: str,
    output_file: str,
    bam_files: List[str],
    threads: int = 4,
    feature_type: str = "exon",
    attribute: str = "gene_id",
    paired_end: bool = True,
    largest_overlap: bool = True,
    strand: int = 0,
) -> List[str]:
    """
    Build featureCounts command for counting reads per gene.

    Parameters
    ----------
    gtf : str
        Path to GTF annotation file.
    output_file : str
        Path for the output counts table.
    bam_files : list of str
        Paths to sorted BAM files.
    threads : int
        Number of threads for parallel counting.
    feature_type : str
        GTF feature type to count (default "exon").
    attribute : str
        GTF attribute to group by (default "gene_id").
    paired_end : bool
        If True, count fragments (read pairs) instead of individual reads.
    largest_overlap : bool
        If True, assign ambiguous reads to the gene with largest overlap.
    strand : int
        Strandedness: 0 = unstranded, 1 = forward, 2 = reverse.

    Returns
    -------
    list of str
        Command ready for runners.local.run().
    """
    cmd = [
        "featureCounts",
        "-a", gtf,
        "-o", output_file,
        "-t", feature_type,
        "-g", attribute,
        "-T", str(threads),
        "-s", str(strand),
    ]

    if paired_end:
        cmd.extend(["-p", "--countReadPairs"])

    if largest_overlap:
        cmd.append("--largestOverlap")

    cmd.extend(bam_files)
    return cmd
