"""
HMMER wrapper — protein domain identification using Pfam HMMs.

hmmscan searches protein sequences against a profile HMM database (Pfam)
to identify conserved protein domains and families.
"""

from typing import List
import subprocess


def check_hmmer() -> bool:
    """Return True if hmmscan is available on PATH."""
    try:
        subprocess.run(
            ["hmmscan", "-h"],
            capture_output=True, check=False
        )
        return True
    except FileNotFoundError:
        return False


def build_hmmpress_cmd(pfam_hmm: str) -> List[str]:
    """
    Build hmmpress command to prepare a Pfam database for searching.

    This must be run once before hmmscan can use the database. It creates
    binary index files (.h3m, .h3i, .h3f, .h3p) alongside the .hmm file.

    Parameters
    ----------
    pfam_hmm : str
        Path to Pfam-A.hmm file.

    Returns
    -------
    list of str
        Command ready for runners.local.run().
    """
    return ["hmmpress", pfam_hmm]


def build_hmmscan_cmd(
    pfam_db: str,
    query_fasta: str,
    domtblout: str,
    full_output: str,
    evalue: float = 1e-5,
    threads: int = 4,
    noali: bool = True,
) -> List[str]:
    """
    Build hmmscan command for Pfam domain scanning.

    Parameters
    ----------
    pfam_db : str
        Path to the pressed Pfam-A.hmm database.
    query_fasta : str
        Path to protein FASTA file to scan.
    domtblout : str
        Path for the per-domain tabular output (easy to parse).
    full_output : str
        Path for the full human-readable output.
    evalue : float
        E-value threshold for reporting.
    threads : int
        Number of CPU threads.
    noali : bool
        If True, skip alignment output (faster, smaller files).

    Returns
    -------
    list of str
        Command ready for runners.local.run().
    """
    cmd = [
        "hmmscan",
        "--cpu", str(threads),
        "--domtblout", domtblout,
        "-E", str(evalue),
    ]

    if noali:
        cmd.append("--noali")

    cmd.extend([pfam_db, query_fasta])

    return cmd
