"""
PROSITE wrapper — short functional motif and active site scanning.

Uses EMBOSS patmatmotifs to scan protein sequences against the PROSITE
database for conserved motifs (active sites, binding sites, post-translational
modification sites).
"""

from typing import List
import subprocess


def check_prosite() -> bool:
    """Return True if EMBOSS patmatmotifs is available on PATH."""
    try:
        subprocess.run(
            ["patmatmotifs", "--help"],
            capture_output=True, check=False
        )
        return True
    except FileNotFoundError:
        return False


def check_prosextract() -> bool:
    """Return True if EMBOSS prosextract is available on PATH."""
    try:
        subprocess.run(
            ["prosextract", "--help"],
            capture_output=True, check=False
        )
        return True
    except FileNotFoundError:
        return False


def build_prosextract_cmd(prosite_dir: str) -> List[str]:
    """
    Build prosextract command to prepare the PROSITE database.

    Must be run once before patmatmotifs. It processes the raw PROSITE
    data files (prosite.dat, prosite.doc) into EMBOSS-readable format.

    Parameters
    ----------
    prosite_dir : str
        Directory containing prosite.dat and prosite.doc files.

    Returns
    -------
    list of str
        Command ready for runners.local.run().
    """
    return [
        "prosextract",
        "-prositedir", prosite_dir,
        "-auto",
    ]


def build_patmatmotifs_cmd(
    query_fasta: str,
    output_file: str,
    full: bool = True,
) -> List[str]:
    """
    Build patmatmotifs command for PROSITE motif scanning.

    Parameters
    ----------
    query_fasta : str
        Path to protein FASTA file to scan.
    output_file : str
        Path for the raw patmatmotifs output.
    full : bool
        If True, report full motif information (recommended).

    Returns
    -------
    list of str
        Command ready for runners.local.run().
    """
    cmd = [
        "patmatmotifs",
        "-sequence", query_fasta,
        "-outfile", output_file,
        "-auto",
    ]

    if full:
        cmd.append("-full")

    return cmd
