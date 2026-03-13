"""
BLAST wrapper — sequence similarity search (BLASTp and BLASTx).

BLASTp searches protein queries against a protein database.
BLASTx translates nucleotide queries in all six frames and searches
against a protein database.
"""

from typing import List, Optional
import os
import subprocess


OUTFMT_COLS = (
    "qseqid sseqid pident length mismatch gapopen "
    "qstart qend sstart send evalue bitscore stitle"
)

OUTFMT = f"6 {OUTFMT_COLS}"


def check_blast() -> bool:
    """Return True if blastp is available on PATH."""
    try:
        subprocess.run(
            ["blastp", "-version"],
            capture_output=True, check=False
        )
        return True
    except FileNotFoundError:
        return False


def build_blastp_cmd(
    query: str,
    db: str,
    output: str,
    evalue: str = "1e-3",
    max_target_seqs: int = 5,
    threads: int = 4,
    qcov_hsp_perc: int = 0,
    outfmt: str = OUTFMT,
) -> List[str]:
    """
    Build BLASTp command for protein-vs-protein search.

    Parameters
    ----------
    query : str
        Path to query protein FASTA file.
    db : str
        Database name or path (e.g. "swissprot" or full path).
    output : str
        Path for the output TSV file.
    evalue : str
        E-value threshold. "1e-3" = discovery, "1e-10" = strict.
    max_target_seqs : int
        Maximum number of hits per query.
    threads : int
        Number of threads.
    qcov_hsp_perc : int
        Minimum query coverage per HSP (0-100). 0 = no filter.
    outfmt : str
        BLAST output format string.

    Returns
    -------
    list of str
        Command ready for runners.local.run().
    """
    cmd = [
        "blastp",
        "-query", query,
        "-db", db,
        "-out", output,
        "-evalue", str(evalue),
        "-max_target_seqs", str(max_target_seqs),
        "-num_threads", str(threads),
        "-outfmt", outfmt,
    ]

    if qcov_hsp_perc > 0:
        cmd.extend(["-qcov_hsp_perc", str(qcov_hsp_perc)])

    return cmd


def build_blastx_cmd(
    query: str,
    db: str,
    output: str,
    evalue: str = "1e-3",
    max_target_seqs: int = 5,
    threads: int = 4,
    qcov_hsp_perc: int = 0,
    outfmt: str = OUTFMT,
) -> List[str]:
    """
    Build BLASTx command for translated-nucleotide-vs-protein search.

    Parameters are identical to build_blastp_cmd; the only difference
    is the executable (blastx translates nucleotides in all six frames).
    """
    cmd = [
        "blastx",
        "-query", query,
        "-db", db,
        "-out", output,
        "-evalue", str(evalue),
        "-max_target_seqs", str(max_target_seqs),
        "-num_threads", str(threads),
        "-outfmt", outfmt,
    ]

    if qcov_hsp_perc > 0:
        cmd.extend(["-qcov_hsp_perc", str(qcov_hsp_perc)])

    return cmd


def get_completed_query_ids(output_file: str) -> set:
    """
    Parse an existing BLAST output to find which query IDs already have
    results. Used for resume logic when a BLAST job times out.
    """
    completed = set()
    if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
        with open(output_file) as fh:
            for line in fh:
                if line.strip():
                    completed.add(line.split("\t")[0])
    return completed
