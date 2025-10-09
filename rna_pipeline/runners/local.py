
# rna_pipeline/runners/local.py
# ----------------------------------------------------------
# PURPOSE:
# A tiny "runner" that executes shell commands for our pipeline
# (e.g., STAR or Trinity) from inside Python.
#
# DESIGN:
# - Print the exact command we're going to run (good for logs).
# - Support a "dry run" mode that only prints the command.
# - Run the command; return 0 on success, or the tool's error code on failure.
#
# COMPATIBILITY NOTE:
# Python 3.8 doesn't support the modern "list[str]" type hint syntax.
# We use typing.List[str] instead so it works on ARC's python 3.8.

from typing import List   # <-- old-style generics that 3.8 understands
import subprocess         # <-- lets Python run external programs (STAR/Trinity)
import sys                # <-- lets us print errors to stderr (log-friendly)


def run(cmd: List[str], dry: bool = False) -> int:
    """
    Run a shell command (like STAR or Trinity) safely.

    Parameters
    ----------
    cmd : List[str]
        The command *split into tokens*. Example:
        ["STAR", "--runThreadN", "8", "--genomeDir", "results/star_index"]
        Why a list? It avoids shell parsing issues and is safer than one big string.
    dry : bool, default False
        If True, we DO NOT execute the command. We only print it.
        Use this to see what would happen without consuming compute.

    Returns
    -------
    int
        Exit status code: 0 means success, anything else means the tool failed.
        We forward the tool's own return code so upstream code (or SLURM)
        can detect failures correctly.

    Behavior
    --------
    1) Print the command we will run (human/audit-friendly).
    2) If dry==True, return 0 immediately (pretend it succeeded).
    3) Otherwise, execute the command. If it fails, catch the error,
       print a readable message to *stderr*, and return the same error code.
    """

    # 1) Show exactly what we are about to run.
    #    " ".join(cmd) turns ["STAR","--runThreadN","8"] into
    #    "STAR --runThreadN 8", which is easier for humans to read in logs.
    print("[run] " + " ".join(cmd), flush=True)

    # 2) Dry-run mode: NO execution; just return success (0).
    #    This is useful for "pipeline smoke tests" and debugging flags/paths.
    if dry:
        return 0

    # 3) Real execution path: try to run the command.
    try:
        # check_call runs the command and raises CalledProcessError if non-zero exit.
        # It streams the tool's stdout/stderr to your terminal/logs live.
        subprocess.check_call(cmd)
        return 0  # success

    except subprocess.CalledProcessError as e:
        # If the tool fails (e.g., STAR can't find a file), we land here.
        # e.returncode is the tool's exit status. By printing to stderr and
        # returning that code, we make failures obvious in SLURM and logs.
        print(
            f"[error] command failed with code {e.returncode}",
            file=sys.stderr,  # stderr = error stream (good for .err files)
            flush=True
        )
        return e.returncode
