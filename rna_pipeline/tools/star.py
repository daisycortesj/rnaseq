
from . import _helpers

def sjdb_overhang_from_readlen(readlen: int) -> int:
    # STAR best practice: overhang = readlen - 1
    return max(1, int(readlen) - 1)

def build_star_index_cmd(fasta: str, gtf: str, outdir: str, threads: int, readlen: int):
    sj = sjdb_overhang_from_readlen(readlen)
    return [
        "STAR",
        "--runThreadN", str(threads),
        "--runMode", "genomeGenerate",
        "--genomeDir", outdir,
        "--genomeFastaFiles", fasta,
        "--sjdbGTFfile", gtf,
        "--sjdbOverhang", str(sj),
    ]
