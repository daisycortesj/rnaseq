#!/usr/bin/env python3
"""
check_mdb_site.py
For every P450 in a protein FASTA (e.g. the TransDecoder .pep file), report:
  1. whether it has the P450 heme-binding signature (PROSITE PS00086)
  2. the residue that lines up with sesamin synthase (SS) Ala308, the
     methylenedioxy-bridge (MDB) site in the I-helix motif (SS: G-T-D-[A]-S)
  3. the I-helix window around that site, and % identity to SS

Interpretation of the SS-308 position:
  A or S  -> MDB-like (SS has A; the nutmeg candidates Mf30795/31663/41180 have S)
  T       -> the usual P450 residue (Thr); in SS, A308T reduced activity
  other   -> unusual; look at it by eye
  -       -> the protein does not reach this region (partial ORF)

Usage:
  python3 check_mdb_site.py --pep proteins.pep --out mdb_site.tsv
  python3 check_mdb_site.py --pep proteins.pep --ids ids.txt --out mdb_site.tsv
  (--ids: one protein ID per line; only those are checked, heme motif or not)
Needs Biopython:  conda install -c conda-forge biopython
"""
import argparse
import re
from Bio import Align
from Bio.Align import substitution_matrices

SS = ("MEAEMLYSALALTFAIFMVYRILSNSQDKRSLTKLPPSPPGWLPVIGHAHLMKNLLHRTLYDFSQKLGPIFSIRFGSRLVVVVSSSSLVEECFTKYDIVLANRPQASVDRR"
      "SLGFSTTSVIGAPYGDHWRNLRKLCDLEVFAPTRLASFLSIRLDERDRMISALYKISSAGFAKVNLEAKIVELTFNNIMRMVAAKRYYGEEAEDDEEAKRFRDLTKEALELT"
      "SASNPGEIFPILRWLGCNGLEKKLAVHSRKTDEFMQGLLDEHRRGERQNTMVDHLLSLQESQPEYYTDEIITGLIVALIIAGTDASVVTTEWAMSLLLNHPKVLEKARKELDTL"
      "VGHERMVDEHDLPKLRYLHCIVLETLRLFPSVPTLVPHEPSEDCKIGGYNVPKGTMVLVNAWAIHRDPKVWDDPLSFKPDRFEIMEVETHKLLPFGMGRRACPGAGLAQKFVG"
      "LALGSLIQCFDWERTSPEKIDLNEGSGITLPKAKTLEAMCKPRHVMEKVLRQVSNV")
SITE = 307                 # 0-based index of SS Ala308
WINDOW = range(302, 313)   # SS residues 303-313: IAGTD[A]SVVTT
assert SS[SITE] == "A" and SS[304:309] == "GTDAS"

HEME = re.compile(r"[FW][SGNH].[GD][^F][RKHPT][^P]C[LIVMFAP][GAD]")  # PROSITE PS00086

parser = argparse.ArgumentParser()
parser.add_argument("--pep", required=True, help="protein FASTA, e.g. TransDecoder .pep")
parser.add_argument("--out", required=True, help="output table (TSV)")
parser.add_argument("--ids", help="optional file with protein IDs to check (one per line)")
parser.add_argument("--min-len", type=int, default=200, help="skip proteins shorter than this (aa)")
args = parser.parse_args()

wanted = None
if args.ids:
    with open(args.ids) as handle:
        wanted = {line.split()[0] for line in handle if line.strip()}

aligner = Align.PairwiseAligner()
aligner.mode = "local"
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5


def read_fasta(path):
    name, chunks = None, []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line)
    if name:
        yield name, "".join(chunks)


def map_ss_positions(seq):
    """Align seq to SS; return {SS index: seq index} and % identity over aligned columns."""
    aln = aligner.align(SS, seq)[0]
    mapping, same, total = {}, 0, 0
    for (s0, s1), (t0, t1) in zip(*aln.aligned):
        for k in range(s1 - s0):
            mapping[s0 + k] = t0 + k
            total += 1
            same += SS[s0 + k] == seq[t0 + k]
    return mapping, (100.0 * same / total if total else 0.0)


def classify(res):
    if res in ("A", "S"):
        return "MDB-like (A/S)"
    if res == "T":
        return "typical P450 (T)"
    if res == "-":
        return "region missing"
    return "other"


checked = 0
class_counts = {}
with open(args.out, "w") as out:
    out.write("protein_id\tlength\theme_motif\tSS308_residue\tclass\tIhelix_window(SS303-313)\tpct_id_to_SS\n")
    for name, seq in read_fasta(args.pep):
        seq = seq.rstrip("*").upper()
        heme = HEME.search(seq)
        if wanted is not None:
            if name not in wanted:
                continue
        elif not heme or len(seq) < args.min_len:
            continue
        mapping, pid = map_ss_positions(seq)
        res = seq[mapping[SITE]] if SITE in mapping else "-"
        window = "".join(seq[mapping[i]] if i in mapping else "-" for i in WINDOW)
        label = classify(res)
        out.write("\t".join([name, str(len(seq)), heme.group(0) if heme else "none",
                             res, label, window, format(pid, ".1f")]) + "\n")
        checked += 1
        class_counts[label] = class_counts.get(label, 0) + 1

print("Checked " + str(checked) + " proteins. Wrote " + args.out)
for label in sorted(class_counts):
    print("  " + str(class_counts[label]) + "  " + label)
