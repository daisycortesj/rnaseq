#!/usr/bin/env python3
"""
summarize_read_hits.py
Summarize DIAMOND blastx read hits in ONE comparable unit, so different
searches and datasets can be compared fairly.

Rules (the same for every dataset):
  1. Only alignments of at least --min-len amino acids count (default 30).
     Short 13-15 aa matches in conserved regions are ignored.
  2. The unit is a FRAGMENT, not a read or a hit row. Read 1 and read 2 of a
     pair share one name, so a pair is counted once.
  3. Each fragment is counted ONCE, for the single protein it matches best
     (highest bit score), at that alignment's % identity.

Identity bins:
  95-100%  the gene itself
  70-95%   close relatives
  35-70%   more distant relatives (same P450 family)
  <35%     weak / distant

Usage:
  python3 summarize_read_hits.py --hits-dir DIR --queries proteins.faa --label "Our reads"
  python3 summarize_read_hits.py --hits-dir DIR --queries proteins.faa --label "Ibdah" --tissue leaf
"""
import argparse
import glob
import os
import re

parser = argparse.ArgumentParser()
parser.add_argument("--hits-dir", required=True, help="folder with *_hits.tsv files")
parser.add_argument("--queries", required=True, help="the protein FASTA that was searched")
parser.add_argument("--label", required=True, help="name of the dataset for the report")
parser.add_argument("--tissue", default="auto",
                    help="auto: MFF = fruit, MFL = leaf. Or a fixed name, e.g. leaf")
parser.add_argument("--min-len", type=int, default=30, help="minimum alignment length (aa)")
args = parser.parse_args()

# Protein names, in FASTA order (first word after ">")
proteins = []
with open(args.queries) as handle:
    for line in handle:
        if line.startswith(">"):
            proteins.append(line[1:].split()[0])


def tissue_of(sample):
    if args.tissue != "auto":
        return args.tissue
    if "MFF" in sample:
        return "fruit"
    if "MFL" in sample:
        return "leaf"
    return "other"


def library_of(sample):
    # MFF1_1_filtered -> MFF1_filtered ; DRR513061_2 -> DRR513061
    return re.sub(r"_[12](_|$)", r"\1", sample)


# best[(library, fragment)] = (bitscore, protein, identity, tissue)
best = {}
hit_paths = sorted(glob.glob(os.path.join(args.hits_dir, "*_hits.tsv")))
rows_read = 0
rows_kept = 0

for path in hit_paths:
    sample = os.path.basename(path)[: -len("_hits.tsv")]
    tissue = tissue_of(sample)
    library = library_of(sample)
    with open(path) as handle:
        for line in handle:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                continue
            rows_read += 1
            # columns: read, protein, %identity, length, evalue, bitscore
            if float(cols[3]) < args.min_len:
                continue
            rows_kept += 1
            fragment = re.sub(r"/[12]$", "", cols[0])
            key = (library, fragment)
            record = (float(cols[5]), cols[1], float(cols[2]), tissue)
            if key not in best or record[0] > best[key][0]:
                best[key] = record

bins = [("95_100", 95, 101), ("70_95", 70, 95), ("35_70", 35, 70), ("below35", 0, 35)]
tissues = sorted({rec[3] for rec in best.values()})
if len(tissues) > 1:
    tissues.append("all")

lines = []
lines.append("Dataset: " + args.label)
lines.append("Hit files: " + str(len(hit_paths)))
lines.append("Hit rows read: " + str(rows_read) + "   rows >= " + str(args.min_len) + " aa: " + str(rows_kept))
lines.append("Unit: unique FRAGMENTS (read pairs), each counted once for its best protein,")
lines.append("      alignments >= " + str(args.min_len) + " aa only.")
lines.append("")
header = "protein".ljust(14) + "tissue".ljust(8) + "fragments".rjust(11)
for name, low, high in bins:
    header += ("f_" + name).rjust(11)
header += "max_id".rjust(9)
lines.append(header)

table_rows = []
for protein in proteins:
    for tissue in tissues:
        ids = [rec[2] for rec in best.values()
               if rec[1] == protein and (tissue == "all" or rec[3] == tissue)]
        row = protein.ljust(14) + tissue.ljust(8) + str(len(ids)).rjust(11)
        counts = []
        for name, low, high in bins:
            n = sum(1 for x in ids if low <= x < high)
            counts.append(n)
            row += str(n).rjust(11)
        max_id = max(ids) if ids else 0.0
        row += format(max_id, ".1f").rjust(9)
        lines.append(row)
        table_rows.append([protein, tissue, len(ids)] + counts + [round(max_id, 1)])

lines.append("")
lines.append("f_95_100 = the gene itself. 70-95 = close relatives. 35-70 = P450 family relatives.")

text = "\n".join(lines) + "\n"
with open(os.path.join(args.hits_dir, "SUMMARY_fragments.txt"), "w") as handle:
    handle.write(text)
with open(os.path.join(args.hits_dir, "summary_fragments.tsv"), "w") as handle:
    handle.write("protein\ttissue\tfragments\tf_95_100\tf_70_95\tf_35_70\tf_below35\tmax_identity\n")
    for row in table_rows:
        handle.write("\t".join(str(x) for x in row) + "\n")
print(text)
print("Wrote " + os.path.join(args.hits_dir, "SUMMARY_fragments.txt"))
