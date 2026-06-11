#!/usr/bin/env python3
# =============================================================================
# SCRIPT:  extract_longest_isoform.py
# PURPOSE: Keep only the longest transcript per Trinity gene locus.
#
# INPUT:
#   --fasta      Pooled Trinity FASTA (all isoforms)
#   --gene-map   Trinity .gene_trans_map file (gene_id <tab> transcript_id)
#
# OUTPUT:
#   --output     FASTA with one transcript per gene (the longest one)
#
# WHY:
#   Trinity makes many isoforms per gene. BUSCO "Duplicated (D)" counts genes
#   found more than once — isoforms inflate D. This file lets you run a third
#   BUSCO benchmark with isoform diversity removed (gene families still count).
#
# EXAMPLE:
#   python scripts/04_busco/extract_longest_isoform.py \
#       --fasta  MF_trinity_pooled.Trinity.fasta \
#       --gene-map MF_trinity_pooled.Trinity.fasta.gene_trans_map \
#       --output MF_longest_isoform.fasta
# =============================================================================

import argparse
import sys


def read_fasta_lengths(fasta_path):
    """
    Read a FASTA file and return a dictionary:
      transcript_id -> sequence length (number of bases)
    """
    lengths = {}
    current_id = None
    current_length = 0

    with open(fasta_path, "r") as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                # Save the previous transcript before starting a new one
                if current_id is not None:
                    lengths[current_id] = current_length

                # Header example: >TRINITY_DN0_c0_g1_i1 len=1234 ...
                # We only need the ID (first word after >)
                header = line[1:]
                current_id = header.split()[0]
                current_length = 0
            else:
                current_length = current_length + len(line)

    # Don't forget the last transcript in the file
    if current_id is not None:
        lengths[current_id] = current_length

    return lengths


def read_gene_trans_map(map_path):
    """
    Read Trinity gene_trans_map and return:
      gene_id -> list of transcript IDs belonging to that gene
    """
    gene_to_transcripts = {}

    with open(map_path, "r") as map_file:
        for line in map_file:
            line = line.strip()
            if not line:
                continue

            parts = line.split("\t")
            if len(parts) < 2:
                print(f"WARNING: skipping malformed line: {line}", file=sys.stderr)
                continue

            gene_id = parts[0]
            transcript_id = parts[1]

            if gene_id not in gene_to_transcripts:
                gene_to_transcripts[gene_id] = []

            gene_to_transcripts[gene_id].append(transcript_id)

    return gene_to_transcripts


def pick_longest_transcript_per_gene(gene_to_transcripts, transcript_lengths):
    """
    For each gene, find the transcript with the greatest length.
    Returns the set of transcript IDs to keep.
    """
    keep_ids = set()
    missing_count = 0

    for gene_id, transcript_list in gene_to_transcripts.items():
        best_transcript = None
        best_length = -1

        for transcript_id in transcript_list:
            length = transcript_lengths.get(transcript_id)
            if length is None:
                missing_count = missing_count + 1
                continue

            if length > best_length:
                best_length = length
                best_transcript = transcript_id

        if best_transcript is not None:
            keep_ids.add(best_transcript)

    return keep_ids, missing_count


def write_filtered_fasta(fasta_path, output_path, keep_ids):
    """
    Write a new FASTA containing only transcripts whose IDs are in keep_ids.
    """
    written = 0
    writing = False
    current_id = None

    with open(fasta_path, "r") as in_file, open(output_path, "w") as out_file:
        for line in in_file:
            if line.startswith(">"):
                header = line[1:].strip()
                current_id = header.split()[0]
                writing = current_id in keep_ids
                if writing:
                    out_file.write(line)
                    written = written + 1
            elif writing:
                out_file.write(line)

    return written


def main():
    parser = argparse.ArgumentParser(
        description="Extract longest Trinity isoform per gene for BUSCO benchmarking."
    )
    parser.add_argument("--fasta", required=True, help="Input Trinity FASTA")
    parser.add_argument("--gene-map", required=True, help="Trinity gene_trans_map file")
    parser.add_argument("--output", required=True, help="Output FASTA path")
    args = parser.parse_args()

    print("Reading transcript lengths from FASTA...")
    transcript_lengths = read_fasta_lengths(args.fasta)
    print(f"  Transcripts in FASTA: {len(transcript_lengths)}")

    print("Reading gene-to-transcript map...")
    gene_to_transcripts = read_gene_trans_map(args.gene_map)
    print(f"  Genes in map: {len(gene_to_transcripts)}")

    print("Picking longest isoform per gene...")
    keep_ids, missing = pick_longest_transcript_per_gene(
        gene_to_transcripts, transcript_lengths
    )
    print(f"  Transcripts to keep: {len(keep_ids)}")
    if missing > 0:
        print(f"  WARNING: {missing} map entries had no matching FASTA sequence")

    print(f"Writing output: {args.output}")
    written = write_filtered_fasta(args.fasta, args.output, keep_ids)
    print(f"  Wrote {written} sequences")
    print("Done.")


if __name__ == "__main__":
    main()
