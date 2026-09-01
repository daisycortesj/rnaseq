#!/usr/bin/env python3
"""
Keep one representative gene per Trinity locus for heatmaps.

WHAT THIS DOES:
  Reads a DE candidate TSV (e.g. cyp_candidates_MF.tsv).
  Looks up every isoform in the Trinity FASTA.
  For each locus group, keeps the gene whose LONGEST isoform is the longest.
  Writes a NEW smaller TSV. The original candidate file is not changed.

WHY:
  The heatmap is already gene-level (no _i in the row names). Isoform counts
  were summed into one gene before DESeq2. Extra rows are usually Trinity
  splitting one biological gene into several IDs:

      TRINITY_DN1245_c0_g1     \
      TRINITY_DN1245_c0_g3      }  same component (c0) → keep 1
      TRINITY_DN1245_c1_g1     }  different component (c1) → keep 1 more
      TRINITY_DN1245_c1_g4     /

  Default grouping is "component" (TRINITY_DN###_c#). That is the usual
  "longest isoform as the representative" step for de novo assemblies.

INPUTS:
  --candidates   TSV from run_filter_genelist (must have gene_id or first column)
  --fasta        Trinity / CD-HIT assembly FASTA (isoform headers with _i)

OUTPUT (replacement heatmap table — original TSV is not overwritten):
  cyp_candidates_MF_longest.tsv           ← USE THIS for plots
  cyp_candidates_MF_longest_dropped.tsv   genes removed (shorter fragments)
  cyp_candidates_MF.tsv                   unchanged full DE list (~304 genes)

HOW TO RUN (on ARC, from the repo):
  sbatch scripts/08_gene_families/run_collapse_longest_isoform.sbatch CYP
  sbatch scripts/08_gene_families/run_collapse_longest_isoform.sbatch OMT
  sbatch scripts/08_gene_families/run_collapse_longest_isoform.sbatch both

  # or on a login node:
  conda activate rnaseq
  python scripts/08_gene_families/collapse_longest_isoform.py --family CYP

  Then replot with the REPLACEMENT file (not the 304-gene TSV):
    sbatch scripts/05_pydeseq2/run_step3_plots.sbatch MF \\
      /projects/tholl_lab_1/daisy_analysis/07_NRdatabase/00_2_MF/cyp_candidates_MF_longest.tsv
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd


DEFAULT_BASE = Path("/projects/tholl_lab_1/daisy_analysis")
DEFAULT_CAND_DIR = DEFAULT_BASE / "07_NRdatabase" / "00_2_MF"
DEFAULT_FASTA = (
    DEFAULT_BASE / "01_processed" / "00_6_cdhit" / "MF_trinity_cdhit95.fasta"
)


def to_gene_id(raw_id):
    """Turn a FASTA header ID into a Trinity GENE id (drop _i and .p).

    Examples:
      TRINITY_DN123_c0_g1_i2.p1  ->  TRINITY_DN123_c0_g1
      TRINITY_DN123_c0_g1_i2     ->  TRINITY_DN123_c0_g1
      TRINITY_DN123_c0_g1        ->  TRINITY_DN123_c0_g1
    """
    text = str(raw_id).strip().split()[0]
    # Strip TransDecoder ORF suffix (.p1, .p2, ...)
    text = re.sub(r"\.p[0-9]+$", "", text)
    # Strip Trinity isoform suffix (_i1, _i2, ...)
    text = re.sub(r"_i[0-9]+$", "", text)
    return text


def group_key(gene_id, mode):
    """Decide which genes are treated as the same locus.

    gene:       TRINITY_DN1245_c0_g1  (no collapse — 1 gene = 1 row)
    component:  TRINITY_DN1245_c0     (collapse g1/g2/g3 in the same c#)
    dn:         TRINITY_DN1245        (collapse every gene under that DN)
    """
    gene_id = str(gene_id)
    if mode == "gene":
        return gene_id
    if mode == "dn":
        match = re.match(r"(TRINITY_DN[0-9]+)", gene_id)
        if match:
            return match.group(1)
        return gene_id
    # default: component
    match = re.match(r"(TRINITY_DN[0-9]+_c[0-9]+)", gene_id)
    if match:
        return match.group(1)
    return gene_id


def load_longest_isoform_per_gene(fasta_path):
    """Walk the FASTA once. Return {gene_id: (isoform_id, length)}.

    Trinity headers usually look like:
      >TRINITY_DN1245_c0_g1_i1 len=5299 path=[...]

    We trust len= when it is present (faster). Otherwise we count A/T/G/C.
    If several isoforms exist, we keep the longest one per gene.
    """
    fasta_path = Path(fasta_path)
    if not fasta_path.is_file():
        print(f"ERROR: FASTA not found: {fasta_path}", file=sys.stderr)
        sys.exit(1)

    best = {}
    current_id = None
    current_len_from_header = None
    current_seq_len = 0

    def finish_record():
        # Save the record we just finished (called when a new '>' starts,
        # and again at the end of the file).
        if current_id is None:
            return
        if current_len_from_header is not None:
            seq_len = current_len_from_header
        else:
            seq_len = current_seq_len
        gene_id = to_gene_id(current_id)
        old = best.get(gene_id)
        # Keep this isoform if we have never seen the gene, or it is longer.
        if old is None or seq_len > old[1]:
            best[gene_id] = (current_id, seq_len)

    print(f"  Reading FASTA (this can take a minute): {fasta_path}")
    with open(fasta_path) as handle:
        for line in handle:
            if line.startswith(">"):
                finish_record()
                header = line[1:].strip()
                current_id = header.split()[0]
                current_seq_len = 0
                # Pull len=5299 from the Trinity header if it is there.
                len_match = re.search(r"\blen=([0-9]+)\b", header)
                if len_match:
                    current_len_from_header = int(len_match.group(1))
                else:
                    current_len_from_header = None
            else:
                # Only count bases if the header did not give us a length.
                if current_len_from_header is None:
                    current_seq_len += len(line.strip())
        finish_record()

    print(f"  Genes with at least one isoform in FASTA: {len(best)}")
    return best


def collapse_table(candidates_path, fasta_path, out_path, dropped_path, mode):
    """Read candidates, keep one gene per group, write two TSVs."""
    candidates_path = Path(candidates_path)
    if not candidates_path.is_file():
        print(f"ERROR: candidate TSV not found: {candidates_path}", file=sys.stderr)
        sys.exit(1)

    table = pd.read_csv(candidates_path, sep="\t")
    if "gene_id" not in table.columns:
        table = table.rename(columns={table.columns[0]: "gene_id"})
    table["gene_id"] = table["gene_id"].astype(str)

    longest = load_longest_isoform_per_gene(fasta_path)

    # Add isoform info next to each candidate gene.
    iso_ids = []
    iso_lens = []
    for gene_id in table["gene_id"]:
        hit = longest.get(gene_id)
        if hit is None:
            iso_ids.append("")
            iso_lens.append(0)
        else:
            iso_ids.append(hit[0])
            iso_lens.append(hit[1])
    table["longest_isoform"] = iso_ids
    table["longest_isoform_len"] = iso_lens
    table["locus_group"] = table["gene_id"].map(lambda g: group_key(g, mode))

    missing = (table["longest_isoform_len"] == 0).sum()
    if missing:
        print(f"  WARNING: {missing} genes had no sequence in the FASTA")
        print("           (they lose length ties and may be dropped)")

    # Tie-break: longer isoform first, then smaller padj, then larger |log2FC|.
    table["_abs_lfc"] = 0.0
    if "log2FoldChange" in table.columns:
        table["_abs_lfc"] = table["log2FoldChange"].abs().fillna(0)
    if "padj" in table.columns:
        table["_padj_sort"] = table["padj"].fillna(1.0)
    else:
        table["_padj_sort"] = 1.0

    table = table.sort_values(
        ["locus_group", "longest_isoform_len", "_padj_sort", "_abs_lfc"],
        ascending=[True, False, True, False],
    )
    # First row in each group is the winner (already sorted).
    keep_mask = ~table.duplicated(subset=["locus_group"], keep="first")
    kept = table[keep_mask].copy()
    dropped = table[~keep_mask].copy()

    kept = kept.drop(columns=["_abs_lfc", "_padj_sort"])
    dropped = dropped.drop(columns=["_abs_lfc", "_padj_sort"])

    out_path = Path(out_path)
    dropped_path = Path(dropped_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    kept.to_csv(out_path, sep="\t", index=False)
    dropped.to_csv(dropped_path, sep="\t", index=False)

    print(f"  Grouping mode: {mode}")
    print(f"  Before: {len(table)} DE genes")
    print(f"  After:  {len(kept)} genes  (one per locus group)")
    print(f"  Dropped: {len(dropped)} genes  →  {dropped_path.name}")
    print(f"  Kept:    {out_path}")
    print()
    print("  Example (first 5 kept genes):")
    show = kept.head(5)
    for _, row in show.iterrows():
        print(
            f"    {row['gene_id']}  "
            f"iso={row['longest_isoform'] or 'NA'}  "
            f"len={int(row['longest_isoform_len'])}  "
            f"group={row['locus_group']}"
        )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Collapse DE candidate genes to the longest-isoform representative "
            "per Trinity locus (original TSV is not modified)"
        )
    )
    parser.add_argument(
        "--family",
        choices=["CYP", "OMT", "both"],
        default="CYP",
        help="Which candidate table to collapse (default: CYP)",
    )
    parser.add_argument(
        "--candidates",
        default="",
        help="Path to a candidate TSV (overrides --family default)",
    )
    parser.add_argument(
        "--fasta",
        default=str(DEFAULT_FASTA),
        help="Trinity/CD-HIT FASTA used to measure isoform lengths",
    )
    parser.add_argument(
        "--group",
        choices=["gene", "component", "dn"],
        default="component",
        help=(
            "How to group genes: "
            "gene = no collapse; "
            "component = TRINITY_DN###_c# (default); "
            "dn = whole TRINITY_DN### cluster"
        ),
    )
    parser.add_argument(
        "--out-dir",
        default=str(DEFAULT_CAND_DIR),
        help="Where to write *_longest.tsv (default: 07_NRdatabase/00_2_MF)",
    )
    parser.add_argument(
        "--output",
        default="",
        help=(
            "Optional path for the replacement TSV (kept genes). "
            "Default: <stem>_longest.tsv next to the input. "
            "The dropped-gene table is written beside this file."
        ),
    )
    args = parser.parse_args()

    cand_dir = Path(args.out_dir)
    families = ["CYP", "OMT"] if args.family == "both" else [args.family]

    print("=" * 60)
    print("  Collapse to longest isoform (original TSV is not modified)")
    print("=" * 60)
    print(f"  FASTA: {args.fasta}")
    print(f"  Group: {args.group}")
    print()

    def replacement_paths(src):
        """Name the replacement TSV and the dropped-gene sidecar."""
        src = Path(src)
        if args.output:
            kept = Path(args.output)
        else:
            kept = cand_dir / f"{src.stem}_longest.tsv"
        dropped = kept.with_name(kept.stem + "_dropped.tsv")
        return kept, dropped

    if args.candidates:
        # One custom file — ignore --family loop.
        src = Path(args.candidates)
        kept, dropped = replacement_paths(src)
        collapse_table(src, args.fasta, kept, dropped, args.group)
    else:
        for family in families:
            if family == "CYP":
                src = cand_dir / "cyp_candidates_MF.tsv"
            else:
                src = cand_dir / "omt_candidates_MF.tsv"
            print(f"### {family}")
            kept, dropped = replacement_paths(src)
            collapse_table(src, args.fasta, kept, dropped, args.group)
            print()

    print("Done. Original candidate TSVs are unchanged.")
    print("Replacement heatmap file: *_longest.tsv")
    print("Next: sbatch plots using that replacement file.")
    print("=" * 60)


if __name__ == "__main__":
    main()
