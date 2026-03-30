#!/usr/bin/env python3
"""
Verify pipeline results against original gene lists and previous student data.

Six checks:
  1. Gene list membership — are all result genes in P450_list_RefSeq.txt?
  2. Previous DE overlap — which DE genes overlap with P450_list_RefSeq_log2fold.txt?
  3. Expression comparison — compare counts against P450_expression_refseq.txt
  4. Log2fold comparison — compare direction/magnitude against P450_expression_refseq_logfold2.txt
  5. Upregulated agreement — are previously upregulated genes also up in yours?
  6. Downregulated agreement — are previously downregulated genes also down in yours?

Previous student used R vs L contrast (positive log2FC = up in Root).
Your contrast can differ — pass --contrast-a and --contrast-b to set yours.
  Default:  --contrast-a R --contrast-b L  (positive log2FC = up in Root)
  DG leaf:  --contrast-a L --contrast-b R  (positive log2FC = up in Leaf)

Usage:
  python scripts/verify_genelist.py \
      --results .../cyp_gene_list.tsv \
      --gene-list .../P450_list_RefSeq.txt \
      --contrast-a L --contrast-b R \
      --prev-de .../P450_list_RefSeq_log2fold.txt \
      --prev-expression .../P450_expression_refseq.txt \
      --prev-log2fold .../P450_expression_refseq_logfold2.txt \
      --prev-upregulated .../P450_upregulated.txt \
      --prev-downregulated .../P450_downregulated_list.txt
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


# ═══════════════════════════════════════════════════════════════════════
# LOADERS — flexible parsers for various file formats
# ═══════════════════════════════════════════════════════════════════════

def load_ids_from_list(filepath):
    """Load gene IDs from .txt (one per line), .csv, or .tsv."""
    path_str = str(filepath)

    if path_str.endswith(".txt"):
        with open(filepath) as fh:
            ids = [line.strip() for line in fh
                   if line.strip() and not line.startswith("#")]
        return set(ids)

    sep = "\t" if path_str.endswith((".tsv", ".tab")) else ","
    df = pd.read_csv(filepath, sep=sep)

    if "gene_id" in df.columns:
        return set(df["gene_id"].astype(str))
    if df.columns[0].startswith("LOC") or df.shape[1] == 1:
        return set(df.iloc[:, 0].astype(str))

    raise ValueError(f"Cannot find gene_id column in {filepath}: {list(df.columns)}")


def load_ids_from_results(filepath):
    """Load gene IDs from a results TSV (gene_id as index or first column)."""
    df = pd.read_csv(filepath, sep="\t", comment="#")

    if "gene_id" in df.columns:
        return set(df["gene_id"].astype(str))
    if df.columns[0] == "Unnamed: 0" or df.iloc[:, 0].astype(str).str.startswith("LOC").any():
        return set(df.iloc[:, 0].astype(str))

    df2 = pd.read_csv(filepath, sep="\t", index_col=0, comment="#")
    if df2.index.astype(str).str.startswith("LOC").any():
        return set(df2.index.astype(str))

    raise ValueError(f"Cannot find gene IDs in {filepath}: {list(df.columns)}")


def load_results_df(filepath):
    """Load results as a DataFrame with gene_id index."""
    df = pd.read_csv(filepath, sep="\t", comment="#")

    if "gene_id" not in df.columns:
        if df.columns[0] == "Unnamed: 0" or df.iloc[:, 0].astype(str).str.startswith("LOC").any():
            df = df.rename(columns={df.columns[0]: "gene_id"})
        else:
            df2 = pd.read_csv(filepath, sep="\t", index_col=0, comment="#")
            df2.index.name = "gene_id"
            df = df2.reset_index()

    df["gene_id"] = df["gene_id"].astype(str)
    return df.set_index("gene_id")


def load_prev_data(filepath):
    """Load previous student's data file flexibly.

    Returns a DataFrame indexed by gene_id with whatever columns are present.
    Handles:
      - TSV with unnamed index col (DESeq2 output saved as .txt)
      - Plain text: one gene ID per line
      - CSV / TSV with header
    """
    path_str = str(filepath)

    # Peek at the file to decide parsing strategy
    with open(filepath) as fh:
        first_line = fh.readline()

    has_tabs = "\t" in first_line
    looks_like_header = not first_line.strip().startswith("LOC")

    # If it has tabs and a header, treat it as a TSV (even if .txt extension)
    if has_tabs and looks_like_header:
        df = pd.read_csv(filepath, sep="\t", index_col=0)
        # gene_id may be the index, a column, or both
        if "gene_id" in df.columns:
            if df.index.astype(str).str.startswith("LOC").any():
                df.index = df.index.astype(str)
                df.index.name = "gene_id"
                df = df.drop(columns=["gene_id"], errors="ignore")
            else:
                df["gene_id"] = df["gene_id"].astype(str)
                df = df.set_index("gene_id")
        else:
            df.index = df.index.astype(str)
            df.index.name = "gene_id"
        return df

    # Plain text: one gene ID per line (possibly with space-separated values)
    if path_str.endswith(".txt") and not has_tabs:
        rows = []
        with open(filepath) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                row = {"gene_id": parts[0]}
                for j, val in enumerate(parts[1:], 1):
                    try:
                        row[f"col_{j}"] = float(val)
                    except ValueError:
                        row[f"col_{j}"] = val
                rows.append(row)

        if not rows:
            return pd.DataFrame(columns=["gene_id"]).set_index("gene_id")
        df = pd.DataFrame(rows)
        df["gene_id"] = df["gene_id"].astype(str)
        return df.set_index("gene_id")

    # CSV / TSV fallback
    sep = "\t" if path_str.endswith((".tsv", ".tab")) else ","
    df = pd.read_csv(filepath, sep=sep)
    if "gene_id" not in df.columns:
        if df.columns[0].startswith("LOC") or df.shape[1] == 1:
            df = df.rename(columns={df.columns[0]: "gene_id"})
        elif df.iloc[:, 0].astype(str).str.startswith("LOC").any():
            df = df.rename(columns={df.columns[0]: "gene_id"})
    df["gene_id"] = df["gene_id"].astype(str)
    return df.set_index("gene_id")


def find_column(df, candidates):
    """Find the first matching column name from a list of candidates (case-insensitive)."""
    lower_map = {c.lower(): c for c in df.columns}
    for name in candidates:
        if name.lower() in lower_map:
            return lower_map[name.lower()]
    return None


# ═══════════════════════════════════════════════════════════════════════
# CHECK FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════

def check1_gene_list(results_path, list_path):
    """Check 1: verify all result genes are in the original gene list."""
    results_ids = load_ids_from_results(results_path)
    list_ids = load_ids_from_list(list_path)

    in_both = results_ids & list_ids
    results_only = results_ids - list_ids
    list_only = list_ids - results_ids

    print("=" * 60)
    print("  CHECK 1: Gene List Membership")
    print("=" * 60)
    print(f"  Results file: {results_path}")
    print(f"    Genes: {len(results_ids)}")
    print(f"  Gene list:    {list_path}")
    print(f"    Genes: {len(list_ids)}")
    print()
    print(f"  In BOTH (verified):       {len(in_both)}")
    print(f"  In results ONLY:          {len(results_only)}")
    print(f"  In gene list ONLY:        {len(list_only)}")
    print()

    if results_only:
        print("  *** UNEXPECTED: genes in results but NOT in gene list ***")
        for g in sorted(results_only):
            print(f"    - {g}")
        print()

    if list_only:
        print(f"  Genes in gene list but NOT in results ({len(list_only)}):")
        print("  (filtered out by DE cutoffs, or not in count matrix)")
        sorted_missing = sorted(list_only)
        if len(sorted_missing) <= 30:
            for g in sorted_missing:
                print(f"    - {g}")
        else:
            for g in sorted_missing[:15]:
                print(f"    - {g}")
            print(f"    ... and {len(sorted_missing) - 15} more")
        print()

    if not results_only:
        pct = len(in_both) / len(list_ids) * 100 if list_ids else 0
        print(f"  PASS: All {len(results_ids)} genes in results are in the gene list")
        print(f"        ({pct:.1f}% of the gene list made it through filtering)")
    else:
        print(f"  FAIL: {len(results_only)} gene(s) in results are NOT in the gene list")
        print("        Something unexpected happened — investigate above genes")
    print()
    return len(results_only) == 0


def check2_prev_de_overlap(results_path, prev_de_path):
    """Check 2: overlap with previous student's DE gene list."""
    print("=" * 60)
    print("  CHECK 2: Previous Student DE Gene Overlap")
    print(f"  File: {prev_de_path}")
    print("=" * 60)

    prev_df = load_prev_data(prev_de_path)
    results_df = load_results_df(results_path)

    prev_ids = set(prev_df.index)
    results_ids = set(results_df.index)
    shared = sorted(prev_ids & results_ids)
    prev_only = sorted(prev_ids - results_ids)
    yours_only = sorted(results_ids - prev_ids)

    print(f"  Previous student DE genes: {len(prev_ids)}")
    print(f"  Your DE genes:             {len(results_ids)}")
    print(f"  In common:                 {len(shared)}")
    print()

    if shared:
        pct = len(shared) / len(prev_ids) * 100
        print(f"  {pct:.0f}% of previous student's DE genes reproduced in your results")
        print()

    if prev_only:
        print(f"  In previous student ONLY ({len(prev_only)}):")
        print("  (may not pass your DE cutoffs, or not in count matrix)")
        for g in prev_only[:20]:
            print(f"    - {g}")
        if len(prev_only) > 20:
            print(f"    ... and {len(prev_only) - 20} more")
        print()

    if yours_only:
        print(f"  In YOUR results ONLY ({len(yours_only)}):")
        print("  (new DE genes not in previous student's list)")
        for g in yours_only[:20]:
            print(f"    - {g}")
        if len(yours_only) > 20:
            print(f"    ... and {len(yours_only) - 20} more")
        print()


def check3_expression(results_path, prev_expr_path):
    """Check 3: compare expression levels with previous student's data."""
    print("=" * 60)
    print("  CHECK 3: Expression Level Comparison")
    print(f"  File: {prev_expr_path}")
    print("=" * 60)

    prev_df = load_prev_data(prev_expr_path)
    results_df = load_results_df(results_path)

    print(f"  Previous student genes: {len(prev_df)}")
    print(f"  Previous student columns: {list(prev_df.columns)}")
    print(f"  Your genes:             {len(results_df)}")
    print()

    prev_ids = set(prev_df.index)
    results_ids = set(results_df.index)
    shared = sorted(prev_ids & results_ids)

    print(f"  Genes in common: {len(shared)}")
    print()

    if not shared:
        print("  No overlapping genes to compare.")
        print()
        return

    # Look for expression/count columns in previous data
    expr_col = find_column(prev_df, [
        "baseMean", "mean_counts", "total_counts", "expression",
        "counts", "RPKM", "FPKM", "TPM", "count", "mean",
    ])
    your_expr_col = find_column(results_df, [
        "baseMean", "mean_counts", "total_counts",
    ])

    if expr_col and your_expr_col:
        prev_vals = prev_df.loc[shared, expr_col].astype(float)
        your_vals = results_df.loc[shared, your_expr_col].astype(float)

        valid = prev_vals.notna() & your_vals.notna()
        if valid.sum() >= 3:
            corr = prev_vals[valid].corr(your_vals[valid])
            print(f"  Expression correlation (Pearson): {corr:.3f}")
            print(f"    Compared: prev '{expr_col}' vs yours '{your_expr_col}'")
            print(f"    Genes with data: {valid.sum()}")
            if corr > 0.7:
                print("    Strong positive correlation — expression patterns agree")
            elif corr > 0.3:
                print("    Moderate correlation")
            else:
                print("    Weak/no correlation — different expression patterns")
            print()

    # Compare which genes are highly expressed in both
    if your_expr_col:
        your_top = results_df.loc[shared].nlargest(10, your_expr_col)
        print(f"  Top 10 shared genes by your {your_expr_col}:")
        for gid in your_top.index:
            your_val = your_top.at[gid, your_expr_col]
            prev_val = ""
            if expr_col and gid in prev_df.index and pd.notna(prev_df.at[gid, expr_col]):
                prev_val = f"  prev {expr_col}={prev_df.at[gid, expr_col]:,.1f}"
            print(f"    {gid}  yours={your_val:,.1f}{prev_val}")
        print()

    # If previous data has sample columns, show overlap summary
    numeric_cols = prev_df.select_dtypes(include=[np.number]).columns.tolist()
    if len(numeric_cols) > 2:
        print(f"  Previous data has {len(numeric_cols)} numeric columns:")
        for col in numeric_cols[:10]:
            mean_val = prev_df[col].mean()
            print(f"    {col}: mean={mean_val:,.1f}")
        if len(numeric_cols) > 10:
            print(f"    ... and {len(numeric_cols) - 10} more")
        print()


    # CHANGED: added contrast_a, contrast_b parameters (was hardcoded R vs L only).
    # Previously positive log2FC was always interpreted as "up_root" regardless of
    # which contrast the user ran. Now flips interpretation when contrast_a="L".
def check4_log2fold(results_path, prev_lfc_path, contrast_a="R", contrast_b="L"):
    """Check 4: compare log2FoldChange direction and magnitude.

    Previous student always used R vs L (positive = up_root).
    Your contrast is set by contrast_a / contrast_b (positive = up in A).
    """
    print("=" * 60)
    print("  CHECK 4: Log2FoldChange Direction & Magnitude")
    print(f"  File: {prev_lfc_path}")
    label_a = {"R": "root", "L": "leaf"}.get(contrast_a, contrast_a)
    label_b = {"R": "root", "L": "leaf"}.get(contrast_b, contrast_b)
    print(f"  Your contrast: {contrast_a} vs {contrast_b} "
          f"(positive log2FC = up in {label_a})")
    print(f"  Previous student contrast: R vs L (positive log2FC = up in root)")
    print("=" * 60)

    prev_df = load_prev_data(prev_lfc_path)
    results_df = load_results_df(results_path)

    print(f"  Previous student genes: {len(prev_df)}")
    print(f"  Previous student columns: {list(prev_df.columns)}")
    print(f"  Your genes:             {len(results_df)}")
    print()

    prev_ids = set(prev_df.index)
    results_ids = set(results_df.index)
    shared = sorted(prev_ids & results_ids)

    print(f"  Genes in common: {len(shared)}")
    print()

    if not shared:
        print("  No overlapping genes to compare.")
        print()
        return

    # Find log2FC column in previous data
    prev_lfc_col = find_column(prev_df, [
        "log2FoldChange", "log2fold", "log2fc", "lfc", "FC",
        "log2_fold_change", "foldchange", "fold_change",
    ])
    your_lfc_col = find_column(results_df, [
        "log2FoldChange", "log2fold", "log2fc",
    ])

    # If prev file has only 1 numeric column, assume it's the log2FC
    if prev_lfc_col is None:
        numeric_cols = prev_df.select_dtypes(include=[np.number]).columns.tolist()
        if len(numeric_cols) == 1:
            prev_lfc_col = numeric_cols[0]
            print(f"  (Assuming '{prev_lfc_col}' is the log2FC column)")
            print()

    if not prev_lfc_col:
        print("  Cannot find log2FoldChange column in previous data.")
        print(f"  Available columns: {list(prev_df.columns)}")
        print()
        _report_overlap(shared, prev_ids, results_ids)
        return

    if not your_lfc_col:
        print("  Cannot find log2FoldChange column in your results.")
        print(f"  Available columns: {list(results_df.columns)}")
        print()
        _report_overlap(shared, prev_ids, results_ids)
        return

    # Previous student: always R vs L (positive = up_root)
    # Your contrast: positive = up in contrast_a
    # To compare tissue direction, convert both to up_root / up_leaf
    # CHANGED: was always your_dir = "up_root" if positive. Now checks if
    # your contrast is flipped (L vs R) and reverses the mapping so that
    # e.g. your positive log2FC (= up in leaf) correctly maps to "up_leaf".
    your_contrast_flipped = (contrast_a == "L")

    # ── Direction comparison ──
    agree = 0
    disagree = 0
    no_data = 0
    disagree_genes = []
    agree_genes = []

    for g in shared:
        prev_lfc = prev_df.at[g, prev_lfc_col] if pd.notna(prev_df.at[g, prev_lfc_col]) else None
        your_lfc = results_df.at[g, your_lfc_col] if pd.notna(results_df.at[g, your_lfc_col]) else None

        if prev_lfc is None or your_lfc is None:
            no_data += 1
            continue

        prev_lfc = float(prev_lfc)
        your_lfc = float(your_lfc)
        prev_dir = "up_root" if prev_lfc > 0 else "up_leaf"
        if your_contrast_flipped:
            your_dir = "up_leaf" if your_lfc > 0 else "up_root"
        else:
            your_dir = "up_root" if your_lfc > 0 else "up_leaf"

        if prev_dir == your_dir:
            agree += 1
            agree_genes.append((g, prev_lfc, your_lfc, your_dir))
        else:
            disagree += 1
            disagree_genes.append((g, prev_lfc, your_lfc, prev_dir, your_dir))

    total = agree + disagree
    pct = agree / total * 100 if total else 0

    print(f"  Direction comparison ({total} genes with log2FC in both):")
    print(f"    AGREE  (same direction):   {agree:>4}  ({pct:.0f}%)")
    print(f"    DISAGREE (opposite):       {disagree:>4}")
    if no_data:
        print(f"    Missing data:              {no_data:>4}")
    print()

    if pct >= 80:
        print(f"  Strong agreement — {pct:.0f}% of shared genes have same direction")
    elif pct >= 50:
        print(f"  Moderate agreement — {pct:.0f}% match, review disagreements")
    else:
        print(f"  Low agreement — only {pct:.0f}% match, check experimental conditions")
    print()

    # ── Magnitude correlation ──
    prev_vals = []
    your_vals = []
    for g in shared:
        pv = prev_df.at[g, prev_lfc_col]
        yv = results_df.at[g, your_lfc_col]
        if pd.notna(pv) and pd.notna(yv):
            prev_vals.append(float(pv))
            your_vals.append(float(yv))

    if len(prev_vals) >= 3:
        prev_s = pd.Series(prev_vals)
        your_s = pd.Series(your_vals)
        corr = prev_s.corr(your_s)
        print(f"  Log2FC magnitude correlation (Pearson): {corr:.3f}")
        if corr > 0.7:
            print("    Strong — similar fold-change magnitudes")
        elif corr > 0.3:
            print("    Moderate — some consistency in magnitudes")
        else:
            print("    Weak — fold-change magnitudes differ")
        print()

    # ── Show disagreeing genes ──
    if disagree_genes:
        print(f"  Genes with OPPOSITE direction ({len(disagree_genes)}):")
        for g, plfc, ylfc, pdir, ydir in disagree_genes[:20]:
            print(f"    {g}  prev={plfc:+.2f} ({pdir})  yours={ylfc:+.2f} ({ydir})")
        if len(disagree_genes) > 20:
            print(f"    ... and {len(disagree_genes) - 20} more")
        print()

    # ── Show top agreeing genes by magnitude ──
    if agree_genes:
        top_agree = sorted(agree_genes, key=lambda x: abs(x[2]), reverse=True)[:10]
        print(f"  Top 10 AGREEING genes (by your |log2FC|):")
        for g, plfc, ylfc, direction in top_agree:
            print(f"    {g}  prev={plfc:+.2f}  yours={ylfc:+.2f}  ({direction})")
        print()

    # Look for padj in both datasets
    prev_padj = find_column(prev_df, ["padj", "p_adj", "adjusted_pvalue", "fdr", "qvalue"])
    your_padj = find_column(results_df, ["padj"])

    if prev_padj and your_padj:
        both_sig = 0
        prev_sig_only = 0
        your_sig_only = 0
        for g in shared:
            pp = prev_df.at[g, prev_padj] if pd.notna(prev_df.at[g, prev_padj]) else 1.0
            yp = results_df.at[g, your_padj] if pd.notna(results_df.at[g, your_padj]) else 1.0
            p_sig = float(pp) < 0.05
            y_sig = float(yp) < 0.05
            if p_sig and y_sig:
                both_sig += 1
            elif p_sig:
                prev_sig_only += 1
            elif y_sig:
                your_sig_only += 1
        print(f"  Significance comparison (padj < 0.05):")
        print(f"    Both significant:          {both_sig}")
        print(f"    Previous only:             {prev_sig_only}")
        print(f"    Yours only:                {your_sig_only}")
        print()


    # CHANGED: added contrast_a, contrast_b parameters (was hardcoded R vs L only).
    # Previously assumed positive log2FC = "upregulated in root". Now if contrast_a="L",
    # a gene up in root has NEGATIVE log2FC, so the check flips accordingly.
def check5_upregulated(results_path, upreg_path, contrast_a="R", contrast_b="L"):
    """Check 5: are previously reported upregulated genes also upregulated in your data?

    Previous student: "upregulated" = positive log2FC in R vs L = higher in Root.
    If your contrast is L vs R, a gene that's higher in Root has NEGATIVE log2FC.
    """
    print("=" * 60)
    print("  CHECK 5: Upregulated Gene Agreement")
    print(f"  File: {upreg_path}")
    print(f"  Previous student: upregulated = up in root (R vs L, positive log2FC)")
    label_a = {"R": "root", "L": "leaf"}.get(contrast_a, contrast_a)
    print(f"  Your contrast:    {contrast_a} vs {contrast_b} "
          f"(positive log2FC = up in {label_a})")
    print("=" * 60)

    # load_ids_from_list handles .txt, .csv, and .tsv formats automatically
    prev_upreg = load_ids_from_list(upreg_path)

    # Plain-text lists may use bare numeric IDs (e.g. "108196352" instead of
    # "LOC108196352"), so normalise any that are missing the LOC prefix.
    prev_upreg = {
        f"LOC{gid}" if gid.isdigit() else gid
        for gid in prev_upreg
    }

    results_df = load_results_df(results_path)
    results_ids = set(results_df.index)
    shared = sorted(prev_upreg & results_ids)
    prev_only = sorted(prev_upreg - results_ids)

    print(f"  Previous upregulated genes: {len(prev_upreg)}")
    print(f"  Your DE genes:              {len(results_ids)}")
    print(f"  In common:                  {len(shared)}")
    print()

    if not shared:
        print("  No overlapping genes to compare.")
        if prev_only:
            print(f"  All {len(prev_only)} previous upregulated genes are missing from your results")
            print("  (may not pass your DE cutoffs)")
        print()
        return

    lfc_col = find_column(results_df, ["log2FoldChange", "log2fold", "log2fc"])
    if not lfc_col:
        print("  Cannot find log2FoldChange in your results.")
        print()
        return

    # Previous "upregulated" = up in root = positive in R vs L.
    # In your data, "up in root" means:
    #   contrast R vs L → positive log2FC
    #   contrast L vs R → negative log2FC (flipped)
    your_contrast_flipped = (contrast_a == "L")

    also_up = []
    flipped_down = []
    no_data = 0

    for g in shared:
        lfc = results_df.at[g, lfc_col]
        if pd.isna(lfc):
            no_data += 1
            continue
        lfc = float(lfc)
        # "up in root" in your data: positive if R vs L, negative if L vs R
        is_up_root = (lfc < 0) if your_contrast_flipped else (lfc > 0)
        if is_up_root:
            also_up.append((g, lfc))
        else:
            flipped_down.append((g, lfc))

    total = len(also_up) + len(flipped_down)
    pct = len(also_up) / total * 100 if total else 0

    print(f"  Direction check ({total} genes with log2FC):")
    print(f"    CONFIRMED upregulated (root):  {len(also_up):>4}  ({pct:.0f}%)")
    print(f"    FLIPPED to downregulated:      {len(flipped_down):>4}")
    if no_data:
        print(f"    Missing log2FC:                {no_data:>4}")
    if your_contrast_flipped:
        print(f"    (your contrast is L vs R, so 'up in root' = negative log2FC)")
    print()

    if pct >= 80:
        print(f"  Strong agreement — {pct:.0f}% confirmed as upregulated in root")
    elif pct >= 50:
        print(f"  Moderate — {pct:.0f}% confirmed, review flipped genes")
    else:
        print(f"  Low agreement — only {pct:.0f}% confirmed, check contrast direction")
    print()

    if also_up:
        top = sorted(also_up, key=lambda x: abs(x[1]), reverse=True)[:10]
        print(f"  Top confirmed upregulated in root (by |log2FC|):")
        for g, lfc in top:
            print(f"    {g}  log2FC={lfc:+.2f}")
        print()

    if flipped_down:
        print(f"  FLIPPED genes (prev=up in root, yours=up in leaf):")
        for g, lfc in sorted(flipped_down, key=lambda x: abs(x[1]), reverse=True):
            print(f"    {g}  log2FC={lfc:+.2f}")
        print()

    if prev_only:
        print(f"  Previous upregulated NOT in your results ({len(prev_only)}):")
        print("  (may not pass your DE cutoffs or not in count matrix)")
        for g in prev_only[:20]:
            print(f"    - {g}")
        if len(prev_only) > 20:
            print(f"    ... and {len(prev_only) - 20} more")
        print()


    # CHANGED: added contrast_a, contrast_b parameters (was hardcoded R vs L only).
    # Previously assumed negative log2FC = "up in leaf". Now if contrast_a="L",
    # a gene up in leaf has POSITIVE log2FC, so the check flips accordingly.
def check6_downregulated(results_path, downreg_path, contrast_a="R", contrast_b="L"):
    """Check 6: are previously reported downregulated genes also downregulated in your data?

    Previous student: "downregulated" = negative log2FC in R vs L = higher in Leaf.
    If your contrast is L vs R, higher-in-leaf genes have POSITIVE log2FC.
    """
    print("=" * 60)
    print("  CHECK 6: Downregulated Gene Agreement")
    print(f"  File: {downreg_path}")
    print(f"  Previous student: downregulated = up in leaf (R vs L, negative log2FC)")
    label_a = {"R": "root", "L": "leaf"}.get(contrast_a, contrast_a)
    print(f"  Your contrast:    {contrast_a} vs {contrast_b} "
          f"(positive log2FC = up in {label_a})")
    print("=" * 60)

    # load_ids_from_list handles .txt, .csv, and .tsv formats automatically
    prev_downreg = load_ids_from_list(downreg_path)

    # Plain-text lists may use bare numeric IDs — add LOC prefix if needed
    prev_downreg = {
        f"LOC{gid}" if gid.isdigit() else gid
        for gid in prev_downreg
    }

    results_df = load_results_df(results_path)
    results_ids = set(results_df.index)
    shared = sorted(prev_downreg & results_ids)
    prev_only = sorted(prev_downreg - results_ids)

    print(f"  Previous downregulated genes: {len(prev_downreg)}")
    print(f"  Your DE genes:                {len(results_ids)}")
    print(f"  In common:                    {len(shared)}")
    print()

    if not shared:
        print("  No overlapping genes to compare.")
        if prev_only:
            print(f"  All {len(prev_only)} previous downregulated genes are missing from your results")
            print("  (may not pass your DE cutoffs)")
        print()
        return

    lfc_col = find_column(results_df, ["log2FoldChange", "log2fold", "log2fc"])
    if not lfc_col:
        print("  Cannot find log2FoldChange in your results.")
        print()
        return

    # Previous "downregulated" = up in leaf = negative in R vs L.
    # In your data, "up in leaf" means:
    #   contrast R vs L → negative log2FC
    #   contrast L vs R → positive log2FC (flipped)
    your_contrast_flipped = (contrast_a == "L")

    also_down = []
    flipped_up = []
    no_data = 0

    for g in shared:
        lfc = results_df.at[g, lfc_col]
        if pd.isna(lfc):
            no_data += 1
            continue
        lfc = float(lfc)
        # "up in leaf" in your data: negative if R vs L, positive if L vs R
        is_up_leaf = (lfc > 0) if your_contrast_flipped else (lfc < 0)
        if is_up_leaf:
            also_down.append((g, lfc))
        else:
            flipped_up.append((g, lfc))

    total = len(also_down) + len(flipped_up)
    pct = len(also_down) / total * 100 if total else 0

    print(f"  Direction check ({total} genes with log2FC):")
    print(f"    CONFIRMED up in leaf (downreg in root): {len(also_down):>4}  ({pct:.0f}%)")
    print(f"    FLIPPED to up in root:                  {len(flipped_up):>4}")
    if no_data:
        print(f"    Missing log2FC:                         {no_data:>4}")
    if your_contrast_flipped:
        print(f"    (your contrast is L vs R, so 'up in leaf' = positive log2FC)")
    print()

    if pct >= 80:
        print(f"  Strong agreement — {pct:.0f}% confirmed as up in leaf")
    elif pct >= 50:
        print(f"  Moderate — {pct:.0f}% confirmed, review flipped genes")
    else:
        print(f"  Low agreement — only {pct:.0f}% confirmed, check contrast direction")
    print()

    if also_down:
        top = sorted(also_down, key=lambda x: abs(x[1]), reverse=True)[:10]
        print(f"  Top confirmed up in leaf (by |log2FC|):")
        for g, lfc in top:
            print(f"    {g}  log2FC={lfc:+.2f}")
        print()

    if flipped_up:
        print(f"  FLIPPED genes (prev=up in leaf, yours=up in root):")
        for g, lfc in sorted(flipped_up, key=lambda x: abs(x[1]), reverse=True):
            print(f"    {g}  log2FC={lfc:+.2f}")
        print()

    if prev_only:
        print(f"  Previous downregulated NOT in your results ({len(prev_only)}):")
        print("  (may not pass your DE cutoffs or not in count matrix)")
        for g in prev_only[:20]:
            print(f"    - {g}")
        if len(prev_only) > 20:
            print(f"    ... and {len(prev_only) - 20} more")
        print()


def _report_overlap(shared, prev_ids, results_ids):
    """Fallback: just report gene overlap when columns can't be compared."""
    prev_only = sorted(prev_ids - results_ids)
    yours_only = sorted(results_ids - prev_ids)

    if shared:
        print(f"  Shared genes ({len(shared)}):")
        for g in shared[:30]:
            print(f"    - {g}")
        if len(shared) > 30:
            print(f"    ... and {len(shared) - 30} more")
        print()

    if prev_only:
        print(f"  In previous ONLY ({len(prev_only)}):")
        for g in prev_only[:15]:
            print(f"    - {g}")
        if len(prev_only) > 15:
            print(f"    ... and {len(prev_only) - 15} more")
        print()

    if yours_only:
        print(f"  In yours ONLY ({len(yours_only)}):")
        for g in yours_only[:15]:
            print(f"    - {g}")
        if len(yours_only) > 15:
            print(f"    ... and {len(yours_only) - 15} more")
        print()


# ═══════════════════════════════════════════════════════════════════════
# OUTPUT TABLE — combined comparison file
# ═══════════════════════════════════════════════════════════════════════

    # CHANGED: added contrast_a, contrast_b parameters so the "your_direction"
    # column in the output table correctly says "up_root" or "up_leaf" based on
    # your actual contrast, not always assuming R vs L.
def build_comparison_table(results_path, list_path,
                           prev_de_path=None, prev_expr_path=None,
                           prev_lfc_path=None, prev_upreg_path=None,
                           prev_downreg_path=None,
                           contrast_a="R", contrast_b="L"):
    """Build a single DataFrame with all genes and side-by-side comparison."""

    results_df = load_results_df(results_path)
    gene_list_ids = load_ids_from_list(list_path)

    # Start with all genes from the gene list (union with results)
    all_ids = sorted(gene_list_ids | set(results_df.index))

    your_contrast_flipped = (contrast_a == "L")

    rows = []
    for gid in all_ids:
        row = {"gene_id": gid}

        # ── Your results ──
        in_results = gid in results_df.index
        row["in_your_results"] = "yes" if in_results else "no"
        row["in_gene_list"] = "yes" if gid in gene_list_ids else "no"

        if in_results:
            for col in ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]:
                if col in results_df.columns:
                    val = results_df.at[gid, col]
                    row[f"your_{col}"] = val if pd.notna(val) else ""

            your_lfc = results_df.at[gid, "log2FoldChange"] if "log2FoldChange" in results_df.columns else None
            if pd.notna(your_lfc):
                lfc_val = float(your_lfc)
                # CHANGED: was always "up_root" if positive. Now flips when
                # contrast_a="L" so the direction column matches biology.
                if your_contrast_flipped:
                    row["your_direction"] = "up_leaf" if lfc_val > 0 else "up_root"
                else:
                    row["your_direction"] = "up_root" if lfc_val > 0 else "up_leaf"
            else:
                row["your_direction"] = ""

            if "mean_counts" in results_df.columns:
                row["your_mean_counts"] = results_df.at[gid, "mean_counts"]
            if "total_counts" in results_df.columns:
                row["your_total_counts"] = results_df.at[gid, "total_counts"]

        rows.append(row)

    table = pd.DataFrame(rows)
    table = table.set_index("gene_id")

    # ── Previous student DE list (Check 2) ──
    if prev_de_path and Path(prev_de_path).exists():
        prev_de_df = load_prev_data(prev_de_path)
        table["in_prev_DE_list"] = table.index.map(
            lambda g: "yes" if g in prev_de_df.index else "no"
        )

    # ── Previous student expression (Check 3) ──
    if prev_expr_path and Path(prev_expr_path).exists():
        prev_expr_df = load_prev_data(prev_expr_path)
        table["in_prev_expression"] = table.index.map(
            lambda g: "yes" if g in prev_expr_df.index else "no"
        )
        for col in ["baseMean", "log2FoldChange", "padj"]:
            if col in prev_expr_df.columns:
                table[f"prev_expr_{col}"] = table.index.map(
                    lambda g, c=col: prev_expr_df.at[g, c]
                    if g in prev_expr_df.index and pd.notna(prev_expr_df.at[g, c]) else ""
                )

    # ── Previous student log2fold (Check 4) ──
    if prev_lfc_path and Path(prev_lfc_path).exists():
        prev_lfc_df = load_prev_data(prev_lfc_path)
        table["in_prev_log2fold"] = table.index.map(
            lambda g: "yes" if g in prev_lfc_df.index else "no"
        )
        for col in ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]:
            if col in prev_lfc_df.columns:
                table[f"prev_lfc_{col}"] = table.index.map(
                    lambda g, c=col: prev_lfc_df.at[g, c]
                    if g in prev_lfc_df.index and pd.notna(prev_lfc_df.at[g, c]) else ""
                )

        lfc_col = find_column(prev_lfc_df, ["log2FoldChange", "log2fold", "log2fc"])
        if lfc_col:
            table["prev_lfc_direction"] = table.index.map(
                lambda g: (
                    "up_root" if g in prev_lfc_df.index
                    and pd.notna(prev_lfc_df.at[g, lfc_col])
                    and float(prev_lfc_df.at[g, lfc_col]) > 0
                    else ("up_leaf" if g in prev_lfc_df.index
                          and pd.notna(prev_lfc_df.at[g, lfc_col])
                          else "")
                )
            )

    # ── Previous upregulated (Check 5) ──
    if prev_upreg_path and Path(prev_upreg_path).exists():
        upreg_ids = load_ids_from_list(prev_upreg_path)
        upreg_ids = {f"LOC{gid}" if gid.isdigit() else gid for gid in upreg_ids}
        table["prev_upregulated"] = table.index.map(
            lambda g: "yes" if g in upreg_ids else "no"
        )

    # ── Previous downregulated (Check 6) ──
    if prev_downreg_path and Path(prev_downreg_path).exists():
        downreg_ids = load_ids_from_list(prev_downreg_path)
        downreg_ids = {f"LOC{gid}" if gid.isdigit() else gid for gid in downreg_ids}
        table["prev_downregulated"] = table.index.map(
            lambda g: "yes" if g in downreg_ids else "no"
        )

    # ── Direction agreement column ──
    if "your_direction" in table.columns:
        prev_dir_col = None
        if "prev_lfc_direction" in table.columns:
            prev_dir_col = "prev_lfc_direction"

        if prev_dir_col:
            def _agree(g):
                yd = table.at[g, "your_direction"]
                pd_ = table.at[g, prev_dir_col]
                if not yd or not pd_:
                    return ""
                return "agree" if yd == pd_ else "DISAGREE"
            table["direction_agreement"] = [_agree(g) for g in table.index]

    # ── Sort: matched genes first, then by your padj ──
    sort_cols = []
    if "in_your_results" in table.columns:
        sort_cols.append("in_your_results")
    if "your_padj" in table.columns:
        table["_sort_padj"] = pd.to_numeric(table["your_padj"], errors="coerce").fillna(1)
        sort_cols.append("_sort_padj")

    if sort_cols:
        table = table.sort_values(sort_cols, ascending=[False] + [True] * (len(sort_cols) - 1))
        table = table.drop(columns=["_sort_padj"], errors="ignore")

    return table


# ═══════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Verify pipeline results with 6 checks + output comparison table",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--results", required=True,
                        help="Results file to verify (TSV from any path or step 3)")
    parser.add_argument("--gene-list", required=True,
                        help="Original gene list (P450_list_RefSeq.txt)")
    parser.add_argument("--prev-de", default=None,
                        help="Check 2: Previous student's DE gene list "
                             "(P450_list_RefSeq_log2fold.txt)")
    parser.add_argument("--prev-expression", default=None,
                        help="Check 3: Previous student's expression data "
                             "(P450_expression_refseq.txt)")
    parser.add_argument("--prev-log2fold", default=None,
                        help="Check 4: Previous student's log2fold data "
                             "(P450_expression_refseq_logfold2.txt)")
    parser.add_argument("--prev-upregulated", default=None,
                        help="Check 5: Previous student's upregulated gene list "
                             "(P450_upregulated.txt)")
    parser.add_argument("--prev-downregulated", default=None,
                        help="Check 6: Previous student's downregulated gene list "
                             "(P450_downregulated_list.txt)")
    # CHANGED: new arguments — tell the verify script which contrast YOUR data used.
    # Previous student always used R vs L. If you used L vs R, pass --contrast-a L --contrast-b R
    # so checks 4/5/6 and the comparison table interpret directions correctly.
    parser.add_argument("--contrast-a", default="R",
                        help="YOUR contrast A (numerator). Default: R. "
                             "Set to L if your pipeline used L vs R.")
    parser.add_argument("--contrast-b", default="L",
                        help="YOUR contrast B (denominator). Default: L. "
                             "Set to R if your pipeline used L vs R.")
    parser.add_argument("-o", "--output", default=None,
                        help="Output comparison TSV (side-by-side table)")
    # CHANGED: added --database-label so the summary file ends with
    # _SUMMARY_<label>.txt (e.g. _SUMMARY_ahmed.txt).  When not provided,
    # the old naming (_SUMMARY.txt) is used for backwards compatibility.
    parser.add_argument("--database-label", default=None,
                        help="Name of the comparison database (e.g. 'sukman' or "
                             "'ahmed').  Used in the summary filename so you can "
                             "tell which database was compared.")
    args = parser.parse_args()

    results_path = Path(args.results)
    list_path = Path(args.gene_list)

    if not results_path.exists():
        print(f"ERROR: results file not found: {results_path}")
        sys.exit(1)
    if not list_path.exists():
        print(f"ERROR: gene list not found: {list_path}")
        sys.exit(1)

    label_a = {"R": "root", "L": "leaf"}.get(args.contrast_a, args.contrast_a)
    label_b = {"R": "root", "L": "leaf"}.get(args.contrast_b, args.contrast_b)
    print(f"  Your contrast: {args.contrast_a} vs {args.contrast_b} "
          f"(positive log2FC = up in {label_a})")
    if args.contrast_a != "R":
        print(f"  NOTE: Previous student used R vs L — direction comparisons "
              f"will account for the flip")
    print()

    # ── Check 1 ──
    passed = check1_gene_list(results_path, list_path)

    # ── Check 2 ──
    if args.prev_de:
        p = Path(args.prev_de)
        if p.exists():
            check2_prev_de_overlap(results_path, p)
        else:
            print(f"  WARNING: {p} not found, skipping Check 2")
            print()

    # ── Check 3 ──
    if args.prev_expression:
        p = Path(args.prev_expression)
        if p.exists():
            check3_expression(results_path, p)
        else:
            print(f"  WARNING: {p} not found, skipping Check 3")
            print()

    # ── Check 4 ──
    if args.prev_log2fold:
        p = Path(args.prev_log2fold)
        if p.exists():
            check4_log2fold(results_path, p,
                            contrast_a=args.contrast_a,
                            contrast_b=args.contrast_b)
        else:
            print(f"  WARNING: {p} not found, skipping Check 4")
            print()

    # ── Check 5 ──
    if args.prev_upregulated:
        p = Path(args.prev_upregulated)
        if p.exists():
            check5_upregulated(results_path, p,
                               contrast_a=args.contrast_a,
                               contrast_b=args.contrast_b)
        else:
            print(f"  WARNING: {p} not found, skipping Check 5")
            print()

    if args.prev_downregulated:
        p = Path(args.prev_downregulated)
        if p.exists():
            check6_downregulated(results_path, p,
                                 contrast_a=args.contrast_a,
                                 contrast_b=args.contrast_b)
        else:
            print(f"  WARNING: {p} not found, skipping Check 6")
            print()

    # ── Build and save comparison table + summary ──
    if args.output:
        print("=" * 60)
        print("  Building comparison table...")
        print("=" * 60)
        table = build_comparison_table(
            results_path, list_path,
            prev_de_path=args.prev_de,
            prev_expr_path=args.prev_expression,
            prev_lfc_path=args.prev_log2fold,
            prev_upreg_path=args.prev_upregulated,
            prev_downreg_path=args.prev_downregulated,
            contrast_a=args.contrast_a,
            contrast_b=args.contrast_b,
        )

        out_path = Path(args.output)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(out_path, sep="\t")

        n_yours = (table["in_your_results"] == "yes").sum()
        n_total = len(table)
        print(f"  Rows: {n_total} genes ({n_yours} in your results)")
        print()
        print(f"  Saved: {out_path}")
        print()

        # ── Write summary text file ──
        # CHANGED: when --database-label is provided (e.g. "ahmed"), the
        # summary file ends with _SUMMARY_ahmed.txt instead of _SUMMARY.txt.
        # This strips the trailing database name from the stem first so it
        # doesn't appear twice (verify_DC_filtered_ahmed → verify_DC_filtered).
        if args.database_label:
            label = args.database_label
            stem = out_path.stem
            if stem.endswith(f"_{label}"):
                stem = stem[: -len(f"_{label}")]
            summary_path = out_path.with_name(
                f"{stem}_SUMMARY_{label}.txt"
            )
        else:
            summary_path = out_path.with_name(
                out_path.stem + "_SUMMARY.txt"
            )
        write_summary(table, summary_path, args)
        print(f"  Saved: {summary_path}")
        print()

    print("=" * 60)
    print("  VERIFICATION COMPLETE")
    print("=" * 60)
    sys.exit(0 if passed else 1)


def write_summary(table, summary_path, args):
    """Write a human-readable summary text file."""
    from datetime import datetime
    lines = []
    w = lines.append

    n_total = len(table)
    in_yours = table["in_your_results"] == "yes"
    n_yours = in_yours.sum()
    in_list = table["in_gene_list"] == "yes"
    n_list = in_list.sum()

    w("=" * 70)
    w("  P450 VERIFICATION SUMMARY")
    w(f"  Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    w("=" * 70)
    w("")

    # ── Overview ──
    w("OVERVIEW")
    w("-" * 70)
    w(f"  Your results file:          {args.results}")
    w(f"  Original gene list:         {args.gene_list}")
    w(f"  Total unique genes:         {n_total}")
    w(f"  Genes in your results:      {n_yours}")
    w(f"  Genes in P450 list:         {n_list}")
    w("")

    # ── Check 1: Gene list membership ──
    both = (in_yours & in_list).sum()
    yours_only = (in_yours & ~in_list).sum()
    list_only = (~in_yours & in_list).sum()

    w("CHECK 1: GENE LIST MEMBERSHIP")
    w("-" * 70)
    w(f"  In both (verified):         {both}")
    w(f"  In your results only:       {yours_only}  {'*** UNEXPECTED ***' if yours_only else ''}")
    w(f"  In gene list only:          {list_only}  (filtered out by DE cutoffs)")
    pct = both / n_list * 100 if n_list else 0
    w(f"  Retention rate:             {pct:.1f}% of gene list passed your filters")
    w(f"  Result:                     {'PASS' if yours_only == 0 else 'FAIL'}")
    w("")

    # ── Check 2: Previous DE overlap ──
    if "in_prev_DE_list" in table.columns:
        prev_de = table["in_prev_DE_list"] == "yes"
        n_prev_de = prev_de.sum()
        overlap = (in_yours & prev_de).sum()
        prev_de_only = (prev_de & ~in_yours).sum()
        yours_new = (in_yours & ~prev_de).sum()

        w("CHECK 2: PREVIOUS STUDENT DE OVERLAP")
        w("-" * 70)
        w(f"  Previous student DE genes:   {n_prev_de}")
        w(f"  Overlap with your results:   {overlap}")
        w(f"  Previous DE only:            {prev_de_only}  (didn't pass your cutoffs)")
        w(f"  Your new DE genes:           {yours_new}  (not in previous student's DE)")
        if n_prev_de:
            w(f"  Reproducibility:             {overlap/n_prev_de*100:.0f}% of previous DE genes reproduced")
        w("")

    # ── Check 3: Expression correlation ──
    if "in_prev_expression" in table.columns:
        prev_expr = table["in_prev_expression"] == "yes"
        n_prev_expr = prev_expr.sum()
        expr_overlap = (in_yours & prev_expr).sum()

        w("CHECK 3: EXPRESSION DATA COMPARISON")
        w("-" * 70)
        w(f"  Previous student total genes: {n_prev_expr}")
        w(f"  Overlap with your DE genes:   {expr_overlap}")

        if "your_baseMean" in table.columns and "prev_expr_baseMean" in table.columns:
            shared_mask = in_yours & prev_expr
            yours_bm = pd.to_numeric(table.loc[shared_mask, "your_baseMean"], errors="coerce")
            prev_bm = pd.to_numeric(table.loc[shared_mask, "prev_expr_baseMean"], errors="coerce")
            valid = yours_bm.notna() & prev_bm.notna()
            if valid.sum() >= 3:
                corr = yours_bm[valid].corr(prev_bm[valid])
                w(f"  baseMean correlation:         {corr:.3f}  ({'strong' if corr > 0.7 else 'moderate' if corr > 0.3 else 'weak'})")
        w("")

    # ── Check 4: Log2fold direction ──
    if "direction_agreement" in table.columns:
        agree = (table["direction_agreement"] == "agree").sum()
        disagree = (table["direction_agreement"] == "DISAGREE").sum()
        total_compared = agree + disagree
        pct_agree = agree / total_compared * 100 if total_compared else 0

        w("CHECK 4: LOG2FOLDCHANGE DIRECTION")
        w("-" * 70)
        w(f"  Genes compared:              {total_compared}")
        w(f"  Same direction (agree):      {agree}  ({pct_agree:.0f}%)")
        w(f"  Opposite direction:          {disagree}")
        if pct_agree >= 80:
            w(f"  Assessment:                  STRONG agreement")
        elif pct_agree >= 50:
            w(f"  Assessment:                  Moderate agreement")
        else:
            w(f"  Assessment:                  Low agreement — check contrast direction")
        w("")

        if "your_log2FoldChange" in table.columns and "prev_lfc_log2FoldChange" in table.columns:
            yours_lfc = pd.to_numeric(table["your_log2FoldChange"], errors="coerce")
            prev_lfc = pd.to_numeric(table["prev_lfc_log2FoldChange"], errors="coerce")
            valid = yours_lfc.notna() & prev_lfc.notna()
            if valid.sum() >= 3:
                corr = yours_lfc[valid].corr(prev_lfc[valid])
                w(f"  log2FC magnitude correlation: {corr:.3f}  ({'strong' if corr > 0.7 else 'moderate' if corr > 0.3 else 'weak'})")
                w("")

        if disagree:
            w("  Genes with OPPOSITE direction:")
            dis_mask = table["direction_agreement"] == "DISAGREE"
            for gid in table[dis_mask].index[:15]:
                ylfc = table.at[gid, "your_log2FoldChange"] if "your_log2FoldChange" in table.columns else ""
                plfc = table.at[gid, "prev_lfc_log2FoldChange"] if "prev_lfc_log2FoldChange" in table.columns else ""
                w(f"    {gid}  yours={ylfc}  prev={plfc}")
            if dis_mask.sum() > 15:
                w(f"    ... and {dis_mask.sum() - 15} more (see comparison TSV)")
            w("")

    # ── Check 5: Upregulated ──
    if "prev_upregulated" in table.columns:
        prev_up = table["prev_upregulated"] == "yes"
        n_prev_up = prev_up.sum()
        up_in_yours = (in_yours & prev_up).sum()

        w("CHECK 5: UPREGULATED GENE AGREEMENT")
        w("-" * 70)
        w(f"  Previous upregulated genes:  {n_prev_up}")
        w(f"  Found in your results:       {up_in_yours}")

        if up_in_yours and "your_log2FoldChange" in table.columns:
            shared_up = in_yours & prev_up
            yours_lfc = pd.to_numeric(table.loc[shared_up, "your_log2FoldChange"], errors="coerce")
            confirmed_up = (yours_lfc > 0).sum()
            flipped_down = (yours_lfc <= 0).sum()
            pct_confirmed = confirmed_up / (confirmed_up + flipped_down) * 100 if (confirmed_up + flipped_down) else 0
            w(f"  Confirmed upregulated:       {confirmed_up}  ({pct_confirmed:.0f}%)")
            w(f"  Flipped to downregulated:    {flipped_down}")
            if pct_confirmed >= 80:
                w(f"  Assessment:                  STRONG upregulation agreement")
            elif pct_confirmed >= 50:
                w(f"  Assessment:                  Moderate — review flipped genes")
            else:
                w(f"  Assessment:                  Low — check contrast direction (R vs L)")
        w("")

    # ── Check 6: Downregulated ──
    down_in_yours = 0
    confirmed_down = 0
    pct_confirmed_dn = 0
    if "prev_downregulated" in table.columns:
        prev_down = table["prev_downregulated"] == "yes"
        n_prev_down = prev_down.sum()
        down_in_yours = (in_yours & prev_down).sum()

        w("CHECK 6: DOWNREGULATED GENE AGREEMENT")
        w("-" * 70)
        w(f"  Previous downregulated genes: {n_prev_down}")
        w(f"  Found in your results:       {down_in_yours}")

        if down_in_yours and "your_log2FoldChange" in table.columns:
            shared_down = in_yours & prev_down
            yours_lfc_dn = pd.to_numeric(table.loc[shared_down, "your_log2FoldChange"], errors="coerce")
            confirmed_down = (yours_lfc_dn < 0).sum()
            flipped_up_dn = (yours_lfc_dn >= 0).sum()
            pct_confirmed_dn = confirmed_down / (confirmed_down + flipped_up_dn) * 100 if (confirmed_down + flipped_up_dn) else 0
            w(f"  Confirmed downregulated:     {confirmed_down}  ({pct_confirmed_dn:.0f}%)")
            w(f"  Flipped to upregulated:      {flipped_up_dn}")
            if pct_confirmed_dn >= 80:
                w(f"  Assessment:                  STRONG downregulation agreement")
            elif pct_confirmed_dn >= 50:
                w(f"  Assessment:                  Moderate — review flipped genes")
            else:
                w(f"  Assessment:                  Low — check contrast direction (R vs L)")
        w("")

    # ── Contrast verification note ──
    # CHANGED: was hardcoded "Both analyses use R vs L". Now reads your actual
    # contrast and warns if it's flipped relative to the previous student.
    contrast_a = getattr(args, 'contrast_a', 'R')
    contrast_b = getattr(args, 'contrast_b', 'L')
    label_a = {"R": "ROOT", "L": "LEAF"}.get(contrast_a, contrast_a)
    label_b = {"R": "ROOT", "L": "LEAF"}.get(contrast_b, contrast_b)

    w("CONTRAST VERIFICATION")
    w("-" * 70)
    w(f"  Previous student: contrast = R vs L (positive log2FC = up in ROOT)")
    w(f"  Your pipeline:    contrast = {contrast_a} vs {contrast_b} "
      f"(positive log2FC = up in {label_a})")
    if contrast_a != "R":
        w(f"  NOTE: Your contrast is FLIPPED relative to previous student!")
        w(f"    Your positive log2FC = up in {label_a}")
        w(f"    Previous student positive log2FC = up in ROOT")
        w(f"    Direction comparisons above account for this flip.")
    else:
        w("  Contrasts match — directions are directly comparable.")
    w("")

    # ── Key genes: top DE in both studies ──
    if "your_padj" in table.columns:
        w("TOP DE GENES IN YOUR RESULTS")
        w("-" * 70)
        yours_padj = pd.to_numeric(table["your_padj"], errors="coerce")
        yours_lfc = pd.to_numeric(table.get("your_log2FoldChange", pd.Series(dtype=float)), errors="coerce")
        top_mask = in_yours & yours_padj.notna()
        top_df = table[top_mask].copy()
        top_df["_padj"] = pd.to_numeric(top_df["your_padj"], errors="coerce")
        top_df = top_df.nsmallest(15, "_padj")

        header = f"  {'gene_id':<18} {'your_log2FC':>12} {'your_padj':>12} {'your_dir':>10}"
        if "prev_lfc_log2FoldChange" in table.columns:
            header += f" {'prev_log2FC':>12} {'agree?':>10}"
        w(header)
        w("  " + "-" * (len(header) - 2))

        for gid in top_df.index:
            ylfc = top_df.at[gid, "your_log2FoldChange"] if "your_log2FoldChange" in top_df.columns else ""
            ypadj = top_df.at[gid, "your_padj"] if "your_padj" in top_df.columns else ""
            ydir = top_df.at[gid, "your_direction"] if "your_direction" in top_df.columns else ""

            try:
                ylfc_str = f"{float(ylfc):+.2f}" if ylfc != "" else ""
                ypadj_str = f"{float(ypadj):.2e}" if ypadj != "" else ""
            except (ValueError, TypeError):
                ylfc_str = str(ylfc)
                ypadj_str = str(ypadj)

            line = f"  {gid:<18} {ylfc_str:>12} {ypadj_str:>12} {ydir:>10}"

            if "prev_lfc_log2FoldChange" in table.columns:
                plfc = top_df.at[gid, "prev_lfc_log2FoldChange"] if "prev_lfc_log2FoldChange" in top_df.columns else ""
                agr = top_df.at[gid, "direction_agreement"] if "direction_agreement" in top_df.columns else ""
                try:
                    plfc_str = f"{float(plfc):+.2f}" if plfc != "" else ""
                except (ValueError, TypeError):
                    plfc_str = str(plfc)
                line += f" {plfc_str:>12} {agr:>10}"

            w(line)
        w("")

    # ── Final verdict ──
    w("=" * 70)
    w("FINAL VERDICT")
    w("=" * 70)
    w(f"  Check 1 (gene list):       {'PASS' if yours_only == 0 else 'FAIL'}")
    if "in_prev_DE_list" in table.columns:
        w(f"  Check 2 (DE overlap):      {overlap}/{n_prev_de} reproduced")
    if "direction_agreement" in table.columns:
        w(f"  Check 4 (direction):       {agree}/{total_compared} agree ({pct_agree:.0f}%)")
    if "prev_upregulated" in table.columns and up_in_yours:
        w(f"  Check 5 (upregulated):     {confirmed_up}/{up_in_yours} confirmed ({pct_confirmed:.0f}%)")
    if "prev_downregulated" in table.columns and down_in_yours:
        w(f"  Check 6 (downregulated):   {confirmed_down}/{down_in_yours} confirmed ({pct_confirmed_dn:.0f}%)")
    if contrast_a == "R":
        w(f"  Contrast match:            CONFIRMED (both R vs L)")
    else:
        w(f"  Your contrast:             {contrast_a} vs {contrast_b} (flipped from prev student's R vs L)")
        w(f"                             Direction comparisons above account for the flip")
    w("")
    w("FILES PRODUCED")
    w(f"  Comparison table:  {args.output}")
    w(f"  This summary:      {summary_path}")
    w("=" * 70)

    with open(summary_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")

    for line in lines:
        print(line)


if __name__ == "__main__":
    main()
