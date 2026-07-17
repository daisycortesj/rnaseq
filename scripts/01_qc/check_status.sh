#!/usr/bin/env bash
# check_status.sh — where am I in the MF nutmeg pipeline?
#
# Trusts the filesystem, not the README, not memory.
# Safe to run any time. Reads only; changes nothing.
#
# Usage:  bash check_status.sh
#         bash check_status.sh DC      # other species code

CODE="${1:-MF}"
BASE=/projects/tholl_lab_1/daisy_analysis

# ---- paths (MF defaults; adjust if CODE != MF) -------------------------
TRIN=$BASE/01_processed/00_4_${CODE}_trinity/${CODE}_trinity_pooled.Trinity.fasta
GMAP=$BASE/01_processed/00_4_${CODE}_trinity/${CODE}_trinity_pooled.Trinity.fasta.gene_trans_map
CDHIT=$BASE/01_processed/00_6_cdhit/${CODE}_trinity_cdhit95.fasta
PEP=$BASE/01_processed/00_6_cdhit/${CODE}_trinity_cdhit95.fasta.transdecoder.pep
RSEMMAT=$BASE/01_processed/00_7_RSEM/RSEM.gene.counts.matrix
CNT=$BASE/03_count_tables/00_5_${CODE}_trinity/gene_count_matrix.tsv
META=$BASE/03_count_tables/00_5_${CODE}_trinity/sample_metadata.tsv
DE=$BASE/06_analysis/pydeseq2_${CODE}_step1_unfiltered/pydeseq2_results_UNFILTERED.tsv
CYP=$BASE/06_analysis/hmmer_genefinder_${CODE}/cyp450_hmmer_ids.txt
OMT=$BASE/06_analysis/hmmer_genefinder_${CODE}/omt_hmmer_ids.txt

# ---- helper ------------------------------------------------------------
# prints  DONE/MISSING, size, mtime, and a count of records
chk () {
  local label="$1" path="$2" kind="$3"
  if [ -f "$path" ]; then
    local size mtime n
    size=$(du -h "$path" 2>/dev/null | cut -f1)
    mtime=$(date -r "$path" "+%Y-%m-%d %H:%M" 2>/dev/null)
    case "$kind" in
      fasta) n="$(grep -c '^>' "$path" 2>/dev/null) seqs" ;;
      lines) n="$(wc -l < "$path" 2>/dev/null) lines" ;;
      table) n="$(( $(wc -l < "$path" 2>/dev/null) - 1 )) rows" ;;
      *)     n="" ;;
    esac
    printf "  [DONE]    %-22s %8s  %16s  %s\n" "$label" "$size" "$mtime" "$n"
  else
    printf "  [MISSING] %-22s %s\n" "$label" "$path"
  fi
}

echo "======================================================================"
echo " PIPELINE STATUS — species code: $CODE"
echo " $(date)"
echo "======================================================================"

echo
echo "--- ASSEMBLY -------------------------------------------------------"
chk "Trinity pooled"    "$TRIN"    fasta
chk "gene_trans_map"    "$GMAP"    lines
chk "CD-HIT 95%"        "$CDHIT"   fasta

echo
echo "--- QUANTIFICATION -------------------------------------------------"
chk "RSEM matrix"       "$RSEMMAT" table
chk "count matrix"      "$CNT"     table
chk "sample metadata"   "$META"    lines

echo
echo "--- ANNOTATION -----------------------------------------------------"
chk "TransDecoder .pep" "$PEP"     fasta
chk "HMMER CYP ids"     "$CYP"     lines
chk "HMMER OMT ids"     "$OMT"     lines

echo
echo "--- DIFFERENTIAL EXPRESSION ----------------------------------------"
chk "PyDESeq2 unfiltered" "$DE"    table
echo "  (searching for any other pydeseq2 outputs...)"
find "$BASE/06_analysis" -maxdepth 2 \
     \( -iname "*UNFILTERED*" -o -iname "*FILTERED*" \) 2>/dev/null \
  | sed 's/^/    found: /' | head -10

echo
echo "======================================================================"
echo " THE ID QUESTION — gene-level or isoform-level?"
echo "======================================================================"
if [ -f "$CNT" ]; then
  echo
  echo "  Count matrix, first 2 row IDs:"
  tail -n +2 "$CNT" | cut -f1 | head -2 | sed 's/^/    /'
  echo
  first=$(tail -n +2 "$CNT" | cut -f1 | head -1)
  case "$first" in
    *_i[0-9]*) echo "    ==> ISOFORM-level  (ends in _iN)" ;;
    TRINITY_DN*_c*_g*) echo "    ==> GENE-level    (no _iN)" ;;
    *) echo "    ==> UNRECOGNIZED format — paste this to Claude" ;;
  esac
elif [ -f "$RSEMMAT" ]; then
  echo
  echo "  (no count matrix; showing RSEM matrix instead)"
  tail -n +2 "$RSEMMAT" | cut -f1 | head -2 | sed 's/^/    /'
else
  echo "  (no matrix found yet)"
fi

if [ -f "$CYP" ]; then
  echo
  echo "  HMMER CYP ids, first 2:"
  head -2 "$CYP" | sed 's/^/    /'
fi

echo
echo "======================================================================"
echo " THE OVERLAP — do DE ids and HMMER ids actually join?"
echo "======================================================================"
if [ -f "$DE" ] && [ -f "$CYP" ]; then
  tail -n +2 "$DE" | cut -f1 | sort -u > /tmp/_de_ids.$$
  echo "    CYP in DE table: $(comm -12 /tmp/_de_ids.$$ <(sort -u "$CYP") | wc -l)"
  [ -f "$OMT" ] && echo "    OMT in DE table: $(comm -12 /tmp/_de_ids.$$ <(sort -u "$OMT") | wc -l)"
  rm -f /tmp/_de_ids.$$
else
  echo "    (skipped — need both PyDESeq2 results and HMMER ids)"
fi

echo
echo "======================================================================"
echo " RECENT JOBS"
echo "======================================================================"
echo "  Still running:"
squeue -u "$USER" 2>/dev/null | sed 's/^/    /'
echo
echo "  Last 8 log files touched:"
find "$BASE" -maxdepth 3 -name "*.out" -newermt "-14 days" 2>/dev/null \
  | xargs ls -lt 2>/dev/null | head -8 | awk '{print "    "$6" "$7" "$8"  "$9}'

echo
echo "======================================================================"
echo " Paste this whole output back to Claude."
echo "======================================================================"