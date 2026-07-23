# Known CYP450 and OMT query proteins (for tblastn gene finder)

These FASTA files supply **characterized protein sequences** that BLAST searches
against your Trinity assembly. The files below are already filled with real
UniProt plant sequences (no PLACEHOLDERs).

## Files

| File | Family | Used by | Contents |
|------|--------|---------|----------|
| `known_cyp450.fasta` | Cytochrome P450 | `run_blast_genefinder.sbatch` | CYP71A13, CYP76C1, CYP82G1, CYP81D1, CYP73A27 |
| `known_omt.fasta` | O-methyltransferase | `run_blast_genefinder.sbatch` | OMT1, COMT1, CCoAOMT (parsley + alfalfa) |

## Why these sequences?

Nutmeg has no close annotated genome, so we use conserved plant CYP/OMT proteins
from model / crop species. BLAST will find *similar* nutmeg transcripts.
HMMER (domain search) remains the main safety net for divergent family members;
combine HMMER ∪ BLAST ID lists afterward.

## How to run after this

```bash
cd /projects/tholl_lab_1/daisy_analysis/05_rnaseq-code
# Sync these FASTA files to ARC if you edited them locally
sbatch scripts/07_domains/run_blast_genefinder.sbatch MF
```

## Format

Protein FASTA only (amino acids). Example header:

```
>sp|O49342|C71AD_ARATH Cytochrome P450 71A13 OS=Arabidopsis thaliana GN=CYP71A13
MEMILSISLCLTTLITLLLLRRFLKRTATKVNLPPSPWRLPVIGNLHQLSLHPHRSLRSL
...
```

## Optional: add more queries

You can append extra UniProt plant CYP/OMT sequences to either file if you want
broader coverage (e.g. more CYP71/CYP76 members from mint, basil, or Piper).
