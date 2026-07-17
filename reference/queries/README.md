# Known CYP450 and OMT query proteins (for tblastn gene finder)

These FASTA files supply **characterized protein sequences** that BLAST searches
against your Trinity assembly. You need 3–5 sequences per family.

## Where to get sequences

1. Go to [UniProt](https://www.uniprot.org/) or [NCBI Protein](https://www.ncbi.nlm.nih.gov/protein/)
2. Search for characterized CYP450 or OMT from a **related plant** (basal angiosperm,
   Laurales, or well-studied terpenoid pathway species)
3. Download as FASTA and paste into the files below

## Files to edit

| File | Family | Used by |
|------|--------|---------|
| `known_cyp450.fasta` | Cytochrome P450 | `run_blast_genefinder.sbatch` |
| `known_omt.fasta` | O-methyltransferase | `run_blast_genefinder.sbatch` |

## FASTA format example

Each sequence needs a header line starting with `>` followed by amino acids on
the next lines:

```
>sp|P04798|CP51A_HUMAN Cytochrome P450 51A1 OS=Homo sapiens
MELSVLLFLALLTGLLLLLVQAYRFVQKKLQLQKELANTSSKDLTTNHNLQKYGPSYFA
...
```

Use **protein** sequences only (amino acids, not nucleotides). The script runs
`tblastn`, which translates your assembly on the fly and compares it to these
proteins.

## Tips for nutmeg (MF)

- Prefer sequences from plants with similar secondary metabolism (e.g. *Piper*,
  *Ocimum*, *Salvia*, or basal angiosperms if available)
- Include at least one CYP from the CYP71 or CYP76 clades (common in terpenoid biosynthesis)
- Include SAM-dependent methyltransferases for OMT (Pfam Methyltransf_2 family)

Replace the placeholder sequences in `known_cyp450.fasta` and `known_omt.fasta`
before submitting the BLAST job.
