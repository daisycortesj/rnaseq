
# RNA Pipeline (VT ARC – TinkerCliffs)

This is a minimal, beginner-friendly driver that **chooses between STAR genome indexing** (if a reference genome + annotation exist) **or Trinity de novo assembly** (if no suitable reference is found). It’s written to be clear and easy to modify.

## TL;DR quick start

```bash
# 1) Log into ARC (TinkerCliffs) and allocate resources (interactive) or use sbatch.
#    If you use Conda, activate it (adjust path if needed).
source ~/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq

# 2) Basic environment sanity check
bash scripts/inspect_env.sh

# 3) Run in "reference" mode (STAR index) when FASTA + GTF exist
python -m rna_pipeline.cli   --fasta data/reference/genome.fa   --gtf   data/reference/genes.gtf   --outdir results/star_index   --threads 16   --readlen 150

# 4) Run in "de novo" mode (Trinity) when reference is missing
python -m rna_pipeline.cli   --reads-left  data/reads/sample_R1.fq.gz   --reads-right data/reads/sample_R2.fq.gz   --outdir results/trinity   --threads 16   --mem-gb 64
```

## ARC/SLURM example (sbatch)

Submit the included example job (edit partition, mem, time to your needs):

```bash
sbatch scripts/star_index_example.sbatch
```

Inside, it activates your Conda env and calls the CLI with the correct flags. Replace paths as needed. For Trinity, copy the file and swap the command shown below.

## What the pipeline does

1. **Detect inputs**: if `--fasta` and `--gtf` exist and are non-empty → **STAR index**. Otherwise → **Trinity** de novo assembly.
2. **Check tools**: verifies `STAR` (always) and `Trinity` (only if needed) are in your `PATH`.
3. **Run** with clear console logging and fail-fast errors.
4. **Validate** by checking output directories.
5. **Be teachable**: small functions, heavily commented code, and a `--dry` mode so you can preview commands.

## Expected outputs

- **STAR index**: files inside `results/star_index/` (e.g., `Genome`, `SA`, `SAindex`, etc.).
- **Trinity**: directory `results/trinity/` containing Trinity subfolders and the assembled transcripts (`Trinity.fasta`).

## Common ARC notes

- If `conda activate` fails on ARC, run `source ~/miniconda3/etc/profile.d/conda.sh` first.
- Match SLURM resources to tool flags:
  - `--cpus-per-task` ↔ `--threads`
  - `--mem` ↔ `--mem-gb` (Trinity) and STAR memory needs.
- Use a fast filesystem (e.g., scratch) for heavy I/O when possible.

## Modify confidently

Look in `rna_pipeline/`:
- `cli.py`: parses your command-line arguments, calls `main.py`.
- `main.py`: decides STAR vs Trinity based on inputs.
- `tools/star.py` & `tools/trinity.py`: build the exact command lists with safe defaults.
- `utils/`: tiny, readable helpers (`io_utils`, `sys_utils`).
- `runners/local.py`: how commands are executed (you can later add a SLURM runner).

Everything is small and commented so you can tweak paths, threads, or flags without fear.

---
## Quick tests

We include tiny toy inputs so you can test the pipeline safely:

```bash
# DRY RUNS (no heavy compute)
bash scripts/test_dry_runs.sh

# REAL tiny STAR index (runs fast; OK on a compute node)
bash scripts/test_star_toy.sh

# Or submit the tiny STAR test to SLURM
sbatch scripts/star_toy_test.sbatch
```
Toy files:
- `data/reference/toy.fa` + `toy.gtf` for STAR index
- `data/reads/toy_R1.fq.gz` + `toy_R2.fq.gz` for Trinity dry runs
