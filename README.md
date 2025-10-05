

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
# test sync
