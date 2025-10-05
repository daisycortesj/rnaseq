
#!/bin/bash
#SBATCH --job-name=star-toy
#SBATCH --partition=normal_q
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq

bash scripts/test_star_toy.sh
