#!/bin/bash
#SBATCH --job-name=mc_fetch
#SBATCH --partition=short
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x.%j.out
# Submit from the repo root: sbatch scripts/fetch_reads.sbatch
bash "$SLURM_SUBMIT_DIR/scripts/fetch_reads.sh"
