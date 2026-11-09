#!/bin/bash
#SBATCH --job-name=mc_fetch
#SBATCH --partition=short
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=04:00:00
#SBATCH --output=%x.%j.out
# Run from the repo root or via sbatch scripts/fetch_reads.sbatch
cd "${SLURM_SUBMIT_DIR:-.}"
bash "$(dirname "$(readlink -f "$0")")/fetch_reads.sh"
