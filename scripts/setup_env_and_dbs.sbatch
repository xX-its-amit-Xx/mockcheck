#!/bin/bash
#SBATCH --job-name=mc_setup
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x.%j.out
# Builds the pinned-by-lockfile environment and downloads databases. Idempotent.
set -euo pipefail
REPO="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}"
source "$REPO/scripts/common.sh"
mkdir -p "$MC_DB" "$MOCKCHECK_SCRATCH"
module load miniforge3/26.7.2-0
eval "$(conda shell.bash hook)"
if [ ! -d "$MC_ENV" ]; then
  conda create -y -p "$MC_ENV" -c conda-forge -c bioconda \
    kraken2 bracken metaphlan seqtk pandas matplotlib
fi
conda activate "$MC_ENV"
conda list --export > "$MOCKCHECK_SCRATCH/environment.lock.txt"
# Kraken2 databases (Langmead indexes, 2025-04-02 build)
for name in k2_standard_08gb_20250402 k2_standard_16gb_20250402; do
  d="$MC_DB/$name"
  if [ ! -f "$d/hash.k2d" ]; then
    mkdir -p "$d"
    curl -fL -C - -o "$MC_DB/$name.tar.gz" "https://genome-idx.s3.amazonaws.com/kraken/$name.tar.gz"
    tar -xzf "$MC_DB/$name.tar.gz" -C "$d"
    rm -f "$MC_DB/$name.tar.gz"
  fi
done
# MetaPhlAn database
if [ ! -d "$MC_DB/metaphlan" ]; then
  mkdir -p "$MC_DB/metaphlan"
  metaphlan --install --bowtie2db "$MC_DB/metaphlan"
fi
ls -la "$MC_DB" "$MC_DB"/metaphlan | head -40
echo SETUP_DONE
