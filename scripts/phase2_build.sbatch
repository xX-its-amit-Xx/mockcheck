#!/bin/bash
#SBATCH --job-name=mc_p2build
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x.%j.out
# Fetch Phase 2 reads, then build the spike-in samples. Submit from the repo root.
set -uo pipefail
REPO="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}"
source "$REPO/scripts/common.sh"
bash "$REPO/scripts/phase2_fetch.sh" || { echo "fetch failed"; exit 1; }
mc_activate
if [ ! -s "$MC_READS/p2/spike_H0_2.fastq.gz" ] || [ ! -s "$MC_READS/p2/spike_P3_2.fastq.gz" ]; then
  python "$REPO/scripts/make_spikein.py" "$MC_READS" "$MC_READS/p2" 4000000 || exit 1
fi
echo BUILD_DONE
