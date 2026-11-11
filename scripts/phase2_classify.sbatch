#!/bin/bash
#SBATCH --job-name=mc_p2
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x.%j.out
# Phase 2: classify constructed-truth spike-ins and real plasma runs.
# Kraken2 (16 GB DB, confidence 0.0 and 0.1) on every sample; MetaPhlAn on spike-ins and the first 6 plasma runs.
# Resumable; submit from the repo root.
set -uo pipefail
REPO="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}"
source "$REPO/scripts/common.sh"
mc_activate
T=${SLURM_CPUS_PER_TASK:-8}
P2R="$MC_READS/p2"; P2O="$MC_OUT/p2"
mkdir -p "$P2O"/{kraken2,metaphlan,timing}
MPA_INDEX=$(basename "$(ls "$MC_DB"/metaphlan/*.pkl | head -1)" .pkl)
timeit() { local label=$1; shift; /usr/bin/time -f "$label\t%e\t%M" -a -o "$P2O/timing/timing.tsv" "$@"; }
DB="$MC_DB/k2_standard_16gb_20250402"
n_mpa_ug=0
for s in spike_H0 spike_P0 spike_P1 spike_P2 spike_P3 $(tail -n +2 "$REPO/config/phase2_uganda_runs.tsv" | cut -f1 | sed 's/^/ug_/'); do
  r1="$P2R/${s}_1.fastq.gz"; r2="$P2R/${s}_2.fastq.gz"
  [ -s "$r2" ] || { echo "missing reads for $s"; continue; }
  for conf in 0.0 0.1; do
    rep="$P2O/kraken2/$s.16gb.c$conf.report"
    [ -s "$rep" ] || timeit "kraken2.c$conf.$s" kraken2 --db "$DB" --paired --confidence $conf --threads $T --report "$rep" --output /dev/null "$r1" "$r2"
  done
  want_mpa=0
  case "$s" in spike_*) want_mpa=1;; ug_*) n_mpa_ug=$((n_mpa_ug+1)); [ $n_mpa_ug -le 6 ] && want_mpa=1;; esac
  mp="$P2O/metaphlan/$s.profile.txt"
  if [ $want_mpa -eq 1 ] && [ ! -s "$mp" ]; then
    timeit "metaphlan.$s" env PATH="$MC_ENV_MPA/bin:$PATH" metaphlan "$r1,$r2" --input_type fastq --nproc $T \
      --db_dir "$MC_DB/metaphlan" -x "$MPA_INDEX" --mapout "$P2O/metaphlan/$s.mapout.bz2" -o "$mp"
  fi
  echo "done $s"
done
echo P2_CLASSIFY_DONE
