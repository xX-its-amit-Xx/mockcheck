#!/bin/bash
#SBATCH --job-name=mc_classify
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=18:00:00
#SBATCH --output=logs/%x.%j.out
# Runs Kraken2 (2 DBs x 2 confidences), Bracken, MetaPhlAn 4 on every sample.
# Resumable: a finished output is skipped. Submit from the repo root.
set -uo pipefail
REPO="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}"
source "$REPO/scripts/common.sh"
mc_activate
T=${SLURM_CPUS_PER_TASK:-8}
mkdir -p "$MC_OUT"/{kraken2,bracken,metaphlan} "$MC_OUT/timing"
timeit() { # label, then command; appends label, wall seconds, max RSS KB to timing.tsv
  local label=$1; shift
  /usr/bin/time -f "$label\t%e\t%M" -a -o "$MC_OUT/timing/timing.tsv" "$@"
}
tail -n +2 "$REPO/config/samples.tsv" | while IFS=$'\t' read -r sample acc role n truth; do
  r1="$MC_READS/${sample}_1.fastq.gz"; r2="$MC_READS/${sample}_2.fastq.gz"
  [ -s "$r2" ] || { echo "missing reads for $sample"; continue; }
  # read length from the first read, used by Bracken
  rl=$(zcat "$r1" | sed -n 2p | awk '{print length($0)}')
  for dbname in k2_standard_08gb_20250402 k2_standard_16gb_20250402; do
    short=${dbname#k2_standard_}; short=${short%%_*}
    for conf in 0.0 0.1; do
      tag="${sample}.${short}.c${conf}"
      rep="$MC_OUT/kraken2/$tag.report"
      if [ ! -s "$rep" ]; then
        timeit "kraken2.${short}.c${conf}.$sample" kraken2 --db "$MC_DB/$dbname" --paired \
          --confidence $conf --threads $T --report "$rep" --output /dev/null "$r1" "$r2"
      fi
      if [ "$short" = "16gb" ] && [ "$conf" = "0.0" ] && [ ! -s "$MC_OUT/bracken/$tag.bracken" ]; then
        # Bracken needs a database100/150mers file; use the nearest available read length
        bl=100; [ "$rl" -ge 125 ] && bl=150
        timeit "bracken.${short}.$sample" bracken -d "$MC_DB/$dbname" -i "$rep" \
          -o "$MC_OUT/bracken/$tag.bracken" -w "$MC_OUT/bracken/$tag.breport" -r $bl -l S -t 0
      fi
    done
  done
  mp="$MC_OUT/metaphlan/$sample.profile.txt"
  if [ ! -s "$mp" ]; then
    timeit "metaphlan.$sample" metaphlan "$r1,$r2" --input_type fastq --nproc $T \
      --bowtie2db "$MC_DB/metaphlan" --bowtie2out "$MC_OUT/metaphlan/$sample.bt2out" \
      --no_map -o "$mp"
  fi
  echo "done $sample"
done
echo CLASSIFY_DONE
