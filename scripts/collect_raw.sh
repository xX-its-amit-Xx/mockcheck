#!/bin/bash
# Run on the cluster after run_classifiers: shrink outputs into results/raw (small, committable).
# Kraken2 reports keep only taxa with clade reads > 0, gzipped.
set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
source "$REPO/scripts/common.sh"
R="$REPO/results/raw"; mkdir -p "$R"/{kraken2,bracken,metaphlan}
for f in "$MC_OUT"/kraken2/*.report; do
  awk -F'\t' '$2>0' "$f" | gzip -9 > "$R/kraken2/$(basename "$f").gz"
done
cp "$MC_OUT"/bracken/*.bracken "$R/bracken/"
cp "$MC_OUT"/metaphlan/*.profile.txt "$R/metaphlan/"
sed -i "s|$MOCKCHECK_SCRATCH|\$MOCKCHECK_SCRATCH|g" "$R"/metaphlan/*.profile.txt   # scrub site-specific paths
cp "$MC_OUT/timing/timing.tsv" "$R/timing.tsv"
cp "$MOCKCHECK_SCRATCH/environment.lock.txt" "$REPO/env/environment.lock.txt" 2>/dev/null || true
du -sh "$R"
