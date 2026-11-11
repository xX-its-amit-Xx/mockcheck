#!/bin/bash
# Run on the cluster after phase2_classify: shrink Phase 2 outputs into results/raw/phase2 and copy the spike manifest.
set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
source "$REPO/scripts/common.sh"
R="$REPO/results/raw/phase2"; mkdir -p "$R"/{kraken2,metaphlan}
for f in "$MC_OUT"/p2/kraken2/*.report; do awk -F'\t' '$2>0' "$f" | gzip -9 > "$R/kraken2/$(basename "$f").gz"; done
cp "$MC_OUT"/p2/metaphlan/*.profile.txt "$R/metaphlan/" 2>/dev/null || true
sed -i "s|$MOCKCHECK_SCRATCH|\$MOCKCHECK_SCRATCH|g" "$R"/metaphlan/*.profile.txt 2>/dev/null || true
cp "$MC_OUT/p2/timing/timing.tsv" "$R/timing.tsv" 2>/dev/null || true
du -sh "$R"
