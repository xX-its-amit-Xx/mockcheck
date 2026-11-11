#!/bin/bash
# Phase 2 reads: human host background (NA12878, ERR194147) and real plasma runs (Uganda sepsis, PRJNA898830).
# Streams only the first N pairs, picks mates by _1/_2 suffix, verifies both mates, resumable.
set -uo pipefail
here="$(cd "$(dirname "$0")/.." && pwd)"
source "$here/scripts/common.sh"
P2="$MC_READS/p2"; mkdir -p "$P2"
HUMAN_PAIRS=${HUMAN_PAIRS:-4000000}
UG_PAIRS=${UG_PAIRS:-500000}

stream_pairs() { # accession, n_pairs, out_prefix
  local acc=$1 n=$2 out=$3
  if [ -s "${out}_2.fastq.gz" ] && [ "$(zcat "${out}_1.fastq.gz" | wc -l)" -eq "$(zcat "${out}_2.fastq.gz" | wc -l)" ] && [ "$(zcat "${out}_1.fastq.gz" | head -c 1 | wc -c)" -gt 0 ]; then echo "skip $out"; return 0; fi
  local urls u1 u2
  urls=$(curl -fsS "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$acc&result=read_run&fields=fastq_ftp&format=tsv" | tail -n 1 | cut -f2)
  u1=$(echo "$urls" | tr ';' '\n' | grep '_1\.fastq\.gz$' | head -1); u2=$(echo "$urls" | tr ';' '\n' | grep '_2\.fastq\.gz$' | head -1)
  [ -n "$u1" ] && [ -n "$u2" ] || { echo "FAILED $acc: no _1/_2 files in $urls" >&2; return 1; }
  for attempt in 1 2 3; do
    local i=1
    for u in "$u1" "$u2"; do
      curl -fsS "https://$u" | gzip -dc | head -n $((4*n)) | gzip -1 > "${out}_${i}.fastq.gz.tmp"
      mv "${out}_${i}.fastq.gz.tmp" "${out}_${i}.fastq.gz"; i=$((i+1))
    done
    local l1 l2; l1=$(zcat "${out}_1.fastq.gz" | wc -l); l2=$(zcat "${out}_2.fastq.gz" | wc -l)
    if [ "$l1" -gt 0 ] && [ "$l1" -eq "$l2" ]; then echo "fetched $acc -> $out: $((l1/4)) pairs (both mates verified)"; return 0; fi
    echo "retry $acc attempt $attempt: mate lines $l1 vs $l2" >&2; sleep 5
  done
  rm -f "${out}_"{1,2}.fastq.gz; echo "FAILED $acc" >&2; return 1
}

stream_pairs ERR194147 "$HUMAN_PAIRS" "$P2/human"
tail -n +2 "$here/config/phase2_uganda_runs.tsv" | while IFS=$'\t' read -r run alias total; do
  stream_pairs "$run" "$UG_PAIRS" "$P2/ug_$run"
done
echo FETCH_DONE
