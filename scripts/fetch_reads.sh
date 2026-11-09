#!/bin/bash
# Stream the first N read pairs of each run in config/samples.tsv from ENA.
# Never downloads a full FASTQ. Output: $MC_READS/<sample>_{1,2}.fastq.gz
set -uo pipefail
here="$(cd "$(dirname "$0")/.." && pwd)"
source "$here/scripts/common.sh"
mkdir -p "$MC_READS"
tail -n +2 "$here/config/samples.tsv" | while IFS=$'\t' read -r sample acc role n truth; do
  if [ -s "$MC_READS/${sample}_2.fastq.gz" ] && [ "$(zcat "$MC_READS/${sample}_2.fastq.gz" | head -c 1 | wc -c)" -gt 0 ]; then echo "skip $sample"; continue; fi
  urls=$(curl -fsS "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$acc&result=read_run&fields=fastq_ftp&format=tsv" | tail -n 1 | cut -f2)
  IFS=';' read -r u1 u2 <<< "$urls"
  ok=0
  for attempt in 1 2 3; do
    i=1
    for u in "$u1" "$u2"; do
      # head closes the pipe early on purpose, so the pipeline status is not checked;
      # correctness is verified below by counting lines in BOTH mates.
      curl -fsS "https://$u" | gzip -dc | head -n $((4*n)) | gzip -1 > "$MC_READS/${sample}_${i}.fastq.gz.tmp"
      mv "$MC_READS/${sample}_${i}.fastq.gz.tmp" "$MC_READS/${sample}_${i}.fastq.gz"
      i=$((i+1))
    done
    l1=$(zcat "$MC_READS/${sample}_1.fastq.gz" | wc -l); l2=$(zcat "$MC_READS/${sample}_2.fastq.gz" | wc -l)
    if [ "$l1" -gt 0 ] && [ "$l1" -eq "$l2" ]; then ok=1; break; fi
    echo "retry $sample attempt $attempt: mate lines $l1 vs $l2" >&2; sleep 5
  done
  if [ "$ok" -ne 1 ]; then echo "FAILED $sample ($acc)" >&2; rm -f "$MC_READS/${sample}_"{1,2}.fastq.gz; continue; fi
  echo "fetched $sample ($acc) $((l1/4)) pairs (both mates verified)"
done
