#!/bin/bash
# Stream the first N read pairs of each run in config/samples.tsv from ENA.
# Never downloads a full FASTQ. Output: $MC_READS/<sample>_{1,2}.fastq.gz
set -uo pipefail
here="$(cd "$(dirname "$0")/.." && pwd)"
source "$here/scripts/common.sh"
mkdir -p "$MC_READS"
tail -n +2 "$here/config/samples.tsv" | while IFS=$'\t' read -r sample acc role n truth; do
  if [ -s "$MC_READS/${sample}_2.fastq.gz" ]; then echo "skip $sample"; continue; fi
  urls=$(curl -fsS "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$acc&result=read_run&fields=fastq_ftp&format=tsv" | tail -n 1 | cut -f2)
  IFS=';' read -r u1 u2 <<< "$urls"
  i=1
  for u in "$u1" "$u2"; do
    # head closes the pipe early on purpose; the pipeline status is not checked
    curl -fsS "https://$u" | gzip -dc | head -n $((4*n)) | gzip -1 > "$MC_READS/${sample}_${i}.fastq.gz.tmp"
    mv "$MC_READS/${sample}_${i}.fastq.gz.tmp" "$MC_READS/${sample}_${i}.fastq.gz"
    i=$((i+1))
  done
  echo "fetched $sample ($acc) $(zcat "$MC_READS/${sample}_1.fastq.gz" | wc -l | awk '{print $1/4}') pairs"
done
