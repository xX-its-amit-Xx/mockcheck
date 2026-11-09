#!/bin/bash
# One-command reproduction on a Slurm cluster (tested on Northeastern's Explorer).
#
#   bash run_all.sh           # submit everything: environment + databases, reads, classifiers
#   bash run_all.sh --local   # no cluster needed: re-score the committed raw reports
#
# Environment: MOCKCHECK_SCRATCH (default /scratch/$USER/mockcheck) holds the big files.
# Set SLURM_ACCOUNT if your site needs one. Public/short partition only, modest footprint.
set -euo pipefail
cd "$(dirname "$0")"

if [ "${1:-}" = "--local" ]; then
  python scripts/evaluate.py     # needs: python3, pandas, matplotlib
  python scripts/build_site.py   # refreshes docs/data/results for the dashboard
  echo "Done. Tables are in results/, figures in results/fig/, site in docs/."
  exit 0
fi

mkdir -p logs
ACC=${SLURM_ACCOUNT:+--account=$SLURM_ACCOUNT}
setup=$(sbatch --parsable $ACC scripts/setup_env_and_dbs.sbatch)
fetch=$(sbatch --parsable $ACC scripts/fetch_reads.sbatch)
classify=$(sbatch --parsable $ACC --dependency=afterok:$setup:$fetch scripts/run_classifiers.sbatch)
echo "submitted: setup=$setup fetch=$fetch classify=$classify (classify waits for both)"
echo "when classify finishes:  bash scripts/collect_raw.sh && bash run_all.sh --local"
