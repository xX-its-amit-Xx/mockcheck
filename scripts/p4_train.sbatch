#!/usr/bin/env bash
#SBATCH --job-name=mockcheck-p4train
#SBATCH --account=b.gyori
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=03:00:00
#SBATCH --output=/scratch/shenoy.am/mockcheck/logs/p4train.%j.out
# Phase 4 training: 13 loss conditions x 3 seeds on one GPU. Resumable: finished runs are skipped, so a re-submit continues.
# Public gpu partition, 1 GPU, started in the weekday night window. Submit from the repo root: sbatch scripts/p4_train.sbatch
set -euo pipefail
S=/scratch/shenoy.am/mockcheck
. "$S/env_ft/bin/activate"
export HF_HOME=$S/hf HF_HUB_OFFLINE=1 TOKENIZERS_PARALLELISM=false
cd "$S/repo"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
python scripts/p4_train.py --data data/groundtruth/p4_groundtruth.tsv.gz --out "$S/out/p4" \
  --base "$S/models/esm2_35m" --conds probe,ce,util_avg,util_worst,util_single --seeds 0,1,2 --steps 200
echo P4_TRAIN_DONE
