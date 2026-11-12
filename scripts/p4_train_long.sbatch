#!/usr/bin/env bash
#SBATCH --job-name=mockcheck-p4long
#SBATCH --account=b.gyori
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/shenoy.am/mockcheck/logs/p4long.%j.out
# Phase 4 budget control (3x steps, 5x encoder lr): 24 loss conditions x 3 seeds (finished runs are skipped) on one GPU. Resumable: finished runs are skipped, so a re-submit continues.
# Public gpu partition, 1 GPU, started in the weekday night window. Submit from the repo root: sbatch scripts/p4_train.sbatch
set -euo pipefail
S=/scratch/shenoy.am/mockcheck
. "$S/env_ft/bin/activate"
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True HF_HOME=$S/hf HF_HUB_OFFLINE=1 TOKENIZERS_PARALLELISM=false
cd "$S/repo"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
python scripts/p4_train.py --data data/groundtruth/p4_groundtruth.tsv.gz --out "$S/out/p4" \
  --base "$S/models/esm2_35m" --conds ce_long --seeds 0,1,2 --steps 200
echo P4_TRAIN_DONE
