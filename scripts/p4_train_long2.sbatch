#!/usr/bin/env bash
#SBATCH --job-name=mockcheck-p4long2
#SBATCH --account=b.gyori
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=03:00:00
#SBATCH --output=/scratch/shenoy.am/mockcheck/logs/p4long2.%j.out
# Phase 4 transfer experiment at the proper training budget (3x steps, 5x encoder lr): warm-started utility conditions.
# Average and worst-stratum: 3 seeds. One-stratum-only (9 strata): 2 seeds. Resumable: finished runs are skipped.
set -euo pipefail
S=/scratch/shenoy.am/mockcheck
. "$S/env_ft/bin/activate"
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True HF_HOME=$S/hf HF_HUB_OFFLINE=1 TOKENIZERS_PARALLELISM=false
cd "$S/repo"
python scripts/p4_train.py --data data/groundtruth/p4_groundtruth.tsv.gz --out "$S/out/p4" --base "$S/models/esm2_35m" --conds wu_avg_long,wu_worst_long --seeds 0,1,2 --steps 200
python scripts/p4_train.py --data data/groundtruth/p4_groundtruth.tsv.gz --out "$S/out/p4" --base "$S/models/esm2_35m" --conds wu_single_long --seeds 0,1 --steps 200
echo P4_LONG2_DONE
