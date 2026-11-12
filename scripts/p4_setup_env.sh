#!/usr/bin/env bash
# Build the Phase 4 training environment on Explorer scratch (login node, CPU only; downloads wheels and model weights).
# Idempotent: skips what exists. Run: nohup bash scripts/p4_setup_env.sh > /scratch/shenoy.am/mockcheck/logs/p4_env.log 2>&1 &
set -euo pipefail
S=/scratch/shenoy.am/mockcheck
module load anaconda3/2024.06 2>/dev/null || true
if [ ! -x "$S/env_ft/bin/python" ]; then
  python -m venv "$S/env_ft"
fi
. "$S/env_ft/bin/activate"
python -c "import torch" 2>/dev/null || pip install --quiet torch --index-url https://download.pytorch.org/whl/cu121
python -c "import transformers, sklearn, numpy" 2>/dev/null || pip install --quiet transformers scikit-learn numpy pandas
export HF_HOME=$S/hf
python - <<'PY'
from transformers import AutoTokenizer, AutoModel
for m in ("facebook/esm2_t12_35M_UR50D",):
    AutoTokenizer.from_pretrained(m); AutoModel.from_pretrained(m); print("cached", m)
import torch; print("torch", torch.__version__, "cuda build", torch.version.cuda)
PY
echo ENV_DONE
