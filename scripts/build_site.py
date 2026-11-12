#!/usr/bin/env python3
"""Copy the small committed result tables into docs/ so the static site (Cloudflare Pages) can show them.
Run after scripts/evaluate.py and scripts/evaluate_phase2.py. Copies every results/*.tsv."""
import glob, os, shutil

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
dst = f"{ROOT}/docs/data/results"
os.makedirs(dst, exist_ok=True)
for src in sorted(glob.glob(f"{ROOT}/results/*.tsv")):
    shutil.copy(src, f"{dst}/{os.path.basename(src)}")
    print("copied", os.path.basename(src))


# agent traces: copy each distilled trace and write an index the dashboard can list
import json, shutil
tdst = f"{ROOT}/docs/data/traces"
shutil.rmtree(tdst, ignore_errors=True)
os.makedirs(tdst)
runs = []
for src in sorted(glob.glob(f"{ROOT}/results/agent_traces/*.json")):
    shutil.copy(src, tdst)
    runs.append(dict(run=os.path.basename(src)[:-5]))
if runs:
    with open(f"{tdst}/index.json", "w", newline="\n") as idxf:
        json.dump(dict(runs=runs), idxf, indent=1)
    print(f"copied {len(runs)} agent traces")
