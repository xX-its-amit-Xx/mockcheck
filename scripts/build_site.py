#!/usr/bin/env python3
"""Copy the small committed result tables and figures into docs/ so the static
site (Cloudflare Pages) can show them. Run after scripts/evaluate.py."""
import os, shutil

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
dst = f"{ROOT}/docs/data/results"
os.makedirs(dst, exist_ok=True)
for name in ("even_summary.tsv", "even_per_organism.tsv", "negative_controls.tsv",
             "log_summary.tsv", "timing_even.tsv"):
    src = f"{ROOT}/results/{name}"
    if os.path.exists(src):
        shutil.copy(src, f"{dst}/{name}")
        print("copied", name)
    else:
        print("missing (skipped)", name)
