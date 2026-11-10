#!/usr/bin/env python3
"""Generate COOKBOOK.md from docs/data/cookbook.json (the dashboard reads the same JSON)."""
import json, os
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
d = json.load(open(f"{root}/docs/data/cookbook.json", encoding="utf-8"))
out = ["# Cookbook", "", d["intro"], "", "_Generated from `docs/data/cookbook.json` by `scripts/build_cookbook.py`; the dashboard shows the same recipes._", ""]
for r in d["recipes"]:
    out += [f"## {r['title']}", "", f"**For:** {r['for']}", "", "```bash"] + r["steps"] + ["```", "", f"**You should see:** {r['expect']}", "", f"**If it goes wrong:** {r['if_wrong']}", ""]
open(f"{root}/COOKBOOK.md", "w", encoding="utf-8", newline="\n").write("\n".join(out))
print("wrote COOKBOOK.md")
