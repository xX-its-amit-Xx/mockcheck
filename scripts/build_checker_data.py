"""Write docs/data/checker_truth.json (mock truth + aliases) for the in-browser checker, from config/."""
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from mockcheck import truth  # noqa: E402
from mockcheck.parsers import _aliases  # noqa: E402

mocks = {}
for key, title, kind in truth.list_mocks():
    m = truth.load_mock(key)
    mocks[key] = dict(title=m.title, kind=m.kind, dna=m.dna, cell=m.cell,
                      species=sorted(m.species), anchors=m.anchors)
out = dict(aliases=dict(sorted(_aliases.items())), mocks=mocks)
path = os.path.join(ROOT, "docs", "data", "checker_truth.json")
with open(path, "w", newline="\n") as f:
    json.dump(out, f, indent=1, sort_keys=True)
    f.write("\n")
print("wrote", path)
