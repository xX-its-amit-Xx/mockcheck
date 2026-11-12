#!/bin/bash
# One command to verify the repo: unit tests, demo over all committed outputs, site JS syntax, generated files in sync.
set -euo pipefail
cd "$(dirname "$0")/.."
python -m unittest discover -s tests
python -m mockcheck demo --out out/demo --site docs/data/reports.json
python scripts/build_cookbook.py
node -e "const s=require('fs').readFileSync('docs/index.html','utf8');new Function(s.match(/<script>([\s\S]*?)<\/script>/)[1]);new Function(require('fs').readFileSync('docs/phase2.js','utf8'));new Function(require('fs').readFileSync('docs/phase3.js','utf8'));for(const f of require('fs').readdirSync('docs/data').filter(x=>/^p[34]_/.test(x)))JSON.parse(require('fs').readFileSync('docs/data/'+f,'utf8'));console.log('site JS syntax ok')"
if [ -n "$(git status --porcelain docs/data COOKBOOK.md)" ]; then echo "generated files changed: review and commit"; git status --short docs/data COOKBOOK.md; else echo "generated files in sync"; fi
echo "check_all: OK"
