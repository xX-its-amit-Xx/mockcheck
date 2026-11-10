#!/bin/bash
# Deploy docs/ to Cloudflare Pages (project already exists; needs a logged-in wrangler).
set -euo pipefail
cd "$(dirname "$0")/.."
bash scripts/check_all.sh
npx wrangler pages deploy docs --project-name mockcheck --branch main --commit-dirty=true
