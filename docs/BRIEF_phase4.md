# Phase 4 brief (from Amit, 2026-11-12)

Go ahead with a small fine-tune plus RL. Most important: a well-curated ground truth that is true to Valthos and related biosecurity companies, with a holdout the models are never trained on but should generalise to across many categories.

Needed
1. Decide the label categories (geographic location, host compatibility, etc.).
2. Find samples that carry labels for all categories.
3. Make sure the labels are actually ground truth.
In some cases create synthetic/derived labels or accept incomplete evals: samples rich in some labels but missing others.

Our guardrails: classification of existing public sequences only; no generation or design; no select-agent viruses; GISAID not used; Explorer fair-use rules for any GPU job.
