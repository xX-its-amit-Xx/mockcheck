# Phase 4 brief (from Amit, 2026-11-12)

Go ahead with a small fine-tune plus RL. Most important: a well-curated ground truth that is true to Valthos and related biosecurity companies, with a holdout the models are never trained on but should generalise to across many categories.

Needed
1. Decide the label categories (geographic location, host compatibility, etc.).
2. Find samples that carry labels for all categories.
3. Make sure the labels are actually ground truth.
In some cases create synthetic/derived labels or accept incomplete evals: samples rich in some labels but missing others.

Our guardrails: classification of existing public sequences only; no generation or design; no select-agent viruses; GISAID not used; Explorer fair-use rules for any GPU job.

## Addendum (from Amit, 2026-11-15): redirect toward clinical endpoints
- Pull in as much data as possible across multiomics or clinical trials, ideally. Aim: a deep understanding of these sequences so we are best prepared.
- Ideally study the endpoints relevant to infectious disease, not pure omics labels such as geographic location or host compatibility.
- Approved: commit/push/deploy of the checker (done, 3cfb9a5) and GPU follow-ups, "really go hard and see what you can break or find".
- Guardrails unchanged: classification of existing public sequences only; no generation or design; no select-agent viruses; no GISAID; Explorer fair-use rules.

## Addendum 2 (from Amit, 2026-11-24)
- "Make genuine evals; exhaustively at every stage from geographic source of origin to host compatibility" -> done as the Phase 4 audit (scripts/p4_eval_suite.py, dashboard #p4eval).
- "Go ahead" on cluster runs: cluster-held-out CV, robustness perturbations, shuffled-label control, more seeds (scripts/p4_robust.py, .sbatch).
- "Try and see if you can use frontier models like Evo or Marin-DNA or Mol-1 etc. to incorporate better inductive bias."
- "BioEmu2 and others are also great to include in a sandbox."
- Notes from checking (2026-11-24): Evo 2 excluded eukaryotic-host viruses from training (our three viruses are all eukaryotic-host) and is DNA, our data are proteins; BioEmu-1 is MIT (BioEmu-2 only on unverified third-party mirrors); marin-dna has a 1B model on Hugging Face, provenance unclear; "Mol-1" not identified. Use as embedding/feature extractors only (classification), never generation.
