# Phase 3 brief (from Amit, 2026-11-12, saved verbatim in substance)

Questions
- Key decisions for evaluations: geographic accuracy, host compatibility. What makes each evaluation biologically useful for a company like Valthos?
- Look at competitors: which evaluation strategies are most useful?
- Add everything to the dashboard (digestible, dropdown, nested dropdown) and to the repo.
- Try several evaluation frameworks and see how each holds up. Track agent traces for each attempt, using real packages, real sandboxes, real tools, to show what real agents try to do and how to measure their performance.
- Which factors matter most for screening proprietary, fine-tuned or suite-of-models protein language models? At each step, how to benchmark efficiently?

Benchmark types to cover
1. Model development: is it working at all (tool-call checks, output checks, agent-trace checks that each step ran).
2. Gold data: does it agree with ground truth we set.
3. Practical: lab or real-world samples; what that looks like for infectious disease.

Context: human genetics predicts molecular mechanism of a variant. Infectious disease adds prediction of how viruses evolve, so we can stop them. Look at Valthos's mission and other companies': proprietary data, how they validate, what makes the difference in evals, which public evals exist.
