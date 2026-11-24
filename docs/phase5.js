"use strict";
/* Phase 5 dashboard: Cas13a guide-target distance and variant-escape flag. Reads data/p5_story.json (narrative),
   data/p5_distance_eval.json and data/p5_escape_flag.json (results). Uses helpers from index.html: $, el, j. */
(async function renderPhase5() {
  const root = $("#p5-root");
  const S = await j("data/p5_story.json"), D = await j("data/p5_distance_eval.json"), E = await j("data/p5_escape_flag.json");
  if (!root || !S) return;
  const f = (x, d = 2) => (x == null || isNaN(Number(x)) ? "n/a" : Number(x).toFixed(d));
  const ci = c => (c && !isNaN(c[0]) ? ` [${f(c[0])}, ${f(c[1])}]` : "");
  const table = (head, rows) => el("div", { class: "tw" }, el("table", { class: "mini" },
    el("thead", {}, el("tr", {}, head.map(h => el("th", {}, h)))),
    el("tbody", {}, rows.map(r => el("tr", {}, r.map(c => el("td", {}, c == null ? "" : String(c))))))));

  root.append(el("p", { class: "mut small" }, "Updated " + S.updated + ". Plain-language summary first; tables below are read live from the committed result files."));
  root.append(el("div", { class: "tldr" }, el("b", {}, "In short"), el("ul", {}, S.tldr.map(t => el("li", {}, t)))));

  if (D) {
    const names = { H: "Mismatch count (baseline)", D1: "Position-weighted", D2: "Substitution type + runs", D3: "k-mer spectrum", ALL: "All three combined" };
    const sp = Object.keys(D.splits).filter(k => !k.includes("leaky"));
    const box = el("div", { class: "card" }, el("h4", {}, "Task 1: how well does each distance predict activity change? (Spearman, higher is better)"),
      table(["Distance", ...sp, "Verdict"], ["H", "D1", "D2", "D3", "ALL"].map(k => [names[k], ...sp.map(s => f(D.splits[s].spearman[k]) + ci(D.splits[s].spearman_ci95[k])),
        k === "H" ? "baseline" : (D.verdict[k] && (D.verdict[k].holds_up_vs_linear) ? "holds up" : "does not hold up")])),
      el("p", { class: "mut small" }, "Brackets are 95% intervals from resampling guides. 'Holds up' = beats plain mismatch count in every held-out split with the paired interval above zero. Random-pair results are omitted because they leak guide identity."));
    root.append(box);
  }
  if (E) {
    const sp = Object.keys(E.splits);
    root.append(el("div", { class: "card" }, el("h4", {}, "Task 2: flag that a variant has pushed an existing assay to background (AUROC, higher is better)"),
      table(["Method", ...sp], ["ref_only", "H", "D1", "ALL"].map(k => [{ ref_only: "Guide's own reference only", H: "Mismatch count", D1: "Position-weighted", ALL: "All three combined" }[k],
        ...sp.map(s => f(E.splits[s].auroc[k]) + ci(E.splits[s].auroc_ci95[k]))])),
      el("p", { class: "mut small" }, "Strictest test: guides and targets both unseen ('held_out_guides_AND_targets')."),
      el("h4", {}, "What a usable alarm costs: set to catch 90% of failures (threshold fixed on training data only)"),
      table(["Split", "Position-weighted: recall", "precision", "Failure base rate in test"],
        sp.map(s => { const v = E.splits[s]; return [s, f(v.at_90pct_train_recall.D1.recall), f(v.at_90pct_train_recall.D1.precision), f(v.n_fail_test / v.n_test)]; })),
      el("p", {}, el("b", {}, "Reading it: "), "ranking quality is good, but at that catch rate roughly 7 of every 10 alarms are false, and in the extrapolation split the alarm flags everything (precision equals the base rate).")));
  }
  for (const sec of S.sections) {
    root.append(el("details", sec.title.startsWith("What failed") ? { open: "" } : {}, el("summary", {}, sec.title), sec.body.map(p => el("p", {}, p))));
  }
})();
