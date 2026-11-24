"use strict";
/* Phase 4 evaluation audit section. Reads data/p4_eval_story.json (narrative) and data/p4_eval_slim.json (results).
   Uses helpers from index.html: $, el, j. */
(async function renderP4Eval() {
  const root = $("#p4eval-root");
  const S = await j("data/p4_eval_story.json"), D = await j("data/p4_eval_slim.json");
  if (!root || !S) return;
  const f = (x, d = 2) => (x == null || isNaN(Number(x)) ? "n/a" : Number(x).toFixed(d));
  const table = (head, rows) => el("div", { class: "tw" }, el("table", { class: "mini" },
    el("thead", {}, el("tr", {}, head.map(h => el("th", {}, h)))),
    el("tbody", {}, rows.map(r => el("tr", {}, r.map(c => el("td", {}, c == null ? "" : String(c))))))));
  root.append(el("p", { class: "mut small" }, "Updated " + S.updated + ". Summary first; tables are read live from the committed result file."));
  root.append(el("div", { class: "tldr" }, el("b", {}, "In short"), el("ul", {}, S.tldr.map(t => el("li", {}, t)))));
  if (D) {
    const NAMES = { majority: "Always the most common class", "virus+length_tree": "Virus identity + length", "3mer_logreg": "3-mer logistic regression (no neural net)",
      CONTROL_3mer_logreg_train_labels_permuted: "Control: 3-mer model, shuffled labels", ce: "Fine-tune (ce)", ce_long: "Fine-tune, longer (ce_long)", hybrid: "Hybrid", knn: "Nearest neighbour (knn)", wu_avg_long: "Utility-trained, longer (answers or abstains)" };
    const cell = c => {
      if (!c) return "";
      let t = `${f(c.acc)} (n=${c.n})`;
      if (c.diff != null && c.diff_ci && c.diff_ci[0] === c.diff_ci[0]) t += ` | vs 3-mer ${c.diff >= 0 ? "+" : ""}${f(c.diff)} [${f(c.diff_ci[0])}, ${f(c.diff_ci[1])}]`;
      return t;
    };
    const BINS = ["0-0.5", "0.5-0.7", "0.7-0.9", "0.9-1.0"];
    for (const h of ["host", "region"]) {
      const H = D.heads[h];
      root.append(el("div", { class: "card" }, el("h4", {}, (h === "host" ? "Host" : "Region") + ": accuracy by split on the rows the models are scored on"),
        table(["Method", ...D.splits], H.rows.map(r => [NAMES[r.name] || r.name, ...D.splits.map(s => cell(r.cells[s]))])),
        el("p", { class: "mut small" }, "Brackets: 95% interval of (model minus 3-mer baseline), resampling whole clusters, first seed. The utility-trained row is not compared because it can abstain."),
        el("h4", {}, "Accuracy against similarity to the nearest training sequence (all non-dev splits pooled)"),
        table(["Method", "<0.5", "0.5-0.7", "0.7-0.9", "0.9-1.0"], [
          ...["majority", "virus+length_tree", "3mer_logreg"].map(n => [NAMES[n], ...BINS.map(k => { const v = H.by_similarity[n][k]; return v && v.acc != null ? `${f(v.acc)} (n=${v.n})` : "n/a"; })]),
          ...["ce", "ce_long", "hybrid", "knn"].map(n => [NAMES[n], ...BINS.map(k => { const v = H.model_by_similarity[n][k]; return v ? `${f(v.acc)} (n=${v.n})` : "n/a"; })])])));
    }
    root.append(el("div", { class: "card" }, el("h4", {}, "Split integrity: how far is each holdout from training?"),
      table(["Split", "Sequences", "Clusters shared with training", "Exact copies of training", "Median similarity", "Share >= 0.9 similar"],
        D.splits.map(s => { const v = D.integrity[s]; return [s, v.n, `${v.clusters_shared_with_train} of ${v.n_clusters}`, v.seqs_exact_in_train, f(v.nn_sim_quantiles && v.nn_sim_quantiles["50"]), v.frac_sim_ge_0_9 == null ? "n/a" : Math.round(v.frac_sim_ge_0_9 * 100) + "%"]; })),
      el("details", {}, el("summary", {}, "All tripwires raised (" + D.tripwires.length + ")"), el("ul", {}, D.tripwires.map(t => el("li", {}, t))))));
  }
  for (const sec of S.sections) root.append(el("details", sec.title.startsWith("Mistakes") || sec.title.startsWith("What each") ? { open: "" } : {}, el("summary", {}, sec.title), sec.body.map(p => el("p", {}, p))));
})();
