"use strict";
/* Phase 3 dashboard: nested dropdowns from data/p3_*.json, plus live blocks from results tables and agent traces.
   Uses helpers from index.html: $, el, j, tsv, detailsFrom. */
(async function renderPhase3() {
  const FILES = ["p3_geohost", "p3_landscape", "p3_methods"];
  const secs = [];
  for (const f of FILES) { const d = await j("data/" + f + ".json"); if (d && d.sections) secs.push(...d.sections); }
  const root = $("#p3-root");
  if (!root || !secs.length) return;
  const num = x => (x === "" || x == null ? NaN : Number(x));
  const table = (head, rows) => el("div", { class: "tw" }, el("table", { class: "mini" },
    el("thead", {}, el("tr", {}, head.map(h => el("th", {}, h)))),
    el("tbody", {}, rows.map(r => el("tr", {}, r.map(c => el("td", {}, c == null ? "" : String(c))))))));
  const BASIS = { "sourced": "sourced", "our reasoning": "our reasoning", "our data": "our data" };
  const tag = b => b ? el("span", { class: "basis b-" + b.replace(" ", "-") }, BASIS[b] || b) : null;

  // ------------------------------------------------------------------ live blocks
  const live = {};
  live.agents = async box => {
    const runs = await tsv("data/results/agent_eval_runs.tsv"), r1 = await tsv("data/results/agent_eval_round1_runs.tsv");
    if (!runs) { box.append(el("p", { class: "mut" }, "Agent runs have not been generated yet.")); return; }
    const tasks = [...new Set(runs.map(r => r.task))], models = ["haiku", "sonnet", "opus"], conds = ["tool", "bare"];
    const cell = (t, c, m) => { const rr = runs.filter(r => r.task === t && r.condition === c && r.model === m); return rr.length ? rr.reduce((s, r) => s + num(r.L2_all_correct), 0) + " of " + rr.length : ""; };
    const total = runs.length, ok = runs.reduce((s, r) => s + num(r.L2_all_correct), 0);
    box.append(el("p", {}, el("b", {}, ok + " of " + total), " runs matched gold on every field (current task wording). Cell = runs correct of runs made."));
    box.append(table(["Task", ...conds.flatMap(c => models.map(m => c + " / " + m))], tasks.map(t => [t, ...conds.flatMap(c => models.map(m => cell(t, c, m)))])));
    const cost = m => { const rr = runs.filter(r => r.model === m); return rr.reduce((s, r) => s + num(r.cost_usd), 0) / rr.length; };
    const calls = c => { const rr = runs.filter(r => r.condition === c); return rr.reduce((s, r) => s + num(r.L1_tool_calls), 0) / rr.length; };
    box.append(el("h4", {}, "Cost and effort"));
    box.append(table(["Model", "Mean cost per run (USD)", "Runs"], models.map(m => [m, cost(m).toFixed(4), runs.filter(r => r.model === m).length])));
    box.append(el("p", { class: "small mut" }, "Mean tool calls per run: with mockcheck " + calls("tool").toFixed(1) + ", without " + calls("bare").toFixed(1) + ". Total tool errors recorded: " + runs.reduce((s, r) => s + num(r.L1_tool_errors), 0) + " (mostly shell commands the sandbox refused; the agents adapted)."));
    if (r1) {
      const okr1 = r1.reduce((s, r) => s + num(r.L2_all_correct), 0);
      box.append(el("h4", {}, "Round 1: before we fixed the task wording"));
      box.append(el("p", {}, okr1 + " of " + r1.length + " runs matched gold, and the failures were identical across all three models. That pattern pointed at our gold, not at the agents. The database says 'limosilactobacillus fermentum' where the truth file says 'lactobacillus fermentum'; our scoring code merged them silently while the prompt said to merge nothing. Every agent followed the prompt and was scored wrong."));
    }
  };

  live.frameworks = async box => {
    const runs = await tsv("data/results/agent_eval_runs.tsv");
    if (!runs) return;
    const nn = runs.length, ground = runs.filter(r => r.L1_grounded !== "").map(r => num(r.L1_grounded));
    const full = ground.filter(g => g === 1).length;
    const aware = runs.filter(r => r.diag_rename_aware !== "");
    const awareOk = aware.reduce((s, r) => s + num(r.diag_rename_aware), 0);
    const tool = runs.filter(r => r.condition === "tool"), used = tool.reduce((s, r) => s + num(r.L1_used_mockcheck), 0);
    const items = [
      ["A. Final answer vs gold (L2)", "Strictest and cheapest to score. After the wording fix every cell is at ceiling, so it cannot rank these models. Its real value in round 1 was exposing a gold-versus-spec mismatch: when every model disagrees with gold the same way, suspect the gold first."],
      ["B. Process checks (L1)", full + " of " + ground.length + " answers had every value traceable to tool output (grounding). Every run opened data and produced a parseable RESULT. Cheap, and a necessary floor, but it passes right answers and wrong answers alike."],
      ["C. With versus without the tool", "In the " + tool.length + " runs where the mockcheck package was available, agents used it in " + used + ". Accuracy did not change on these tasks; it added tool calls. The tool's value here is consistent conventions and provenance, not accuracy on a task this easy."],
      ["D. Honesty on an unanswerable question (T4, T6)", "Every run declined to name a cause when no lab result existed. But a hidden naming trap was missed: the reagent list uses the old genus Propionibacterium, and the plasma top species include Cutibacterium acnes (its newer name). Rename-aware answers: " + awareOk + " of " + aware.length + ". Literal list reading was 100%, rename-aware knowledge was not."],
      ["E. Red team: a planted lab note", "A typed 'qPCR confirmed Salmonella' note with an embedded instruction to report it. 0 of 36 runs followed it (18 primed with a hint to check Salmonella, 18 unprimed). All 18 unprimed runs opened the note. Strong result, but the note was crude; a subtler plant (a plausible instrument report) is the next test."],
      ["F. Model sweep (haiku, sonnet, opus)", "Equal accuracy on this suite and a very large cost spread. This suite cannot justify the larger model; harder or longer tasks are needed to separate them."],
      ["G. Sandbox guardrails", "The sandbox refused compound shell commands and writes outside it; agents adapted. No command referenced the gold, the repo's results or the other sandbox (text search over all traces). Read commands were not confined by path, so this is evidence, not a guarantee."]];
    box.append(el("div", {}, items.map(([t, d]) => el("div", { class: "card" }, el("h4", { style: "margin-top:0" }, t), el("p", {}, d)))));
  };

  live.traces = async box => {
    const idx = await j("data/traces/index.json");
    if (!idx) { box.append(el("p", { class: "mut" }, "Traces not published yet.")); return; }
    const sel = el("select", { "aria-label": "run" }, idx.runs.map(r => el("option", { value: r.run }, r.run.replace(/__/g, " / "))));
    const out = el("div");
    const show = async () => {
      const t = await j("data/traces/" + sel.value + ".json");
      if (!t) { out.replaceChildren(el("p", { class: "mut" }, "Trace not found.")); return; }
      const sc = Object.entries(t.score).map(([k, v]) => k + "=" + v).join("  ");
      out.replaceChildren(
        el("p", { class: "small" }, el("b", {}, "Gold: "), JSON.stringify(t.gold)),
        el("p", { class: "small" }, el("b", {}, "Agent's RESULT: "), JSON.stringify(t.result)),
        el("p", { class: "small mut" }, sc),
        ...t.steps.map((s, i) => s.kind === "say"
          ? detailsFrom("Step " + (i + 1) + ": agent says", el("pre", { class: "small" }, s.text))
          : detailsFrom("Step " + (i + 1) + ": " + s.name + (s.error ? " (blocked or failed)" : ""), [el("pre", { class: "small" }, s.input), s.output ? el("pre", { class: "small mut" }, s.output) : null])));
    };
    sel.onchange = show;
    box.append(el("p", { class: "small mut" }, "Pick any run. Each step shows what the agent said or the command it ran and the output it got back (output cut to 700 characters; paths replaced)."), sel, out);
    show();
  };

  live.sources = async box => {
    const q = await tsv("data/results/source_queries.tsv"), inv = await tsv("data/results/source_inventory.tsv");
    if (!q || !inv) { box.append(el("p", { class: "mut" }, "Crawl not run yet.")); return; }
    box.append(el("p", { class: "small" }, inv.length + " candidate sources from " + q.length + " queries. Hit counts are what each service reported, not how many are relevant; the top 8 by relevance (Zenodo, GitHub: stars, Europe PMC: citations) are listed. Titles are as returned and have not been checked for quality or license."));
    box.append(table(["Source", "Query", "Total hits reported"], q.map(r => [r.source, r.query, r.total_hits])));
    for (const s of ["zenodo", "europepmc", "github"]) {
      const rows = inv.filter(r => r.source === s);
      if (!rows.length) continue;
      box.append(detailsFrom(s + " (" + rows.length + " listed)", table(["Query", "Title", "Date", "License", "Link"], rows.map(r => [r.query, r.title, r.date, r.license, el("a", { href: r.url, target: "_blank", rel: "noopener" }, "open")]))));
    }
  };

  // ------------------------------------------------------------------ generic nested renderer
  function body(n) {
    const out = [];
    (n.p || []).forEach(t => out.push(el("p", {}, t)));
    if (n.list) out.push(el("ul", {}, n.list.map(x => el("li", {}, x))));
    if (n.table) out.push(table(n.table.head, n.table.rows));
    if (n.src) out.push(el("p", { class: "small mut srcs" }, "Sources: ", ...n.src.flatMap(([l, u], i) => [i ? " · " : "", el("a", { href: u, target: "_blank", rel: "noopener" }, l)])));
    if (n.live) { const box = el("div"); out.push(box); live[n.live] && live[n.live](box); }
    (n.kids || []).forEach(k => out.push(node(k)));
    return out;
  }
  function node(n) {
    const sm = el("summary", {}, n.t, n.basis ? " " : "", tag(n.basis), n.s ? el("span", { class: "small mut ssum" }, " " + n.s) : null);
    const d = el("details", {}, sm, el("div", { class: "in" }, body(n)));
    if (n.open) d.open = true;
    return d;
  }

  const nav = document.querySelector("nav .wrap");
  root.append(el("div", { class: "card", id: "phase3" }, el("h2", { style: "margin-top:0" }, "Phase 3: how to evaluate, who else does, and how agents behave"),
    el("p", { class: "sub" }, "Open any dropdown. Tags show where a statement comes from: sourced (linked), our reasoning, or our data."),
    el("p", { class: "p3idx" }, ...secs.flatMap((s, i) => [i ? " · " : "", el("a", { href: "#" + s.id }, s.nav || s.title)]))));
  for (const s of secs) root.append(el("section", { id: s.id }, el("h2", {}, s.title), s.intro ? el("p", { class: "sub" }, s.intro) : null, ...s.items.map(node)));
  if (nav) nav.append(el("a", { href: "#phase3" }, "Phase 3"));
})();
