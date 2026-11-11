"use strict";
/* Phase 2 dashboard: real data, decisions, innovation, data needs, data dictionary.
   Uses the helpers defined in index.html: $, el, j, tsv. Every table is read live from the committed files. */
(async function renderPhase2() {
  const P = await j("data/phase2.json");
  if (!P) return;
  const num = x => (x === "" || x == null ? NaN : Number(x));
  const median = a => { const b = a.filter(v => !isNaN(v)).sort((x, y) => x - y); return b.length ? b[Math.floor(b.length / 2)] : NaN; };
  const f = (x, d = 2) => (isNaN(x) ? "n/a" : Number(x).toFixed(d));
  const badge = t => el("span", { class: "tier t-" + t.toLowerCase() }, t);
  const table = (head, rows) => el("div", { class: "tw" }, el("table", { class: "mini" },
    el("thead", {}, el("tr", {}, head.map(h => el("th", {}, h)))),
    el("tbody", {}, rows.map(r => el("tr", {}, r.map(c => el("td", {}, c)))))));
  const bar = (v, max, label) => el("div", { class: "bar", title: label }, el("i", { style: "width:" + Math.max(0, Math.min(100, 100 * v / (max || 1))) + "%" }));

  // ---- tiers
  $("#rd-tiers").replaceChildren(...P.tiers.map(t => el("div", { class: "card tiercard" },
    el("h3", {}, badge(t.tier), " ", t.name),
    el("p", {}, t.what),
    el("p", {}, el("b", {}, "Truth: "), t.truth),
    el("p", {}, el("b", {}, "Used here: "), t.used),
    el("p", { class: "ok" }, el("b", {}, "Can show: "), t.shows),
    el("p", { class: "bad" }, el("b", {}, "Cannot show: "), t.cannot))));

  // ---- findings (written after the analysis, from docs/data/findings.json)
  const F = await j("data/findings.json");
  const fbox = $("#rd-findings");
  if (F && F.findings) {
    fbox.replaceChildren(...F.findings.map(x => el("div", { class: "card" },
      el("h3", {}, badge(x.tier), " ", x.title),
      el("p", {}, x.what),
      x.caveat ? el("p", { class: "small mut" }, el("b", {}, "Caveat: "), x.caveat) : null,
      x.evidence ? el("p", { class: "small mut" }, "Evidence: ", x.evidence) : null)));
  } else fbox.replaceChildren(el("p", { class: "mut" }, "Findings are written once the Phase 2 runs finish."));

  // ---- live: spike-in sensitivity
  const S = await tsv("data/results/spike_summary.tsv"), sbox = $("#rd-spike");
  if (S && S.length) {
    const sel = el("select", { "aria-label": "minimum read pairs" }, [1, 3, 10].map(t => el("option", { value: t }, "at least " + t + " read pair" + (t > 1 ? "s" : ""))));
    sel.value = "10";
    const out = el("div");
    const draw = () => {
      const t = sel.value, methods = [...new Set(S.map(r => r.method))];
      out.replaceChildren(...methods.map(m => {
        const rows = S.filter(r => r.method === m);
        const levels = [...new Set(rows.map(r => r.fraction))].sort((a, b) => b - a);
        return el("div", { class: "card" }, el("h4", {}, m),
          table(["Added at (fraction of reads)", "Read pairs added", "Bacteria detected", "", "Yeast detected", "", "Median recovery"],
            levels.map(l => {
              const b = rows.find(r => r.fraction === l && r.group === "bacteria") || {}, y = rows.find(r => r.fraction === l && r.group === "yeast") || {};
              const bv = num(b["detected_ge" + t]), yv = num(y["detected_ge" + t]);
              return [Number(l).toExponential(0), b.expected_pairs || "", isNaN(bv) ? "" : bv + " of " + b.n, bar(bv, num(b.n)), isNaN(yv) ? "" : yv + " of " + y.n, bar(yv, num(y.n)), f(num(b.median_recovery))];
            })));
      }));
    };
    sel.addEventListener("change", draw);
    sbox.replaceChildren(el("p", { class: "small mut" }, "Each row is one level at which a microbe was added to a human background. Pick how many read pairs a species needs before you call it detected."),
      el("label", {}, "Call a species when it has ", sel), out);
    draw();
  } else sbox.replaceChildren(el("p", { class: "mut" }, "Spike-in results are not in yet."));

  // ---- live: false positives
  const FP = await tsv("data/results/spike_false_positives.tsv"), pbox = $("#rd-fp");
  if (FP && FP.length) {
    const sel = el("select", { "aria-label": "minimum read pairs" }, [1, 3, 10].map(t => el("option", { value: t }, "at least " + t + " read pair" + (t > 1 ? "s" : ""))));
    sel.value = "10";
    const out = el("div");
    const draw = () => {
      const rows = FP.filter(r => r.min_reads === sel.value);
      out.replaceChildren(table(["Sample", "Tool", "Wrong species", "In pathogen-containing genera", "Largest wrong species (read pairs)"],
        rows.map(r => [r.sample === "spike_H0" ? "human only (no microbe added)" : r.sample, r.method, r.false_species, r.false_species_pathogen_genus, r.top3])));
    };
    sel.addEventListener("change", draw);
    pbox.replaceChildren(el("p", { class: "small mut" }, "Species the tool reported that were not added. The human-only row is the cleanest false-positive floor: nothing microbial was put in."),
      el("label", {}, "Count species with ", sel), out);
    draw();
  } else pbox.replaceChildren(el("p", { class: "mut" }, "False-positive results are not in yet."));

  // ---- live: real plasma
  const U = await tsv("data/results/uganda_summary.tsv"), T = await tsv("data/results/uganda_top_species.tsv"), A = await tsv("data/results/uganda_tool_agreement.tsv");
  const ubox = $("#rd-uganda");
  if (U && U.length) {
    const c1 = U.filter(r => r.confidence === "0.1"), rows = c1.length ? c1 : U;
    const hum = median(rows.map(r => num(r.human_frac))), spc = median(rows.map(r => num(r.species_ge10))), rc = median(rows.map(r => num(r.reagent_contaminant_genus_share)));
    const jac = A && A.length ? median(A.map(r => num(r.jaccard))) : NaN;
    ubox.replaceChildren(
      el("div", { class: "sc" }, [["runs", rows.length], ["median human share", f(100 * hum, 1) + "%"], ["median non-human species (10+ pairs)", f(spc, 0)], ["median share in listed reagent-contaminant genera", f(100 * rc, 0) + "%"], ["Kraken2 vs MetaPhlAn agreement (Jaccard)", f(jac)]].map(([k, v]) => el("div", {}, k + ": ", el("b", {}, String(v))))),
      table(["Run", "Pairs", "Human", "Species 10+ pairs", "Reagent-genus share", "Top species (pairs)"],
        rows.map(r => [r.run, r.pairs, f(100 * num(r.human_frac), 1) + "%", r.species_ge10, f(100 * num(r.reagent_contaminant_genus_share), 0) + "%",
          (T || []).filter(t => t.run === r.run && t.confidence === r.confidence).slice(0, 3).map(t => t.species + " " + t.pairs).join("; ")])));
  } else ubox.replaceChildren(el("p", { class: "mut" }, "Plasma results are not in yet."));

  // ---- decisions, innovation, data needs
  $("#dec-body").replaceChildren(...P.decisions.map(d => detailsFrom(d.id + " (" + d.date + "): " + d.title,
    [el("p", {}, el("b", {}, "Choice: "), d.choice), el("p", {}, el("b", {}, "Why: "), d.why), el("p", { class: "small mut" }, el("b", {}, "Cost: "), d.cost)])));
  $("#inn-body").replaceChildren(...P.innovation.map(x => detailsFrom(x.idea,
    [el("p", {}, x.why), el("p", { class: "small mut" }, el("b", {}, "First step to test it: "), x.first_step)])));
  $("#need-body").replaceChildren(...P.needs_data.map(x => detailsFrom(x.need,
    [el("p", {}, x.why), el("p", { class: "small mut" }, el("b", {}, "Where it could come from: "), x.who)])));

  // ---- data representation: live dictionary
  const D = await j("data/dictionary.json"), rbox = $("#rep-dict");
  const rep = P.representation;
  $("#rep-intro").replaceChildren(el("p", {}, rep.intro), el("ul", {}, rep.principles.map(p => el("li", {}, p))));
  if (D) {
    const sel = el("select", { "aria-label": "table" }, D.tables.map((t, i) => el("option", { value: i }, t.file)));
    const out = el("div", { class: "card" });
    const show = async () => {
      const t = D.tables[sel.value], name = t.file.split("/").pop();
      const rows = await tsv("data/results/" + name);
      out.replaceChildren(el("h4", {}, badge(t.tier), " ", t.file),
        el("p", {}, el("b", {}, "One row is: "), t.row), el("p", {}, el("b", {}, "Answers: "), t.answers),
        el("dl", { class: "dl" }, Object.entries(t.columns).flatMap(([k, v]) => [el("dt", {}, k), el("dd", {}, v)])),
        rows && rows.length ? [el("p", { class: "small mut" }, "First rows of the real file:"), table(Object.keys(rows[0]), rows.slice(0, 5).map(r => Object.values(r)))] : el("p", { class: "small mut" }, "This table has not been generated yet."));
    };
    sel.addEventListener("change", show);
    rbox.replaceChildren(el("label", {}, "Table ", sel), out);
    show();
  }
})();
