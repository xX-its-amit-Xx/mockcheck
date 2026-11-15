"use strict";
/* UI for the in-browser checker. Logic is in checker.js (MockCheck). Runs entirely in the browser. */
(async function renderChecker() {
  const root = document.getElementById("chk-root");
  if (!root || !window.MockCheck) return;
  const T = await j("data/checker_truth.json");
  if (!T) { root.append(el("p", { class: "mut" }, "Checker data could not be loaded.")); return; }
  MockCheck.setAliases(T.aliases);
  const f = (x, d = 2) => (x == null || isNaN(x) ? "n/a" : Number(x).toFixed(d));
  const pct = x => f(x * 100, 0) + "%";

  const mockSel = el("select", {}, ...Object.entries(T.mocks).map(([k, m]) => el("option", { value: k }, m.title)));
  const fmtSel = el("select", {}, ...[["auto", "Auto-detect"], ["kraken2", "Kraken2 report"], ["bracken", "Bracken output"], ["metaphlan", "MetaPhlAn profile"], ["generic", "Two columns: organism, percent"]].map(([v, t]) => el("option", { value: v }, t)));
  const cut = el("input", { type: "number", min: "0", step: "0.01", value: "0.1", style: "width:5em" });
  const ta = el("textarea", { rows: "8", style: "width:100%;font:13px ui-monospace,monospace", placeholder: "Paste a Kraken2 report, Bracken output, MetaPhlAn profile, or two-column TSV here (or choose a file)" });
  const file = el("input", { type: "file" });
  const out = el("div", { "aria-live": "polite" });

  function table(head, rows) {
    return el("table", {}, el("tr", {}, ...head.map(h => el("th", {}, h))), ...rows.map(r => el("tr", {}, ...r.map(c => el("td", {}, c)))));
  }
  function run() {
    out.replaceChildren();
    const text = ta.value;
    if (!text.trim()) return;
    try {
      const prof = MockCheck.readProfile(text, fmtSel.value);
      const mock = T.mocks[mockSel.value], c = Number(cut.value) || 0;
      const r = MockCheck.analyze(prof, mock, c);
      const a = r.at_cutoff;
      out.append(el("div", { class: "card" },
        el("p", {}, el("b", {}, "Read as: "), prof.tool, ` (${r.n_species_reported} species). Scored against: ${mock.title}, cutoff ${c}%.`),
        el("p", {}, `Of ${mock.species.length} organisms really in the tube, ${a.tp} found (recall ${pct(a.recall)}). ${a.fp} extra species called (precision ${pct(a.precision)}, F1 ${f(a.f1)}).`),
        el("p", {}, `Reads assigned to species that are not in the tube: ${f(r.false_positive_mass_pct)}% of the species-level profile.`),
        r.bray_curtis != null ? el("p", {}, `Bray-Curtis distance to the expected mix (${r.basis} basis, 0 = perfect, 1 = nothing shared): ${f(r.bray_curtis, 3)}.`) : null,
        r.spizizenii_note ? el("p", { class: "mut" }, "Note: B. subtilis was reported as B. spizizenii. Databases split these names; the strain is the same.") : null));
      if (r.missed.length) out.append(el("p", {}, el("b", {}, "Missed: "), r.missed.join(", ")));
      if (r.false_positives.length) out.append(el("details", {}, el("summary", {}, `${r.false_positives.length} extra species at or above the cutoff`),
        table(["Species", "% of profile"], r.false_positives.slice(0, 50).map(([k, v]) => [k, f(v, 3)]))));
      if (r.per_organism) out.append(el("details", { open: "" }, el("summary", {}, "Expected vs observed"),
        table(["Organism", "Expected %", "Observed %", "Fold"], r.per_organism.map(o => [o.organism, f(o.expected, 2), f(o.observed, 3), o.fold == null ? "n/a" : f(o.fold, 2) + "x"]))));
      if (r.anchors) out.append(el("details", { open: "" }, el("summary", {}, "Reliably known anchors"),
        table(["Organism", "Expected %", "Observed %"], r.anchors.map(o => [o.organism, o.expected, f(o.observed, 6)]))));
      out.append(el("details", {}, el("summary", {}, "Cutoff sweep"),
        table(["Cutoff %", "Found", "Extra", "Precision", "Recall"], r.sweep.map(s => [s.cutoff, s.tp, s.fp, pct(s.precision), pct(s.recall)]))));
    } catch (e) { out.append(el("p", { class: "tldr" }, String(e.message || e))); }
  }
  file.addEventListener("change", async () => { if (file.files[0]) { ta.value = await file.files[0].text(); run(); } });
  for (const c of [mockSel, fmtSel, cut]) c.addEventListener("change", run);
  ta.addEventListener("input", run);

  root.append(el("div", { class: "card" },
    el("p", {}, el("label", {}, "Mock community: ", mockSel), " ", el("label", {}, "Format: ", fmtSel), " ", el("label", {}, "Call a species present at or above ", cut, " %")),
    ta, el("p", {}, el("label", {}, "Or choose a file: ", file)), out));
})();
