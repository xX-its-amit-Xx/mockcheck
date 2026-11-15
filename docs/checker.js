"use strict";
/* In-browser "paste your output" checker. A port of mockcheck/parsers.py and score.py (same maths, same names).
   Everything runs locally in the page; nothing is uploaded. Pure functions are exported for tests (node). */
(function (root) {
  let ALIASES = {};
  const setAliases = a => { ALIASES = a || {}; };
  const norm = n => { const s = n.replace(/_/g, " ").trim().toLowerCase(); return ALIASES[s] || s; };
  const genus = n => { const g = n.split(/\s+/)[0] || n; return ALIASES[g] || g; };
  const renorm = d => { const s = Object.values(d).reduce((a, b) => a + b, 0); const o = {}; if (s > 0) for (const k in d) o[k] = 100 * d[k] / s; return o; };
  const lines = t => t.split(/\r?\n/);

  function detectFormat(text) {
    const head = lines(text).slice(0, 40);
    const h = head.join("\n");
    if (h.includes("fraction_total_reads")) return "bracken";
    if (h.includes("#mpa_") || h.includes("s__") || h.includes("relative_abundance")) return "metaphlan";
    for (const l of head) {
      const p = l.split("\t");
      if (p.length >= 6 && !isNaN(parseFloat(p[0])) && /^\d+$/.test(p[1].trim()) && /^\d+$/.test(p[2].trim())) return "kraken2";
    }
    if (head.filter(l => l.trim()).every(l => l.split("\t").length >= 2)) return "generic";
    throw new Error("Cannot recognise the format. Pick one from the list.");
  }

  function readKraken(text) {
    const sp = {}; let total = 0, unclass = 0;
    for (const l of lines(text)) {
      const p = l.split("\t"); if (p.length < 6) continue;
      const reads = parseInt(p[1], 10), rank = p[3], name = p[5].trim();
      if (rank === "U") unclass = reads;
      else if (rank === "R" && name === "root") total = reads;
      else if (rank === "S") sp[norm(name)] = (sp[norm(name)] || 0) + reads;
    }
    const all = total + unclass, s = Object.values(sp).reduce((a, b) => a + b, 0);
    return { tool: "kraken2", species: renorm(sp), meta: { unclassified_frac: all ? unclass / all : null, species_level_frac: all ? s / all : null } };
  }

  function readBracken(text) {
    const L = lines(text), header = L[0].split("\t");
    const ni = header.indexOf("name"), fi = header.indexOf("fraction_total_reads");
    const sp = {}; let frac = 0;
    for (const l of L.slice(1)) {
      const p = l.split("\t"); if (p.length <= Math.max(ni, fi)) continue;
      const n = norm(p[ni]), v = parseFloat(p[fi]); if (isNaN(v)) continue;
      sp[n] = (sp[n] || 0) + v; frac += v;
    }
    return { tool: "bracken", species: renorm(sp), meta: { species_level_frac: frac } };
  }

  function readMetaphlan(text) {
    const sp = {};
    for (const l of lines(text)) {
      if (l.startsWith("#") || !l.includes("s__")) continue;
      const p = l.split("\t"); if (p.length < 3) continue;
      const last = p[0].split("|").pop();
      if (last.startsWith("s__") && !p[0].includes("t__")) sp[norm(last.slice(3))] = parseFloat(p[2]);
    }
    return { tool: "metaphlan4", species: renorm(sp), meta: {} };
  }

  function readGeneric(text) {
    const sp = {};
    for (const l of lines(text)) {
      const p = l.split("\t"); if (p.length < 2) continue;
      const v = parseFloat(p[1]); if (isNaN(v) || p[1].trim() === "") continue;
      sp[norm(p[0])] = (sp[norm(p[0])] || 0) + v;
    }
    return { tool: "generic", species: renorm(sp), meta: {} };
  }

  const READERS = { kraken2: readKraken, bracken: readBracken, metaphlan: readMetaphlan, metaphlan4: readMetaphlan, generic: readGeneric };
  function readProfile(text, fmt) {
    fmt = !fmt || fmt === "auto" ? detectFormat(text) : fmt;
    const prof = READERS[fmt](text);
    if (!Object.keys(prof.species).length) throw new Error("No species-level rows found. Is this a " + fmt + " file, and is it empty?");
    return prof;
  }

  const SWEEP = [0, 0.01, 0.05, 0.1, 0.5, 1, 2];
  const bray = (a, b) => {
    const keys = new Set([...Object.keys(a), ...Object.keys(b)]); let num = 0, den = 0;
    for (const k of keys) { const x = a[k] || 0, y = b[k] || 0; num += Math.abs(x - y); den += x + y; }
    return den ? num / den : NaN;
  };
  const called = (p, c) => new Set(Object.keys(p).filter(k => p[k] > 0 && p[k] >= c));
  function prf(pred, truth) {
    let tp = 0; for (const k of pred) if (truth.has(k)) tp++;
    const fp = pred.size - tp, fn = truth.size - tp;
    const pr = tp + fp ? tp / (tp + fp) : 0, rc = tp + fn ? tp / (tp + fn) : 0;
    return { tp, fp, fn, precision: pr, recall: rc, f1: pr + rc ? 2 * pr * rc / (pr + rc) : 0 };
  }

  /* mock: {kind, dna, cell, species:[...], anchors}. Same maths as score.analyze. */
  function analyze(profile, mock, cutoff) {
    const prof = profile.species, truth = new Set(mock.species);
    const hasCell = mock.cell && Object.keys(mock.cell).length > 0;
    const basis = profile.tool === "metaphlan4" && hasCell ? "cell" : "dna";
    const expected = basis === "cell" ? mock.cell : mock.dna;
    const out = { tool: profile.tool, kind: mock.kind, cutoff, basis, n_species_reported: Object.keys(prof).length, meta: profile.meta };
    out.at_cutoff = prf(called(prof, cutoff), truth);
    out.sweep = SWEEP.map(c => Object.assign({ cutoff: c }, prf(called(prof, c), truth)));
    out.false_positive_mass_pct = Object.keys(prof).filter(k => !truth.has(k)).reduce((a, k) => a + prof[k], 0);
    out.false_positives = Object.keys(prof).filter(k => !truth.has(k) && prof[k] >= cutoff).map(k => [k, prof[k]]).sort((a, b) => b[1] - a[1]);
    out.missed = [...truth].filter(k => (prof[k] || 0) < cutoff).sort();
    out.detected = [...truth].filter(k => (prof[k] || 0) >= cutoff).sort();
    if (mock.kind === "quantitative" && expected && Object.keys(expected).length) {
      out.bray_curtis = bray(prof, expected);
      out.per_organism = Object.keys(expected).sort((a, b) => expected[b] - expected[a]).map(o => {
        const ev = expected[o], ob = prof[o] || 0; return { organism: o, expected: ev, observed: ob, fold: ev ? ob / ev : null };
      });
    }
    if (mock.anchors && Object.keys(mock.anchors).length) out.anchors = Object.keys(mock.anchors).map(o => ({ organism: o, expected: mock.anchors[o], observed: prof[o] || 0 }));
    out.spizizenii_note = truth.has("bacillus subtilis") && (prof["bacillus subtilis"] || 0) < cutoff && (prof["bacillus spizizenii"] || 0) >= cutoff;
    return out;
  }

  const api = { setAliases, norm, genus, detectFormat, readProfile, analyze, SWEEP };
  if (typeof module !== "undefined" && module.exports) module.exports = api; else root.MockCheck = api;
})(typeof window !== "undefined" ? window : globalThis);
