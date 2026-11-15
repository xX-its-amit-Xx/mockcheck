// Usage: node checker_parity.js <truth.json> <file>...  -> JSON {file: analysis} for the mock chosen by file name.
const fs = require("fs");
const zlib = require("zlib");
const M = require("../docs/checker.js");
const [truthPath, ...files] = process.argv.slice(2);
const T = JSON.parse(fs.readFileSync(truthPath, "utf8"));
M.setAliases(T.aliases);
const out = {};
for (const f of files) {
  const base = require("path").basename(f);
  const key = base.startsWith("even") ? "zymo-d6300-even" : "zymo-d6310-log";
  const prof = M.readProfile((f.endsWith(".gz") ? zlib.gunzipSync(fs.readFileSync(f)).toString("utf8") : fs.readFileSync(f, "utf8")), "auto");
  const r = M.analyze(prof, T.mocks[key], 0.1);
  out[f] = { tool: prof.tool, key, tp: r.at_cutoff.tp, fp: r.at_cutoff.fp, fn: r.at_cutoff.fn,
    fpmass: r.false_positive_mass_pct, bray: r.bray_curtis == null ? null : r.bray_curtis, nrep: r.n_species_reported };
}
console.log(JSON.stringify(out));
