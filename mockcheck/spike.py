"""Score a Kraken2 report against a constructed spike-in truth, in absolute read pairs.

In a host-dominated sample (blood, plasma) percentages of the profile are meaningless, so a species is
called when at least N read pairs support it. Truth comes from results/spike_manifest.tsv."""
import csv
import os

from .parsers import genus, norm, read_kraken_counts
from .report import WATCH_GENERA, pct, title_case
from .truth import ROOT, load_contaminants, load_mock, load_organisms

SYN = {"bacillus subtilis": {"bacillus spizizenii"}}   # naming trap, see REPORT.md
HUMAN = "homo sapiens"
DEFAULT_MANIFEST = os.path.join(ROOT, "results", "spike_manifest.tsv")


def analyze_spike(path, sample, manifest=DEFAULT_MANIFEST, min_reads=10):
    load_mock("zymo-d6300-even")  # loads aliases
    with open(manifest, newline="") as f:
        rows = [r for r in csv.DictReader(f, delimiter="\t") if r["sample"] == sample]
    if not rows:
        raise ValueError(f"sample '{sample}' not found in {manifest}")
    sp, ge, tot = read_kraken_counts(path)
    truth, items = set(), []
    for r in rows:
        if r["organism"] == "(none)":
            continue
        o = norm(r["organism"])
        names = {o} | SYN.get(o, set())
        truth |= names
        got = sum(sp.get(n, 0) for n in names)
        exp = int(r["spiked_pairs"])
        items.append(dict(organism=o, fraction=float(r["fraction"]), expected=exp, observed=got,
                          genus_reads=ge.get(genus(o), 0), detected=got >= min_reads, recovery=got / exp if exp else None))
    items.sort(key=lambda x: -x["fraction"])
    contam = load_contaminants()
    fps = sorted(((v, k) for k, v in sp.items() if k != HUMAN and k not in truth and v >= min_reads), reverse=True)
    n = tot["total"]
    return dict(sample=sample, min_reads=min_reads, items=items, total_pairs=n,
                human_frac=(sp.get(HUMAN, 0) / n if n else 0), unclassified_frac=(tot["unclassified"] / n if n else 0),
                false_positives=[(k, v) for v, k in fps],
                fp_pathogen_genus=[k for v, k in fps if genus(k) in WATCH_GENERA],
                fp_contaminant_genus=[k for v, k in fps if genus(k) in contam])


def render_spike(a, label=""):
    T = a["min_reads"]
    det = [i for i in a["items"] if i["detected"]]
    miss = [i for i in a["items"] if not i["detected"]]
    orgs = load_organisms()
    first = (f"- Added organisms detected: **{len(det)} of {len(a['items'])}**" if a["items"]
             else "- Human-only control: no microbes were added")
    h = [f"# mockcheck spike-in report{': ' + label if label else ''}", "",
         f"Truth tier **T2** (constructed from real reads). Sample **{a['sample']}**, {a['total_pairs']:,} read pairs, "
         f"{pct(100 * a['human_frac'], 1)} human. A species counts as detected with at least **{T}** read pairs.", "",
         "## Scorecard", "", first,
         f"- Wrong species (not added, not human): **{len(a['false_positives'])}**, of which {len(a['fp_pathogen_genus'])} in pathogen-containing genera and {len(a['fp_contaminant_genus'])} in listed reagent-contaminant genera", ""]
    if a["items"]:
        h += ["| Organism | Added | Fraction | Found at species | Found at genus | Detected |", "|---|---|---|---|---|---|"]
        h += [f"| {title_case(i['organism'])} | {i['expected']} | {i['fraction']:g} | {i['observed']} | {i['genus_reads']} | {'yes' if i['detected'] else 'no'} |" for i in a["items"]]
    S, N, C = [], [], []
    if det:
        lo = min(det, key=lambda i: i["fraction"])
        S.append(f"{len(det)} added organisms were detected; the rarest detected was added at fraction {lo['fraction']:g} ({lo['expected']} read pairs, {lo['observed']} found).")
    for i in miss:
        gnote = f" ({i['genus_reads']} pairs did reach the genus)" if i["genus_reads"] >= T else ""
        N.append(f"{title_case(i['organism'])} added at fraction {i['fraction']:g} ({i['expected']} pairs) was not detected at {T}+ pairs{gnote}. Silence at this level is not evidence of absence.")
        info = orgs.get(i["organism"])
        if info and i["fraction"] >= 1e-3:
            N.append(f"At a generous fraction this miss is likely about the database, not depth: {info['why_it_matters_for_detection']}")
    if a["false_positives"]:
        top = ", ".join(f"{title_case(k)} {v}" for k, v in a["false_positives"][:3])
        N.append(f"Wrong species appeared even though they were not added. Largest: {top} (read pairs).")
        C.append("Run the human-only control through the same pipeline: whatever it reports with no microbe present is your false-alarm floor.")
        if a["fp_contaminant_genus"]:
            N.append("Some wrong species are in listed reagent-contaminant genera, so they may be real background in the human reads or reagents.")
    C.append(f"Thresholds depend on depth: {T} read pairs is a different bar at {a['total_pairs']:,} pairs than at ten times the depth. Repeat at your own depth.")
    C.append("These reads come from clean isolates and a cell-line genome, so detection limits here are optimistic for real low-biomass samples (no kit contamination, no cell-free DNA fragmentation).")
    h += ["", "## What this supports", ""] + [f"- {x}" for x in S or ["Nothing strong."]]
    h += ["", "## What it does not support", ""] + [f"- {x}" for x in N or ["No failure seen at this threshold."]]
    h += ["", "## What to check next", ""] + [f"- {x}" for x in C] + [""]
    return "\n".join(h)
