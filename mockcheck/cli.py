import argparse
import glob
import json
import os
import sys

from . import __version__
from .parsers import read_profile
from .report import render_isolate, render_markdown, build_findings
from .score import analyze, analyze_isolate
from .truth import ROOT, list_mocks, load_mock


def _label(p):
    return os.path.basename(p)


def cmd_mocks(_):
    for k, t, kind in list_mocks():
        print(f"{k:18s} {kind:13s} {t}")
    print("custom: pass --truth FILE.tsv with columns organism, expected_percent [, expected_cell_percent]")


def cmd_report(a):
    prof = read_profile(a.profile, a.format)
    if a.isolate:
        text = render_isolate(analyze_isolate(prof, a.isolate, a.cutoff), _label(a.profile))
    else:
        mock = load_mock(a.truth or a.mock)
        text = render_markdown(analyze(prof, mock, a.cutoff), mock, _label(a.profile))
    if a.out in (None, "-"):
        sys.stdout.write(text)
    else:
        os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
        with open(a.out, "w", encoding="utf-8", newline="\n") as f:
            f.write(text)
        print(f"wrote {a.out}")


def cmd_compare(a):
    mock = load_mock(a.truth or a.mock)
    rows = []
    for p in a.profiles:
        r = analyze(read_profile(p, a.format), mock, a.cutoff)
        rows.append((r.get("bray_curtis", float("nan")), _label(p), r))
    rows.sort(key=lambda t: (t[0] != t[0], t[0]))
    print(f"| Run | Tool | Found | False species | Precision | Recall | Bray-Curtis |\n|---|---|---|---|---|---|---|")
    for bc, lab, r in rows:
        c = r["at_cutoff"]
        print(f"| {lab} | {r['tool']} | {c['tp']}/{len(mock.species)} | {c['fp']} | {c['precision']:.2f} | {c['recall']:.2f} | {bc:.3f} |")


def _demo_runs(raw):
    runs = []
    for p in sorted(glob.glob(f"{raw}/kraken2/*.report.gz")) + sorted(glob.glob(f"{raw}/bracken/*.bracken")) + sorted(glob.glob(f"{raw}/metaphlan/*.profile.txt")):
        b = os.path.basename(p)
        for suf in (".report.gz", ".bracken", ".profile.txt"):
            if b.endswith(suf):
                stem = b[: -len(suf)]
        sample = stem.split(".", 1)[0]
        tool = os.path.basename(os.path.dirname(p))
        method = tool if tool == "metaphlan" else tool + "_" + stem.split(".", 1)[1].replace(".", "_")
        runs.append((sample, method, p))
    return runs


def cmd_demo(a):
    """Run the tool end to end on the committed raw classifier outputs; no cluster needed."""
    raw = os.path.join(ROOT, "results", "raw")
    if not os.path.isdir(raw):
        sys.exit("results/raw not found: run from a repo checkout")
    out = a.out
    os.makedirs(f"{out}/reports", exist_ok=True)
    summary, site = [], {"cutoff": a.cutoff, "version": __version__, "runs": []}
    for sample, method, path in _demo_runs(raw):
        prof = read_profile(path)
        lab = f"{sample} / {method}"
        fname = f"{sample}__{method}.md"
        if sample.startswith("iso_"):
            org = sample[4:].replace("_", " ")
            i = analyze_isolate(prof, org, a.cutoff)
            text = render_isolate(i, lab)
            summary.append((sample, method, "isolate", i["truth_pct"], i["other_ge_cutoff"], None))
            site["runs"].append(dict(sample=sample, method=method, kind="isolate", markdown=text, truth_pct=i["truth_pct"]))
        else:
            mock = load_mock("zymo-d6300-even" if sample == "even" else "zymo-d6310-log")
            r = analyze(prof, mock, a.cutoff)
            S, N, C, M = build_findings(r, mock)
            text = render_markdown(r, mock, lab)
            summary.append((sample, method, mock.kind, r["at_cutoff"]["tp"], r["at_cutoff"]["fp"], r.get("bray_curtis")))
            site["runs"].append(dict(sample=sample, method=method, kind=mock.kind, supports=S, not_supported=N, check_next=C,
                                     microbiology=M, scorecard=dict(found=r["at_cutoff"]["tp"], of=len(mock.species), false_species=r["at_cutoff"]["fp"],
                                                                      precision=round(r["at_cutoff"]["precision"], 3), recall=round(r["at_cutoff"]["recall"], 3),
                                                                      bray_curtis=None if "bray_curtis" not in r else round(r["bray_curtis"], 3),
                                                                      false_mass_pct=round(r["false_positive_mass_pct"], 2)), markdown=text))
        with open(f"{out}/reports/{fname}", "w", encoding="utf-8", newline="\n") as f:
            f.write(text)
    with open(f"{out}/summary.tsv", "w", newline="\n") as f:
        f.write("sample\tmethod\tkind\tfound_or_truth_pct\tfalse_species\tbray_curtis\n")
        for s in summary:
            f.write("\t".join("" if v is None else (f"{v:.3f}" if isinstance(v, float) else str(v)) for v in s) + "\n")
    if a.site:
        os.makedirs(os.path.dirname(os.path.abspath(a.site)), exist_ok=True)
        with open(a.site, "w", encoding="utf-8", newline="\n") as f:
            json.dump(site, f, indent=1)
    print(f"demo: {len(summary)} runs scored, reports in {out}/reports, table in {out}/summary.tsv")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="mockcheck", description="Check a classifier run against a mock community and say what it supports.")
    ap.add_argument("--version", action="version", version=__version__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("mocks", help="list built-in mock communities").set_defaults(fn=cmd_mocks)

    def common(p):
        p.add_argument("--mock", default="zymo-d6300-even", help="built-in mock name (see: mockcheck mocks)")
        p.add_argument("--truth", help="your own truth TSV (organism, expected_percent[, expected_cell_percent])")
        p.add_argument("--format", default="auto", help="auto | kraken2 | bracken | metaphlan | generic")
        p.add_argument("--cutoff", type=float, default=0.1, help="species is 'called' at or above this percent (default 0.1)")

    r = sub.add_parser("report", help="write the supports / does-not-support / check-next report for one run")
    r.add_argument("profile")
    common(r)
    r.add_argument("--isolate", help="the run is a single organism, e.g. 'Escherichia coli' (negative control)")
    r.add_argument("--out", help="markdown file to write (default: print)")
    r.set_defaults(fn=cmd_report)
    c = sub.add_parser("compare", help="rank several runs against the same mock")
    c.add_argument("profiles", nargs="+")
    common(c)
    c.set_defaults(fn=cmd_compare)
    d = sub.add_parser("demo", help="run end to end on the committed raw outputs (no cluster, no downloads)")
    d.add_argument("--out", default="out/demo")
    d.add_argument("--cutoff", type=float, default=0.1)
    d.add_argument("--site", help="also write a JSON bundle for the dashboard, e.g. docs/data/reports.json")
    d.set_defaults(fn=cmd_demo)
    a = ap.parse_args(argv)
    try:
        a.fn(a)
    except (ValueError, FileNotFoundError) as e:
        sys.exit(f"mockcheck: {e}")
    return 0
