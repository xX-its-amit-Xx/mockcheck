"""Numbers only: presence/absence, abundance error, false-positive mass, cutoff sweep."""
from .parsers import genus, norm

SWEEP = [0.0, 0.01, 0.05, 0.1, 0.5, 1.0, 2.0]


def bray(a, b):
    keys = set(a) | set(b)
    num = sum(abs(a.get(k, 0) - b.get(k, 0)) for k in keys)
    den = sum(a.get(k, 0) + b.get(k, 0) for k in keys)
    return num / den if den else float("nan")


def called(prof, cutoff):
    return {k for k, v in prof.items() if v > 0 and v >= cutoff}


def prf(pred, truth):
    tp, fp, fn = len(pred & truth), len(pred - truth), len(truth - pred)
    pr = tp / (tp + fp) if tp + fp else 0.0
    rc = tp / (tp + fn) if tp + fn else 0.0
    return dict(tp=tp, fp=fp, fn=fn, precision=pr, recall=rc, f1=(2 * pr * rc / (pr + rc) if pr + rc else 0.0))


def analyze(profile, mock, cutoff=0.1):
    prof = profile.species
    truth = set(mock.species)
    basis = "cell" if (profile.tool == "metaphlan4" and mock.cell) else "dna"
    expected = mock.cell if basis == "cell" else mock.dna
    out = dict(tool=profile.tool, mock=mock.key, kind=mock.kind, cutoff=cutoff, basis=basis,
               n_species_reported=len(prof), meta=profile.meta)
    out["at_cutoff"] = prf(called(prof, cutoff), truth)
    out["sweep"] = [dict(cutoff=c, **prf(called(prof, c), truth)) for c in SWEEP]
    out["false_positive_mass_pct"] = sum(v for k, v in prof.items() if k not in truth)
    fps = sorted(((v, k) for k, v in prof.items() if k not in truth and v >= cutoff), reverse=True)
    out["false_positives"] = [(k, v) for v, k in fps]
    out["missed"] = sorted(k for k in truth if prof.get(k, 0) < cutoff)
    out["detected"] = sorted(k for k in truth if prof.get(k, 0) >= cutoff)
    if mock.kind == "quantitative" and expected:
        out["bray_curtis"] = bray(prof, expected)
        rows = []
        for org, ev in sorted(expected.items(), key=lambda kv: -kv[1]):
            ob = prof.get(org, 0.0)
            rows.append(dict(organism=org, expected=ev, observed=ob, fold=(ob / ev if ev else None)))
        out["per_organism"] = rows
        out["within_2x"] = [r["organism"] for r in rows if r["fold"] is not None and 0.5 <= r["fold"] <= 2]
        out["off_by_more_than_2x"] = [r for r in rows if r["fold"] is not None and not 0.5 <= r["fold"] <= 2 and r["observed"] > 0]
    if mock.anchors:
        out["anchors"] = [dict(organism=o, expected=ev, observed=prof.get(o, 0.0)) for o, ev in mock.anchors.items()]
    # naming trap: subtilis expected, spizizenii reported
    out["spizizenii_note"] = ("bacillus subtilis" in truth and prof.get("bacillus subtilis", 0) < cutoff
                              and prof.get("bacillus spizizenii", 0) >= cutoff)
    return out


def analyze_isolate(profile, organism, cutoff=0.1):
    prof = profile.species
    tr = norm(organism)
    others = {k: v for k, v in prof.items() if k != tr}
    top = max(others.items(), key=lambda kv: kv[1]) if others else ("", 0.0)
    return dict(tool=profile.tool, truth=tr, truth_pct=prof.get(tr, 0.0), other_mass_pct=100 - prof.get(tr, 0.0),
                other_ge_cutoff=sum(v >= cutoff for v in others.values()), top_other=top[0], top_other_pct=top[1],
                same_genus_other_pct=sum(v for k, v in others.items() if genus(k) == genus(tr)), cutoff=cutoff,
                meta=profile.meta)
