"""Phase 5: do three guide-target 'distance' measures predict Cas13a activity loss on held-out data?

Data: adapt-seq-design (Broad, MIT) data/CCF-curated/CCF_merged_pairs_annotated.curated.tsv (85 guides x ~225 targets).
Target y = out_logk_median minus that guide's own median at its perfect-match (Hamming 0) targets: activity change
caused by divergence. Scored on rows with Hamming >= 1. No design, no generation: prediction on existing measurements.

PRE-REGISTERED PASS RULE (written before any run): a method "holds up" only if, in ALL THREE held-out splits
(held-out guides, held-out position blocks, train Hamming<=3 -> test >=4), its Spearman beats the plain Hamming-count
baseline with a paired guide-bootstrap 95% CI for the difference that lies entirely above 0. The random-pair split
is shown only as a leaky reference and does not count.

AMENDMENT (after run 1, disclosed): run 1 used a linear-in-count Hamming baseline (kept in
results/p5/distance_eval_run1_linear_hamming.json). Run 2 adds Hc = mismatch count as a category (nonlinear), a
leave-overlapping-guides-out split, within-guide Spearman and a within-guide label-permutation control. The STRICT rule
is: beat Hc (not just H) in all held-out splits, same paired-bootstrap criterion.

Distances:
  H  hamming            : mismatch count (baseline)
  D1 position-weighted  : per-position mismatch indicators -> learned weights (seed/centre sensitivity)
  D2 type+clustering    : which nucleotide substitution, plus adjacent-mismatch runs
  D3 k-mer spectrum     : 2/3/4-mer spectrum distance guide-window vs target window, plus flank GC
  ALL                   : D1+D2+D3 stacked
Usage: python p5_distance_eval.py <curated.tsv> <out.json> [--boot 300]
"""
import csv, json, sys, itertools, collections
import numpy as np
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold, KFold

NT = "ACGT"
TYPES = {p: i for i, p in enumerate(itertools.permutations(NT, 2))}  # 12 mismatch types
SEED = 0


def load(path):
    rows = [r for r in csv.DictReader(open(path, newline=""), delimiter="\t") if r["type"] == "exp"]
    ref = collections.defaultdict(list)
    for r in rows:
        if int(r["guide_target_hamming_dist"]) == 0:
            ref[r["guide_seq"]].append(float(r["out_logk_median"]))
    out = []
    for r in rows:
        g = r["guide_seq"]
        out.append(dict(g=g, t=r["target_at_guide"], b=r["target_before"], a=r["target_after"],
                        hd=int(r["guide_target_hamming_dist"]), blk=r["crrna_block"], pos=int(r["guide_pos_nt"]),
                        y=float(r["out_logk_median"]) - float(np.median(ref[g]))))
    return out


def kmer(s, k):
    v = np.zeros(4 ** k)
    for i in range(len(s) - k + 1):
        w = s[i:i + k]
        if all(c in NT for c in w):
            idx = 0
            for c in w:
                idx = idx * 4 + NT.index(c)
            v[idx] += 1
    return v


def gc(s):
    return sum(c in "GC" for c in s) / max(1, len(s))


def feats(rows):
    H, D1, D2, D3 = [], [], [], []
    Hc = np.array([[float(r["hd"] == k) for k in range(1, 10)] for r in rows])
    for r in rows:
        g, t = r["g"], r["t"]
        m = [g[i] != t[i] for i in range(len(g))]
        H.append([sum(m)])
        D1.append([float(x) for x in m])
        ty = np.zeros(12)
        for i in range(len(g)):
            if m[i] and (g[i], t[i]) in TYPES:
                ty[TYPES[(g[i], t[i])]] += 1
        runs = sum(1 for i in range(1, len(m)) if m[i] and m[i - 1])
        D2.append(list(ty) + [runs])
        d3 = []
        for k in (2, 3, 4):
            a, b = kmer(g, k), kmer(t, k)
            d3.append(np.abs(a - b).sum() / (len(g) - k + 1))
        d3 += [abs(gc(g) - gc(t)), gc(r["b"] + r["a"]), abs(gc(r["b"]) - gc(r["a"]))]
        D3.append(d3)
    H, D1, D2, D3 = map(np.array, (H, D1, D2, D3))
    return {"H": H, "Hc": Hc, "D1": D1, "D2": D2, "D3": D3, "ALL": np.hstack([D1, D2, D3])}


def splits(rows):
    n = len(rows)
    g = np.array([r["g"] for r in rows]); blk = np.array([r["blk"] for r in rows]); hd = np.array([r["hd"] for r in rows])
    idx = np.arange(n)
    S = {}
    S["random_pairs (leaky reference)"] = list(KFold(5, shuffle=True, random_state=SEED).split(idx))
    S["held_out_guides"] = list(GroupKFold(5).split(idx, groups=g))
    S["held_out_position_blocks"] = list(GroupKFold(5).split(idx, groups=blk))
    pos = {r["g"]: r["pos"] for r in rows}
    order = sorted(pos, key=pos.get); comp = {}; c = 0
    for i, k in enumerate(order):
        if i and pos[k] - pos[order[i - 1]] >= 28:
            c += 1
        comp[k] = c
    print("overlap clusters (guides within <28 nt chain together):", c + 1, flush=True)
    cl = np.array([comp[x] for x in g])
    S["held_out_overlap_clusters"] = list(GroupKFold(min(5, c + 1)).split(idx, groups=cl))
    S["train_hd<=3_test_hd>=4"] = [(idx[hd <= 3], idx[hd >= 4])]
    return S


def predict(F, y, folds, test_mask):
    """Return out-of-fold predictions (nan where a row is never in a test fold)."""
    pred = {k: np.full(len(y), np.nan) for k in list(F) + ["null"]}
    for tr, te in folds:
        pred["null"][te] = y[tr].mean()
        for k, X in F.items():
            mu, sd = X[tr].mean(0), X[tr].std(0) + 1e-9
            m = Ridge(alpha=1.0).fit((X[tr] - mu) / sd, y[tr])
            pred[k][te] = m.predict((X[te] - mu) / sd)
    return pred


def jsonable(o):
    """Replace NaN/inf with None so output is valid JSON (browsers reject bare NaN)."""
    if isinstance(o, dict):
        return {k: jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, float) and not np.isfinite(o):
        return None
    return o


def rho(a, b):
    return spearmanr(a, b)[0] if len(a) > 2 and np.std(a) > 0 else 0.0


def main():
    path, out = sys.argv[1], sys.argv[2]
    boot = int(sys.argv[sys.argv.index("--boot") + 1]) if "--boot" in sys.argv else 300
    rows = load(path)
    y = np.array([r["y"] for r in rows]); hd = np.array([r["hd"] for r in rows]); g = np.array([r["g"] for r in rows])
    F = feats(rows)
    guides = sorted(set(g)); gidx = {k: np.where(g == k)[0] for k in guides}
    rng = np.random.default_rng(SEED)
    res = {}
    for name, folds in splits(rows).items():
        pred = predict(F, y, folds, None)
        test = np.zeros(len(y), bool)
        for _, te in folds:
            test[te] = True
        test &= hd >= 1
        meths = ["H", "Hc", "D1", "D2", "D3", "ALL", "null"]
        point = {k: rho(pred[k][test], y[test]) for k in meths}
        rmse = {k: float(np.sqrt(np.mean((pred[k][test] - y[test]) ** 2))) for k in meths}
        diffs = {k: [] for k in meths if k != "H"}; diffs_c = {k: [] for k in meths if k not in ("H", "Hc")}
        cis = {k: [] for k in meths}
        for _ in range(boot):
            pick = rng.choice(len(guides), len(guides))
            ii = np.concatenate([gidx[guides[j]] for j in pick])
            ii = ii[test[ii]]
            r = {k: rho(pred[k][ii], y[ii]) for k in meths}
            for k in meths:
                cis[k].append(r[k])
                if k != "H":
                    diffs[k].append(r[k] - r["H"])
                if k not in ("H", "Hc"):
                    diffs_c[k].append(r[k] - r["Hc"])
        wg = {}
        for k in meths:
            vals = [rho(pred[k][ii], y[ii]) for ii in (gidx[x][test[gidx[x]]] for x in guides) if len(ii) > 5]
            wg[k] = float(np.mean(vals))
        q = lambda v: [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]
        res[name] = dict(n_test=int(test.sum()), spearman=point, within_guide_mean_spearman=wg, rmse=rmse,
                         spearman_ci95={k: q(v) for k, v in cis.items()},
                         diff_vs_hamming_ci95={k: q(v) for k, v in diffs.items()},
                         diff_vs_Hc_ci95={k: q(v) for k, v in diffs_c.items()})
    held = [k for k in res if "leaky" not in k]
    verdict = {}
    for k in ["D1", "D2", "D3", "ALL"]:
        wins = [res[s]["diff_vs_hamming_ci95"][k][0] > 0 for s in held]
        sw = [res[s]["diff_vs_Hc_ci95"][k][0] > 0 for s in held]
        verdict[k] = dict(beats_linear_hamming_in=[s for s, w in zip(held, wins) if w], holds_up_vs_linear=all(wins),
                          beats_categorical_hamming_in=[s for s, w in zip(held, sw) if w], holds_up_STRICT=all(sw))
    # permutation control: shuffle y within each guide (breaks feature-label link, keeps guide structure); expect rho ~ 0
    yp = y.copy()
    for k in guides:
        yp[gidx[k]] = rng.permutation(y[gidx[k]])
    fo = splits(rows)["held_out_guides"]
    pp = predict(F, yp, fo, None)
    mask = hd >= 1
    perm = {k: float(rho(pp[k][mask], yp[mask])) for k in ["H", "Hc", "D1", "D2", "D3", "ALL"]}
    print("PERMUTATION CONTROL (held_out_guides, y shuffled within guide), expect ~0:", {k: round(v, 3) for k, v in perm.items()})
    json.dump(jsonable(dict(permutation_control_held_out_guides=perm, splits=res, verdict=verdict, n_rows=len(rows), n_guides=len(guides), boot=boot)), open(out, "w"), indent=1, allow_nan=False)
    for s, v in res.items():
        print(f"\n{s}  n_test={v['n_test']}")
        for k in ["null", "H", "Hc", "D1", "D2", "D3", "ALL"]:
            d = v["diff_vs_Hc_ci95"].get(k)
            print(f"  {k:5s} rho={v['spearman'][k]:+.3f} CI{[round(x,3) for x in v['spearman_ci95'][k]]} rmse={v['rmse'][k]:.3f} within-guide={v['within_guide_mean_spearman'][k]:+.3f}"
                  + (f"  vs Hc diff CI{[round(x,3) for x in d]}" if d else ""))
    print("\nVERDICT", json.dumps(verdict))


if __name__ == "__main__":
    main()
