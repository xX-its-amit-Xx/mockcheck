"""Phase 5b: variant-escape flag for an EXISTING detection assay (classification only; no design).

Question: given a guide whose perfect-match activity is known (a deployed assay's reference reading), and a variant target,
can we flag that the assay has effectively failed (activity indistinguishable from background), using only guide-target distance?
Data: adapt-seq-design Cas13a CCF-curated (MIT): 85 guides x ~225 targets plus negative controls.
Label: fail = measured logk <= 95th percentile of negative-control logk (the assay cannot be told from background).
Score: -(reference logk of the guide + predicted activity change). Models: ref_only (guide's reference, no distance),
H (mismatch count), D1 (position-weighted), ALL (D1+D2+D3 from p5_distance_eval). Ridge, out-of-fold.

PRE-REGISTERED (before running): D1/ALL "hold up" if, in every held-out split below, the paired guide-bootstrap 95% CI of
AUROC(D1 or ALL) - AUROC(H) is entirely above 0 AND AUROC(D1) - AUROC(ref_only) is entirely above 0. Splits:
held-out guides, held-out targets, held-out guides AND targets jointly (the strict generalisation test),
train Hamming<=3 -> test>=4. Random pairs shown only as a leaky reference.
Usage: python p5_escape_flag.py <CCF-curated dir> <out.json> [--boot 300]
"""
import csv, json, os, sys, re, collections
import numpy as np
from sklearn.linear_model import Ridge
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.model_selection import GroupKFold, KFold
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p5_distance_eval import feats, jsonable  # same distance definitions as the distance evaluation

SEED = 0


def load(d):
    T = [(r["name"], r["target_seq"].upper().replace("U", "T"))
         for r in csv.DictReader(open(os.path.join(d, "library.targets.tsv"), newline=""), delimiter="\t") if r["type"] == "exp"]
    allrows = list(csv.DictReader(open(os.path.join(d, "CCF_merged_pairs_annotated.curated.tsv"), newline=""), delimiter="\t"))
    neg = [float(r["out_logk_median"]) for r in allrows if r["type"] == "neg"]
    thr = float(np.percentile(neg, 95))
    rows = [r for r in allrows if r["type"] == "exp"]
    ref = collections.defaultdict(list)
    for r in rows:
        if int(r["guide_target_hamming_dist"]) == 0:
            ref[r["guide_seq"]].append(float(r["out_logk_median"]))
    cache, out = {}, []
    for r in rows:
        w = (r["target_before"] + r["target_at_guide"] + r["target_after"]).upper()
        if w not in cache:
            m = [n for n, s in T if w in s]
            cache[w] = m[0] if len(m) == 1 else None
        if cache[w] is None:
            continue
        g = r["guide_seq"]; y = float(r["out_logk_median"]); rf = float(np.median(ref[g]))
        out.append(dict(g=g, t=r["target_at_guide"], b=r["target_before"], a=r["target_after"], hd=int(r["guide_target_hamming_dist"]),
                        blk=r["crrna_block"], pos=int(r["guide_pos_nt"]), tgt=cache[w], ref=rf, yabs=y, y=y - rf, fail=int(y <= thr)))
    return out, thr, len(rows)


def make_splits(rows):
    n = len(rows); idx = np.arange(n)
    g = np.array([r["g"] for r in rows]); tg = np.array([r["tgt"] for r in rows]); hd = np.array([r["hd"] for r in rows])
    S = {}
    S["random_pairs (leaky reference)"] = [(a, b) for a, b in KFold(5, shuffle=True, random_state=SEED).split(idx)]
    S["held_out_guides"] = list(GroupKFold(5).split(idx, groups=g))
    S["held_out_targets"] = list(GroupKFold(5).split(idx, groups=tg))
    ug, ut = sorted(set(g)), sorted(set(tg))
    rng = np.random.default_rng(SEED)
    gf = dict(zip(ug, rng.permutation(len(ug)) % 5)); tf = dict(zip(ut, rng.permutation(len(ut)) % 5))
    gfa = np.array([gf[x] for x in g]); tfa = np.array([tf[x] for x in tg])
    S["held_out_guides_AND_targets (strict)"] = [(idx[(gfa != i) & (tfa != i)], idx[(gfa == i) & (tfa == i)]) for i in range(5)]
    S["train_hd<=3_test_hd>=4"] = [(idx[hd <= 3], idx[hd >= 4])]
    return S


def auc(y, s):
    return roc_auc_score(y, s) if 0 < y.sum() < len(y) else np.nan


def main():
    d, out = sys.argv[1], sys.argv[2]
    boot = int(sys.argv[sys.argv.index("--boot") + 1]) if "--boot" in sys.argv else 300
    rows, thr, nraw = load(d)
    print(f"rows mapped to a single target: {len(rows)}/{nraw}; fail threshold logk <= {thr:.3f}; base rate {np.mean([r['fail'] for r in rows]):.3f}", flush=True)
    y = np.array([r["y"] for r in rows]); ref = np.array([r["ref"] for r in rows]); fail = np.array([r["fail"] for r in rows])
    hd = np.array([r["hd"] for r in rows]); g = np.array([r["g"] for r in rows]); tg = np.array([r["tgt"] for r in rows])
    F = feats(rows); F = {k: F[k] for k in ("H", "D1", "ALL")}
    guides = sorted(set(g)); gidx = {k: np.where(g == k)[0] for k in guides}
    rng = np.random.default_rng(SEED)
    meths = ["ref_only", "H", "D1", "ALL"]
    res = {}
    for name, folds in make_splits(rows).items():
        score = {k: np.full(len(y), np.nan) for k in meths}
        tp = {k: 0 for k in meths}; fp = {k: 0 for k in meths}; fn = {k: 0 for k in meths}
        for tr, te in folds:
            score["ref_only"][te] = -ref[te]
            for k in ("H", "D1", "ALL"):
                X = F[k]; mu, sd = X[tr].mean(0), X[tr].std(0) + 1e-9
                m = Ridge(alpha=1.0).fit((X[tr] - mu) / sd, y[tr])
                score[k][te] = -(ref[te] + m.predict((X[te] - mu) / sd))
                s_tr = -(ref[tr] + m.predict((X[tr] - mu) / sd))
            for k in meths:  # operating point: threshold giving 90% recall of failures on TRAIN rows only
                s_tr = -ref[tr] if k == "ref_only" else -(ref[tr] + Ridge(alpha=1.0).fit((F[k][tr] - F[k][tr].mean(0)) / (F[k][tr].std(0) + 1e-9), y[tr]).predict((F[k][tr] - F[k][tr].mean(0)) / (F[k][tr].std(0) + 1e-9)))
                pos = s_tr[fail[tr] == 1]
                th = np.percentile(pos, 10) if len(pos) else np.inf
                t2 = te[hd[te] >= 1]
                pr = score[k][t2] >= th
                tp[k] += int((pr & (fail[t2] == 1)).sum()); fp[k] += int((pr & (fail[t2] == 0)).sum()); fn[k] += int((~pr & (fail[t2] == 1)).sum())
        test = ~np.isnan(score["H"]) & (hd >= 1)
        point = {k: float(auc(fail[test], score[k][test])) for k in meths}
        ap = {k: float(average_precision_score(fail[test], score[k][test])) for k in meths}
        op = {k: dict(recall=tp[k] / max(1, tp[k] + fn[k]), precision=tp[k] / max(1, tp[k] + fp[k])) for k in meths}
        dH = {k: [] for k in ("D1", "ALL")}; dR = {k: [] for k in ("D1", "ALL")}; ci = {k: [] for k in meths}
        for _ in range(boot):
            ii = np.concatenate([gidx[guides[j]] for j in rng.choice(len(guides), len(guides))]); ii = ii[test[ii]]
            if len(ii) < 20 or fail[ii].sum() in (0, len(ii)):
                continue
            a = {k: auc(fail[ii], score[k][ii]) for k in meths}
            for k in meths: ci[k].append(a[k])
            for k in ("D1", "ALL"):
                dH[k].append(a[k] - a["H"]); dR[k].append(a[k] - a["ref_only"])
        q = lambda v: [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))] if v else [np.nan, np.nan]
        res[name] = dict(n_test=int(test.sum()), n_fail_test=int(fail[test].sum()), auroc=point, auprc=ap, at_90pct_train_recall=op,
                         auroc_ci95={k: q(v) for k, v in ci.items()}, diff_vs_H_ci95={k: q(v) for k, v in dH.items()},
                         diff_vs_ref_only_ci95={k: q(v) for k, v in dR.items()}, boot_used=len(ci["H"]))
    held = [s for s in res if "leaky" not in s]
    verdict = {k: dict(holds_up=all(res[s]["diff_vs_H_ci95"][k][0] > 0 and res[s]["diff_vs_ref_only_ci95"][k][0] > 0 for s in held),
                       per_split_beats_H=[s for s in held if res[s]["diff_vs_H_ci95"][k][0] > 0],
                       per_split_beats_ref_only=[s for s in held if res[s]["diff_vs_ref_only_ci95"][k][0] > 0]) for k in ("D1", "ALL")}
    json.dump(jsonable(dict(fail_threshold_logk=thr, base_rate=float(fail.mean()), n_rows=len(rows), splits=res, verdict=verdict)), open(out, "w"), indent=1, allow_nan=False)
    for s, v in res.items():
        print(f"\n{s}  n_test={v['n_test']} fails={v['n_fail_test']} boots={v['boot_used']}")
        for k in meths:
            extra = f" vsH{[round(x,3) for x in v['diff_vs_H_ci95'][k]]} vsRef{[round(x,3) for x in v['diff_vs_ref_only_ci95'][k]]}" if k in v["diff_vs_H_ci95"] else ""
            print(f"  {k:8s} AUROC={v['auroc'][k]:.3f} CI{[round(x,3) for x in v['auroc_ci95'][k]]} AUPRC={v['auprc'][k]:.3f} "
                  f"@90%recall: recall={v['at_90pct_train_recall'][k]['recall']:.2f} precision={v['at_90pct_train_recall'][k]['precision']:.2f}{extra}")
    print("\nVERDICT", json.dumps(verdict))


if __name__ == "__main__":
    main()
