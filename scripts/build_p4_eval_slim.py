"""Slim results/p4_eval_suite.json down to what the dashboard shows -> docs/data/p4_eval_slim.json."""
import json, os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
d = json.load(open(os.path.join(ROOT, "results/p4_eval_suite.json")))
S = ["dev", "test_geo", "test_host", "test_time", "blind"]
MODELS = ["ce", "ce_long", "hybrid", "knn", "wu_avg_long"]
BASE = ["majority", "virus+length_tree", "3mer_logreg", "CONTROL_3mer_logreg_train_labels_permuted"]
out = dict(tripwires=d["tripwires"], integrity=d["S1_integrity"], splits=S, heads={})
for h in ("host", "region"):
    rows = []
    for b in BASE:
        x = d["S2_shortcut_baselines"][h][b]
        rows.append(dict(name=b, kind="baseline", cells={s: dict(acc=x[s]["acc"], n=x[s]["n"]) for s in S if s in x}))
    for m in MODELS:
        x = d["S3_models"][m][h]
        rows.append(dict(name=m, kind="model", cells={s: dict(acc=x[s]["acc_mean"], n=x[s]["n"], utility=x[s]["utility_mean"], coverage=x[s]["coverage_mean"],
                         diff=(x[s]["acc_mean"] - x[s]["acc_3mer_same_rows"]) if x[s].get("acc_3mer_same_rows") is not None else None,
                         diff_ci=x[s].get("diff_vs_3mer_ci95")) for s in S if s in x}))
    out["heads"][h] = dict(rows=rows, by_similarity=d["S2_shortcut_baselines"][h]["_by_similarity_pooled_nondev"])
    sim = {}
    for m in ["ce", "ce_long", "hybrid", "knn"]:
        tot = {}
        for sp in ("test_geo", "test_host", "test_time", "blind"):
            x = d["S3_models"][m][h].get(sp)
            if x:
                for k, v in x["acc_by_similarity_to_train"].items():
                    if v["n"]:
                        t = tot.setdefault(k, [0, 0.0]); t[0] += v["n"]; t[1] += v["n"] * v["acc"]
        sim[m] = {k: dict(n=t[0], acc=t[1] / t[0]) for k, t in sorted(tot.items())}
    out["heads"][h]["model_by_similarity"] = sim
json.dump(out, open(os.path.join(ROOT, "docs/data/p4_eval_slim.json"), "w"), indent=1, allow_nan=False)
print("wrote docs/data/p4_eval_slim.json")
