#!/usr/bin/env python3
"""Phase 4 robustness check (CPU): are the distance tiers an artefact of counting 5-mers?

Embeds every sequence with the frozen ESM-2 35M encoder (mean-pooled), finds each non-training sequence's nearest training sequence
of the same virus by cosine similarity, and re-assigns tiers by that measure. To compare fairly, the tier SIZES are matched to the
5-mer tiers (the nearest 'near' fraction, etc.), so only WHICH sequences land in which tier can change.
Writes results/p4_embedding_tiers.tsv (per-condition utility by embedding tier) and prints the agreement with the 5-mer tiers.
Usage: python scripts/p4_embed_tiers.py [--model PATH_OR_HF_NAME]"""
import argparse, collections, csv, gzip, os, sys
import numpy as np
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p4_analyze as A  # noqa: E402
import p4_train as T  # noqa: E402

ROOT = A.ROOT


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--model", default="facebook/esm2_t12_35M_UR50D"); a = ap.parse_args()
    from transformers import AutoModel, AutoTokenizer
    rows = T.load_rows(f"{ROOT}/data/groundtruth/p4_groundtruth.tsv.gz")
    tok, enc = AutoTokenizer.from_pretrained(a.model), AutoModel.from_pretrained(a.model).eval()
    emb = []
    with torch.no_grad():
        for i in range(0, len(rows), 16):
            b = tok([r["seq"] for r in rows[i:i + 16]], padding=True, truncation=True, max_length=T.MAXLEN + 2, return_tensors="pt")
            h = enc(**b).last_hidden_state; m = b["attention_mask"].unsqueeze(-1).float()
            emb.append(((h * m).sum(1) / m.sum(1)).numpy())
            if i % 320 == 0:
                print("embedded", i, flush=True)
    E = np.concatenate(emb); E = E / np.linalg.norm(E, axis=1, keepdims=True)
    acc2i = {r["acc"]: i for i, r in enumerate(rows)}
    train = collections.defaultdict(list)
    for i, r in enumerate(rows):
        if r["split"] == "train":
            train[r["virus"]].append(i)
    cos = {}
    for i, r in enumerate(rows):
        if r["split"] != "train":
            cos[r["acc"]] = float((E[train[r["virus"]]] @ E[i]).max())
    preds = A.read_preds(f"{ROOT}/results/p4_preds/ce__seed0.tsv.gz")
    nontrain = [p for p in preds if p["split"] not in ("train", "blind")]
    kmer_tier = {p["accession"]: A.tier_of(p) for p in nontrain}
    order = sorted(nontrain, key=lambda p: -cos[p["accession"]])
    sizes = collections.Counter(kmer_tier.values())
    emb_tier, k = {}, 0
    for name in ("near", "mid", "far"):
        for p in order[k:k + sizes[name]]:
            emb_tier[p["accession"]] = name
        k += sizes[name]
    agree = np.mean([emb_tier[x] == kmer_tier[x] for x in kmer_tier])
    print(f"tier sizes {dict(sizes)}; agreement between 5-mer and embedding tiers: {agree:.2f}")
    # utility by embedding tier for each condition (same tuned thresholds as the main analysis)
    out = []
    for f in sorted(os.listdir(f"{ROOT}/results/p4_preds")):
        cond, seed = f[:-7].rsplit("__seed", 1)
        if cond not in ("ce", "knn", "hybrid", "probe", "majority"):
            continue
        rows_p = A.read_preds(f"{ROOT}/results/p4_preds/{f}")
        thr = A.tune_threshold(rows_p) if cond in ("ce", "probe", "knn") else None
        for tname in ("near", "mid", "far"):
            sub = [p for p in rows_p if emb_tier.get(p["accession"]) == tname]
            u, _, _ = A.util_vec(sub, thr, thr)
            out.append(dict(condition=cond + ("+abstain" if thr is not None else ""), seed=seed, tier=tname, n=len(sub), utility=round(float(u.mean()), 4)))
    agg = collections.defaultdict(list)
    for o in out:
        agg[(o["condition"], o["tier"])].append(o["utility"])
    with open(f"{ROOT}/results/p4_embedding_tiers.tsv", "w", newline="\n") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n"); w.writerow(["condition", "tier_by_embedding", "seeds", "mean_utility", "tier_agreement_with_5mer"])
        for (c, t), v in sorted(agg.items()):
            w.writerow([c, t, len(v), round(float(np.mean(v)), 4), round(float(agree), 3)])
    for (c, t), v in sorted(agg.items()):
        print(f"{c:14s} {t:5s} {np.mean(v):7.3f}")


if __name__ == "__main__":
    main()
