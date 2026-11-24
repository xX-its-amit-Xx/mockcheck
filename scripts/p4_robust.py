#!/usr/bin/env python3
"""Phase 4 robustness and novelty runs (GPU, resumable): cross-entropy fine-tune of ESM-2 35M, host + region heads.

Part A  cluster-held-out 5-fold CV over ALL tier-A rows: whole sequence clusters are held out together (a novelty test that is
        not just country/year/host). Different random fold assignment per seed.
Part B  original train split; test rows predicted under perturbations: clean, trunc_50, trunc_25 (keep first fraction),
        sub_5, sub_20 (random amino-acid substitutions), shuffle (residue order randomised, composition kept).
Part C  control: train with the TRAINING labels permuted, predict clean test rows (should be near chance).
Outputs (tsv.gz, skipped if present): OUT/cv__seed{s}.tsv.gz, OUT/pert__seed{s}.tsv.gz, OUT/labelshuf__seed{s}.tsv.gz
Classification of existing public sequences only; perturbations are for measuring robustness of a classifier, no sequence is generated for use.
  python p4_robust.py --data data/groundtruth/p4_groundtruth.tsv.gz --out OUT --base MODELDIR --seeds 0,1,2,3,4 [--smoke]
"""
import argparse, csv, gzip, os, random, sys, time
import numpy as np
import torch
import torch.nn.functional as F

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p4_train import HOSTS, REGIONS, Net, load_rows  # same data filter and model as the main experiment

AA = "ACDEFGHIKLMNPQRSTVWY"
PERTS = ["clean", "trunc_50", "trunc_25", "sub_5", "sub_20", "shuffle"]


def perturb(seq, kind, key):
    r = random.Random(key)
    if kind == "clean":
        return seq
    if kind.startswith("trunc_"):
        return seq[: max(10, int(len(seq) * int(kind[6:]) / 100))]
    if kind.startswith("sub_"):
        p = int(kind[4:]) / 100
        return "".join(r.choice(AA) if r.random() < p else c for c in seq)
    if kind == "shuffle":
        l = list(seq); r.shuffle(l); return "".join(l)
    raise ValueError(kind)


def train_ce(rows, idx, labels_h, labels_r, tok, base, dev, seed, steps, lr_enc=2e-5, lr_head=1e-3):
    torch.manual_seed(seed); rng = np.random.default_rng(seed)
    net = Net(base, False, False).to(dev)
    params = [{"params": [p for n, p in net.named_parameters() if n.startswith(("hh", "rh"))], "lr": lr_head},
              {"params": net.enc.parameters(), "lr": lr_enc}]
    opt = torch.optim.AdamW(params, weight_decay=0.01)
    scaler = torch.amp.GradScaler("cuda", enabled=(dev == "cuda"))
    seqs = [rows[i]["seq"] for i in idx]
    for step in range(steps):
        b = rng.choice(len(idx), size=min(32, len(idx)), replace=len(idx) < 32)
        enc = tok([seqs[j] for j in b], padding=True, truncation=True, max_length=602, return_tensors="pt")
        yh = torch.tensor([labels_h[idx[j]] for j in b]).to(dev); yr = torch.tensor([labels_r[idx[j]] for j in b]).to(dev)
        net.train(); lh, lr_ = net(enc["input_ids"].to(dev), enc["attention_mask"].to(dev))
        loss = (F.cross_entropy(lh, yh) + F.cross_entropy(lr_, yr)) / 2
        opt.zero_grad(); scaler.scale(loss).backward(); scaler.unscale_(opt)
        torch.nn.utils.clip_grad_norm_(net.parameters(), 1.0); scaler.step(opt); scaler.update()
    return net


@torch.no_grad()
def predict(net, seqs, tok, dev):
    net.eval(); out = []
    for i in range(0, len(seqs), 64):
        enc = tok(seqs[i:i + 64], padding=True, truncation=True, max_length=602, return_tensors="pt")
        lh, lr_ = net(enc["input_ids"].to(dev), enc["attention_mask"].to(dev))
        ph, pr = F.softmax(lh, -1), F.softmax(lr_, -1)
        for j in range(lh.size(0)):
            out.append((int(ph[j].argmax()), float(ph[j].max()), int(pr[j].argmax()), float(pr[j].max())))
    return out


def write(path, header, rows):
    tmp = path + ".tmp"
    with gzip.open(tmp, "wt", newline="\n") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n"); w.writerow(header); w.writerows(rows)
    os.replace(tmp, path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True); ap.add_argument("--out", required=True); ap.add_argument("--base", required=True)
    ap.add_argument("--seeds", default="0,1,2,3,4"); ap.add_argument("--steps", type=int, default=200)
    ap.add_argument("--folds", type=int, default=5); ap.add_argument("--smoke", action="store_true")
    a = ap.parse_args()
    from transformers import AutoTokenizer
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    rows = load_rows(a.data)
    cl = {}
    with gzip.open(a.data, "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            cl[r["accession"]] = r["cluster"]
    if a.smoke:
        rows = rows[::25]; a.steps = 2
    tok = AutoTokenizer.from_pretrained(a.base)
    yh = [r["y_host"] for r in rows]; yr = [r["y_reg"] for r in rows]
    clusters = sorted({cl[r["acc"]] for r in rows})
    train_idx = [i for i, r in enumerate(rows) if r["split"] == "train"]
    test_idx = [i for i, r in enumerate(rows) if r["split"] != "train"]
    os.makedirs(a.out, exist_ok=True)
    print(f"device {dev}; rows {len(rows)}; clusters {len(clusters)}; train {len(train_idx)}; test {len(test_idx)}; smoke {a.smoke}", flush=True)
    H = ["accession", "split", "cluster", "fold_or_pert", "host_true", "host_top", "host_pmax", "reg_true", "reg_top", "reg_pmax"]
    for seed in map(int, a.seeds.split(",")):
        t0 = time.time()
        # Part A: cluster-held-out CV
        pa = os.path.join(a.out, f"cv__seed{seed}.tsv.gz")
        if os.path.exists(pa):
            print("skip", pa)
        else:
            rs = np.random.default_rng(1000 + seed); order = rs.permutation(len(clusters))
            fold_of = {clusters[j]: k % a.folds for k, j in enumerate(order)}
            res = []
            for k in range(a.folds):
                te = [i for i, r in enumerate(rows) if fold_of[cl[r["acc"]]] == k]; tr = [i for i in range(len(rows)) if i not in set(te)]
                if not te or not tr:
                    continue
                net = train_ce(rows, tr, yh, yr, tok, a.base, dev, seed * 10 + k, a.steps)
                pr = predict(net, [rows[i]["seq"] for i in te], tok, dev)
                for i, (h, hp, g, gp) in zip(te, pr):
                    r = rows[i]; res.append([r["acc"], r["split"], cl[r["acc"]], f"fold{k}", HOSTS[r["y_host"]], HOSTS[h], f"{hp:.4f}", REGIONS[r["y_reg"]], REGIONS[g], f"{gp:.4f}"])
                del net
                if dev == "cuda":
                    torch.cuda.empty_cache()
                print(f"  seed {seed} fold {k} done ({time.time()-t0:.0f}s)", flush=True)
            write(pa, H, res)
        # Part B: perturbations on the original split
        pb = os.path.join(a.out, f"pert__seed{seed}.tsv.gz"); pc = os.path.join(a.out, f"labelshuf__seed{seed}.tsv.gz")
        if os.path.exists(pb) and os.path.exists(pc):
            print("skip", pb, pc)
        else:
            net = train_ce(rows, train_idx, yh, yr, tok, a.base, dev, seed, a.steps)
            res = []
            for kind in PERTS:
                seqs = [perturb(rows[i]["seq"], kind, f"{seed}-{rows[i]['acc']}-{kind}") for i in test_idx]
                for i, (h, hp, g, gp) in zip(test_idx, predict(net, seqs, tok, dev)):
                    r = rows[i]; res.append([r["acc"], r["split"], cl[r["acc"]], kind, HOSTS[r["y_host"]], HOSTS[h], f"{hp:.4f}", REGIONS[r["y_reg"]], REGIONS[g], f"{gp:.4f}"])
            write(pb, H, res); del net
            # Part C: shuffled training labels control
            rp = np.random.default_rng(2000 + seed); perm = rp.permutation(len(train_idx))
            sh, sr = list(yh), list(yr)
            for k, i in enumerate(train_idx):
                sh[i] = yh[train_idx[perm[k]]]; sr[i] = yr[train_idx[perm[k]]]
            net = train_ce(rows, train_idx, sh, sr, tok, a.base, dev, seed + 500, a.steps)
            res = []
            for i, (h, hp, g, gp) in zip(test_idx, predict(net, [rows[i]["seq"] for i in test_idx], tok, dev)):
                r = rows[i]; res.append([r["acc"], r["split"], cl[r["acc"]], "labels_shuffled", HOSTS[r["y_host"]], HOSTS[h], f"{hp:.4f}", REGIONS[r["y_reg"]], REGIONS[g], f"{gp:.4f}"])
            write(pc, H, res); del net
            if dev == "cuda":
                torch.cuda.empty_cache()
            print(f"  seed {seed} perturbation + control done ({time.time()-t0:.0f}s)", flush=True)
    print("P4_ROBUST_DONE", flush=True)


if __name__ == "__main__":
    main()
