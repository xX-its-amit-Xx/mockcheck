#!/usr/bin/env python3
"""Phase 4: fine-tune a small protein language model on the curated ground truth, with the LOSS TIED TO THE EVAL.

Task: from a viral protein sequence predict (a) host class and (b) world region of collection.
Eval metric (per stratum): selective utility U = mean over samples of  +1 if the call is right, -lam if wrong, 0 if the model abstains.
Loss families (all use the same encoder, data, steps and learning rates; only the objective and the batch mix differ):
  probe          frozen embeddings + linear heads, cross-entropy        (the "vanilla" model)
  ce             full fine-tune, cross-entropy                           (standard supervised fine-tuning)
  util_avg       full fine-tune, maximise expected utility, average over all strata equally
  util_worst     same, but softmin-weighted toward the worst stratum (group-robust)
  util_single:S  same, but ONLY stratum S contributes to the loss (batches drawn only from S)
Expected utility is computed EXACTLY from the softmax over {classes, abstain}: it is the expectation that the policy-gradient
(REINFORCE) estimator in reinforcement learning with verifiable rewards approximates by sampling, without the sampling noise.

Strata = virus x region. Only sequences whose host AND region labels are both tier A (independently verified) are used.
Splits come from the ground-truth file; dev is never used to choose anything for the util_* conditions (no early stopping).

  python scripts/p4_train.py --data data/groundtruth/p4_groundtruth.tsv.gz --out OUTDIR --conds ce,util_avg --seeds 0,1,2 [--steps 300]
"""
import argparse, csv, gzip, math, os, random, sys, time

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

HOSTS = ["avian", "human", "swine", "equine", "canine", "bovine", "feline", "other_mammal", "bat"]
REGIONS = ["africa", "asia", "europe", "north_america", "oceania", "south_america"]
STRATA = [f"influenza_a|{r}" for r in REGIONS] + [f"rotavirus_a|{r}" for r in ("africa", "asia", "oceania")]
LAM = 1.0
MAXLEN = 600


def load_rows(path):
    rows = []
    with gzip.open(path, "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            st = f"{r['virus']}|{r['region']}"
            if r["host_tier"] == "A" and r["geo_tier"] == "A" and r["host_class"] in HOSTS and st in STRATA:
                rows.append(dict(acc=r["accession"], seq=r["protein"][:MAXLEN], virus=r["virus"], stratum=st, split=r["split"],
                                 y_host=HOSTS.index(r["host_class"]), y_reg=REGIONS.index(r["region"]), nn=r["nn_sim_train"]))
    return rows


class Net(nn.Module):
    def __init__(self, base, abstain, freeze):
        super().__init__()
        from transformers import AutoModel
        self.enc = AutoModel.from_pretrained(base)
        if freeze:
            for p in self.enc.parameters():
                p.requires_grad = False
        d = self.enc.config.hidden_size
        self.hh = nn.Linear(d, len(HOSTS) + int(abstain))
        self.rh = nn.Linear(d, len(REGIONS) + int(abstain))

    def embed(self, ids, mask):
        h = self.enc(input_ids=ids, attention_mask=mask).last_hidden_state
        m = mask.unsqueeze(-1).to(h.dtype)
        return (h * m).sum(1) / m.sum(1)

    def forward(self, ids, mask):
        # encoder in fp16 autocast on GPU (about 3x faster on V100); pooling, heads and the loss stay in fp32
        with torch.autocast(device_type="cuda", dtype=torch.float16, enabled=ids.is_cuda):
            h = self.enc(input_ids=ids, attention_mask=mask).last_hidden_state
        h = h.float()
        m = mask.unsqueeze(-1).float()
        e = (h * m).sum(1) / m.sum(1)
        return self.hh(e), self.rh(e)


def exp_util(logits, y):
    """Exact expected utility of the stochastic policy softmax(logits) over {K classes, abstain=last}."""
    p = F.softmax(logits, -1)
    k = logits.size(1) - 1
    py = p.gather(1, y[:, None]).squeeze(1)
    return py - LAM * (1 - py - p[:, k])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--conds", default="ce,util_avg,util_worst"); ap.add_argument("--seeds", default="0,1,2")
    ap.add_argument("--steps", type=int, default=300); ap.add_argument("--base", default="facebook/esm2_t12_35M_UR50D")
    ap.add_argument("--lr_enc", type=float, default=2e-5); ap.add_argument("--lr_head", type=float, default=1e-3)
    ap.add_argument("--tau", type=float, default=0.25); ap.add_argument("--limit", type=int, default=0, help="debug: use only this many rows")
    a = ap.parse_args()
    from transformers import AutoTokenizer
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    rows = load_rows(a.data)
    if a.limit:
        rows = rows[:a.limit]
    tok = AutoTokenizer.from_pretrained(a.base)
    enc = tok([r["seq"] for r in rows], padding=True, truncation=True, max_length=MAXLEN + 2, return_tensors="pt")
    ids, mask = enc["input_ids"], enc["attention_mask"]
    yh = torch.tensor([r["y_host"] for r in rows]); yr = torch.tensor([r["y_reg"] for r in rows])
    train_idx = [i for i, r in enumerate(rows) if r["split"] == "train"]
    by_s = {s: [i for i in train_idx if rows[i]["stratum"] == s] for s in STRATA}
    os.makedirs(f"{a.out}/preds", exist_ok=True)
    print(f"device {dev}; rows {len(rows)}; train {len(train_idx)}; per stratum", {s: len(v) for s, v in by_s.items()}, flush=True)

    conds = []
    for c in a.conds.split(","):
        conds += [f"util_single:{s}" for s in STRATA] if c == "util_single" else [c]
    for cond in conds:
        for seed in map(int, a.seeds.split(",")):
            tag = f"{cond.replace('|', '_').replace(':', '-')}__seed{seed}"
            outp = f"{a.out}/preds/{tag}.tsv.gz"
            if os.path.exists(outp):
                print("skip", tag); continue
            t0 = time.time()
            torch.manual_seed(seed); rng = np.random.default_rng(seed)
            abstain = cond.startswith("util")
            net = Net(a.base, abstain, freeze=(cond == "probe")).to(dev)
            params = [{"params": [p for n, p in net.named_parameters() if n.startswith(("hh", "rh"))], "lr": a.lr_head}]
            if cond != "probe":
                params.append({"params": net.enc.parameters(), "lr": a.lr_enc})
            opt = torch.optim.AdamW(params, weight_decay=0.01)
            scaler = torch.amp.GradScaler("cuda", enabled=(dev == "cuda" and cond != "probe"))
            steps = a.steps if cond != "probe" else 200
            if cond == "probe":
                with torch.no_grad():
                    emb = torch.cat([net.embed(ids[i:i + 64].to(dev), mask[i:i + 64].to(dev)).float().cpu() for i in range(0, len(rows), 64)])
            for step in range(steps):
                if cond == "ce" or cond == "probe":
                    b = rng.choice(train_idx, size=min(32, len(train_idx)), replace=False if len(train_idx) >= 32 else True)
                elif cond.startswith("util_single"):
                    pool = by_s[cond.split(":", 1)[1]]; b = rng.choice(pool, size=32, replace=len(pool) < 32)
                else:  # util_avg / util_worst: equal draw from every stratum
                    b = np.concatenate([rng.choice(by_s[s], size=6, replace=len(by_s[s]) < 6) for s in STRATA if by_s[s]])
                b = torch.as_tensor(b)
                if cond == "probe":
                    lh, lr_ = net.hh(emb[b].to(dev)), net.rh(emb[b].to(dev))
                else:
                    net.train(); lh, lr_ = net(ids[b].to(dev), mask[b].to(dev))
                Yh, Yr = yh[b].to(dev), yr[b].to(dev)
                if cond in ("ce", "probe"):
                    loss = (F.cross_entropy(lh, Yh) + F.cross_entropy(lr_, Yr)) / 2
                else:
                    u = (exp_util(lh, Yh) + exp_util(lr_, Yr)) / 2        # per-sample expected utility
                    if cond.startswith("util_single"):
                        loss = -u.mean()
                    else:
                        sv = torch.tensor([STRATA.index(rows[i]["stratum"]) for i in b.tolist()], device=dev)
                        us = torch.stack([u[sv == k].mean() for k in range(len(STRATA)) if (sv == k).any()])
                        if cond == "util_avg":
                            loss = -us.mean()
                        else:   # util_worst: softmin weights put most weight on the lowest-utility strata
                            w = F.softmax(-us.detach() / a.tau, 0)
                            loss = -(w * us).sum()
                opt.zero_grad(); scaler.scale(loss).backward(); scaler.unscale_(opt)
                torch.nn.utils.clip_grad_norm_(net.parameters(), 1.0); scaler.step(opt); scaler.update()
                if step % 100 == 0:
                    print(f"  {tag} step {step} loss {loss.item():.4f}", flush=True)
            # predictions for every row (all splits), deterministic policy = argmax action
            net.eval(); out = []
            with torch.no_grad():
                for i in range(0, len(rows), 64):
                    if cond == "probe":
                        lh, lr_ = net.hh(emb[i:i + 64].to(dev)), net.rh(emb[i:i + 64].to(dev))
                    else:
                        lh, lr_ = net(ids[i:i + 64].to(dev), mask[i:i + 64].to(dev))
                    for j in range(lh.size(0)):
                        ph, pr = F.softmax(lh[j], -1), F.softmax(lr_[j], -1)
                        out.append((int(ph.argmax()), float(ph[:len(HOSTS)].max()), int(pr.argmax()), float(pr[:len(REGIONS)].max()), int(ph[:len(HOSTS)].argmax()), int(pr[:len(REGIONS)].argmax())))
            with gzip.open(outp, "wt", newline="\n") as f:
                w = csv.writer(f, delimiter="\t", lineterminator="\n")
                w.writerow(["accession", "split", "stratum", "virus", "nn_sim_train", "host_true", "host_act", "host_pmax", "host_top", "reg_true", "reg_act", "reg_pmax", "reg_top", "abstain_idx_host", "abstain_idx_reg"])
                for r, (ha, hp, ra, rp, ht, rt) in zip(rows, out):
                    w.writerow([r["acc"], r["split"], r["stratum"], r["virus"], r["nn"], HOSTS[r["y_host"]], ha, f"{hp:.4f}", ht, REGIONS[r["y_reg"]], ra, f"{rp:.4f}", rt, len(HOSTS), len(REGIONS)])
            print(f"done {tag} in {time.time()-t0:.0f}s", flush=True)
            del net, opt
            if dev == "cuda":
                torch.cuda.empty_cache()


if __name__ == "__main__":
    main()
