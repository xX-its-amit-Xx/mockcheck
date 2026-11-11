#!/usr/bin/env python3
"""Build constructed-truth spike-in samples from REAL reads.

Background: real human whole-genome reads (blood-derived NA12878). Spiked in: real reads of the ten isolates
at fixed fractions of the total read pairs. Four panels rotate the organisms through four levels so that every
organism is seen at every level exactly once; a human-only sample is the negative control.
Truth is known by construction and written to results/spike_manifest.tsv.

usage: make_spikein.py READS_DIR OUT_DIR [TOTAL_PAIRS]
"""
import gzip, os, sys

ORGS = ["Pseudomonas_aeruginosa", "Escherichia_coli", "Salmonella_enterica", "Lactobacillus_fermentum",
        "Enterococcus_faecalis", "Staphylococcus_aureus", "Listeria_monocytogenes", "Bacillus_subtilis",
        "Saccharomyces_cerevisiae", "Cryptococcus_neoformans"]
LEVELS = [1e-2, 1e-3, 1e-4, 1e-5]
reads_dir, out_dir = sys.argv[1], sys.argv[2]
TOTAL = int(sys.argv[3]) if len(sys.argv) > 3 else 4_000_000
here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def records(path):
    with gzip.open(path, "rt") as f:
        while True:
            r = [f.readline() for _ in range(4)]
            if not r[0]:
                return
            yield r


def load_pairs(prefix):
    return list(zip(records(f"{prefix}_1.fastq.gz"), records(f"{prefix}_2.fastq.gz")))


def evenly(pairs, n):
    """Deterministic, evenly spaced subset (no random numbers, fully reproducible)."""
    if n > len(pairs):
        raise SystemExit(f"need {n} pairs but only {len(pairs)} available")
    step = len(pairs) / n
    return [pairs[int(i * step)] for i in range(n)]


iso = {o: load_pairs(f"{reads_dir}/iso_{o}") for o in ORGS}
os.makedirs(out_dir, exist_ok=True)
manifest = ["sample\torganism\tfraction\tspiked_pairs\ttotal_pairs"]
panels = [(f"spike_P{k}", {o: LEVELS[(i + k) % 4] for i, o in enumerate(ORGS)}) for k in range(4)] + [("spike_H0", {})]
for name, assign in panels:
    spikes = {o: round(f * TOTAL) for o, f in assign.items()}
    n_h = TOTAL - sum(spikes.values())
    # human reads: the first n_h pairs of the human file, streamed (no need to hold 4M pairs in memory)
    outs = [gzip.open(f"{out_dir}/{name}_{i}.fastq.gz", "wt", compresslevel=1) for i in (1, 2)]
    n = 0
    for a, b in zip(records(f"{reads_dir}/human_1.fastq.gz"), records(f"{reads_dir}/human_2.fastq.gz")):
        outs[0].writelines(a); outs[1].writelines(b); n += 1
        if n >= n_h:
            break
    if n < n_h:
        raise SystemExit(f"human reads short: {n} < {n_h}")
    for o, k in spikes.items():
        for a, b in evenly(iso[o], k):
            outs[0].writelines(a); outs[1].writelines(b)
        manifest.append(f"{name}\t{o.replace('_', ' ')}\t{assign[o]:g}\t{k}\t{TOTAL}")
    if not spikes:
        manifest.append(f"{name}\t(none)\t0\t0\t{TOTAL}")
    for h in outs:
        h.close()
    print(f"built {name}: {n_h} human + {sum(spikes.values())} spiked pairs", flush=True)
open(f"{here}/results/spike_manifest.tsv", "w", newline="\n").write("\n".join(manifest) + "\n")
