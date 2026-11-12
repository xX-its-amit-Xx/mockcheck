# mockcheck

**Reproducible benchmark of taxonomic read classifiers on public mock microbial communities, with a plain-language decision report.**

> **Start here (no biology needed): [the live dashboard](https://mockcheck.pages.dev)**
> A guided story with clickable terms, live charts, a toy classifier you can play with, and a dev log of how this was built, including what broke.
> **Then read [REPORT.md](REPORT.md)**: the decision report.

## The question

You have DNA from a mixed sample (gut, soil, a probiotic, a clinical swab) and a computer program says what organisms are in it. **Can you trust it?** Mock communities are tubes mixed from known organisms in known amounts, so we can check programs against a known answer.

## The answer in one table

| If you need... | Use | Watch out for |
|---|---|---|
| Which species are present and roughly how much | **MetaPhlAn 4**: closest to truth, 1 false species on the even mock | Slowest here (548 s, 34 GB RAM); misses organisms below about 0.2% of cells |
| Fast screening or read-level labels | **Kraken2 + Bracken** | Over 150 false low-level species; cannot see fungi at all |
| Abundances straight from Kraken2 | **Don't** | A 12% organism read as 26% to 69% |
| A shrunk 8 GB Kraken2 database | Only if memory forces it | 73% of reads left unclassified at confidence 0.1 |

The biggest swings came from **what the database contains and what names it uses**, not from algorithm cleverness. Details, numbers and caveats: [REPORT.md](REPORT.md). Each run ends with three blocks: what the result **supports**, what it **does not**, and what to **check next**.

## Use it on your own run (the tool)

`mockcheck` reads a Kraken2, Bracken or MetaPhlAn output (or any two-column organism/percent file), scores it against a mock community with a known answer, and writes the report: **what it supports, what it does not, what to check next**, plus the microbiology behind each surprise. Python 3.9+, standard library only, no install step.

```bash
git clone https://github.com/xX-its-amit-Xx/mockcheck && cd mockcheck
python -m mockcheck demo                                    # scores all 72 committed outputs, writes out/demo/reports/*.md
python -m mockcheck report my_run.report --mock zymo-d6300-even
python -m mockcheck compare run_a.report run_b.bracken run_c.profile.txt
python -m mockcheck report my_isolate.report --isolate "Listeria monocytogenes"   # false-call floor
```

All recipes, with what you should see and what to do if it goes wrong, are in **[COOKBOOK.md](COOKBOOK.md)** (also on the dashboard, with a report explorer for every committed run). The tool's numbers match the benchmark tables exactly; `python -m unittest discover -s tests` checks that.

**Why microbiology is in the report:** a species name from short reads is a lead, not an identification. `config/organisms.tsv` records, for each mock organism, its cell wall, genome size, look-alike relatives and why that matters for detection (for example: *E. coli* and *Shigella* are near-identical by k-mers; the *B. cereus* group, which includes the anthrax agent, cannot be told apart reliably from short reads; yeast has a tough wall and is invisible to bacteria-only databases). The report also flags species from pathogen-containing genera that appear in a harmless mix: in the committed Kraken2 run, *Mycobacterium canetti* reached 0.37% of the profile in a tube that contains no mycobacteria. That is the practical lesson for biosurveillance and biodefense teams (for example a company such as Valthos, which builds software to detect biological threats): classifier output needs corroboration and a negative control before anyone escalates. This repository is a benchmark of general-purpose tools on harmless standards, not a threat-detection system, and it makes no claim about any company's methods.

## What is in the box

| Path | What |
|---|---|
| `REPORT.md` | The decision report |
| `docs/` | The dashboard (a static site, deployed on Cloudflare Pages) |
| `results/*.tsv`, `results/fig/` | Small result tables and figures |
| `results/raw/` | Small raw classifier reports (no reads are in the repo) |
| `config/` | Samples and accessions, expected composition **with its source**, name aliases |
| `env/environment.lock.txt` | Exact tool versions |
| `scripts/` | Fetch, setup, run, collect, evaluate, build site |
| `PLAN.md` | The plan written before any code |

## Phase 2: blood-like samples

Real human reads with real microbe reads spiked in at 1e-2 to 1e-5, plus 16 public sepsis plasma runs. Summary and limits are in [REPORT.md](REPORT.md#phase-2-real-reads-constructed-truth-host-heavy-samples) and on the dashboard ("Real data" section). Score your own: `python -m mockcheck spike results/raw/phase2/kraken2/spike_P0.16gb.c0.1.report.gz --sample spike_P0`.

## Phase 3: evaluation frameworks and agent traces

Decisions on geography and host, what other groups do, public evaluations, and a harness that runs real Claude Code agents on these tasks and records every step. See [REPORT.md](REPORT.md#phase-3-how-to-evaluate-who-else-does-and-how-agents-behave) and the dashboard ("Phase 3"). Run it: `python agent_eval/run_agent_evals.py build`, then `run --models haiku --reps 1`, then `score` (needs the Claude Code CLI; a few cents per run). Traces are in `results/agent_traces/`. Source inventory: `python scripts/crawl_sources.py`.

## Data (public and harmless)

Nicholls et al. 2019, *GigaScience* 8(5):giz043, ENA project PRJEB29504. The even mock **ERR2984773** (ZymoBIOMICS D6300, MiSeq 2x151), the log mock **ERR2935805** (D6310, HiSeq 2x101), and ten single-organism runs **ERR2935848 to ERR2935857** as negative controls. Only the first N read pairs of each run are streamed (2,000,000 / 500,000 / 300,000); nothing large is stored in git.

The expected composition comes from Zymo's D6300 datasheet (12% and 2% are theoretical genomic DNA; the cell-number column is also recorded). Details and caveats in `config/expected.tsv` and the report.

## Tools and versions

Kraken2 2.17.2 (standard 8 GB and 16 GB indexes, 2025-04-02), Bracken 3.1, MetaPhlAn 4.2.6 (index `mpa_vJan26_CHOCOPhlAnSGB_202605`), Bowtie2 2.5.5. Full list: `env/environment.lock.txt`.

## Reproduce

```bash
# On a Slurm cluster (roughly 24 GB of disk for databases; the heaviest step is MetaPhlAn at about 9 minutes per 2M-pair sample):
bash run_all.sh
# when the classifier job finishes:
bash scripts/collect_raw.sh && bash run_all.sh --local

# No cluster? Re-score the committed raw reports on any machine (python3, pandas, matplotlib):
bash run_all.sh --local
```

Big files (reads, databases, scratch) live under `MOCKCHECK_SCRATCH` (default `/scratch/$USER/mockcheck`) and never in git. Heavy runs use the cluster's public `short` partition with a modest footprint (one 8-core job at a time).

## Honest limits

Two mocks, one lab, public data; real samples will differ. The even mock's "truth" is the vendor's theoretical recipe. Some background chapters on the dashboard were drafted by AI research agents, spot-checked by a second agent, and anything unverified is listed on the page rather than stated as fact. The full list of limitations is at the end of [REPORT.md](REPORT.md).

## Deploying the dashboard

`docs/` is a plain static site (no build step). On Cloudflare Pages: `wrangler pages deploy docs --project-name mockcheck`.

## License

See [LICENSE](LICENSE).
