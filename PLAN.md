# mockcheck: plan

Goal: a reproducible benchmark of taxonomic read classifiers on public mock
microbial communities, ending in a plain-language decision report.

## Question the report answers
"I have shotgun reads from a community and need to pick a classifier and
settings. Which one should I use, how far can I trust its output, and what
should I check before acting on it?"

## Data (public, harmless, accessions recorded in `config/samples.tsv`)
Source study: Nicholls et al. 2019, GigaScience 8(5):giz043 (ENA PRJEB29504).
| Run | What | Role |
|---|---|---|
| ERR2984773 | ZymoBIOMICS D6300 "even" mock, Illumina MiSeq 2x151 | main benchmark, quantitative truth |
| ERR2935805 | ZymoBIOMICS D6310 "log" mock, HiSeq 2x101 | presence/absence only (the paper text gives only 2 of 10 expected values) |
| ERR2935848 .. ERR2935857 | the 10 single-isolate runs, HiSeq 2x101 | negative controls: truth is one organism, so every other call is a false positive |

Expected composition of D6300 (from the paper): 8 bacteria at 12% each
(B. subtilis, L. monocytogenes, E. faecalis, S. aureus, S. enterica, E. coli,
P. aeruginosa, L. fermentum), 2 yeasts at 2% each (S. cerevisiae,
C. neoformans). The paper does not say whether this is by cells or DNA; the
report treats that as a stated uncertainty (it matters for read-level truth,
because genome size and copy number shift read share).

Subsampling: only the first N read pairs of each run are streamed from ENA
(N = 2,000,000 for the even mock, 500,000 for the log mock, 300,000 per
isolate). No full FASTQ is ever downloaded and no reads are committed.

## Classifiers and settings
| Id | Tool | Database (pinned) | Notes |
|---|---|---|---|
| k2_std8 | Kraken2 | standard 8 GB capped, 2025-04-02 | confidence 0.0 and 0.1 |
| k2_std16 | Kraken2 | standard 16 GB capped, 2025-04-02 | confidence 0.0 and 0.1 |
| bracken_std16 | Bracken on k2_std16 (confidence 0.0) | same | read-length-aware abundance re-estimation |
| mpa4 | MetaPhlAn 4 | latest mpa_vJan25 index, pinned in `env/` | marker-based, reports relative abundance of classified cells |
Tool versions are pinned in `env/environment.yml`.

## Metrics (species level, genus as secondary)
- Presence: precision, recall, F1 at abundance thresholds 0, 0.01%, 0.1%, 1%.
- Abundance: Bray-Curtis and L1 distance to expected; per-organism error.
- Fraction of reads/abundance assigned to organisms not in the mock (false
  positive mass), and fraction unclassified.
- Negative controls: false-positive species count and mass per isolate.
- Runtime and peak memory per classifier.
Taxonomy names are harmonised to NCBI species by a small, committed,
reviewable alias table (`config/aliases.tsv`); anything not mapped is
reported, not silently dropped.

## Layout
```
config/      samples.tsv, expected.tsv (with source), aliases.tsv
env/         environment.yml (pinned)
scripts/     fetch_reads.sh, setup_env_and_dbs.sh, run_classifiers.sbatch,
             parse_*.py, evaluate.py, make_report.py
results/     small committed tables (tsv) and figures (png)
REPORT.md    the decision report
run_all.sh   one-command reproduction (Explorer/Slurm), with a --local path
             for the parsing/evaluation step from committed raw tables
```

## Compute and storage
- Reads, databases, scratch: `/scratch/$USER/mockcheck` on Explorer
  (set `MOCKCHECK_SCRATCH`). Nothing large on D: or in git.
- Slurm: public `short` partition only, one job at a time, modest memory
  (<= 64 GB, 8 cores), no GPUs. Submitted as `$USER`.
- Small raw classifier outputs (reports only, no per-read files) are
  committed under `results/raw/`; cold large outputs go to OneDrive offload.

## Steps
1. Plan, repo skeleton, config tables (this commit).
2. Fetch reads (streamed subsample) and build env + databases on Explorer.
3. Run classifiers; collect reports.
4. Parse, harmonise, evaluate; commit result tables and figures.
5. Write REPORT.md: for each run, what the result supports, what it does
   not, and what to check next. Write README with one-command reproduction.

## Report layer (the point of the project)
Each section of REPORT.md ends with three fixed blocks:
**Supports**, **Does not support**, **Check next**, in plain language, and the
report opens with a one-paragraph recommendation by use case.

## Known risks
- Capped databases can miss genomes: a miss is a finding, not a bug.
- Species naming differences between tools (e.g. Bacillus/Priestia, Listeria
  synonyms) can look like errors; aliases are explicit and listed.
- First-N-reads subsampling assumes reads are not ordered by taxon; checked
  by comparing a second window for the even mock.
