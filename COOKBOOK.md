# Cookbook

Copy-paste recipes. Everything runs with Python 3.9 or newer and nothing else installed; the classifiers themselves are only needed if you want to make new outputs (recipe 6).

_Generated from `docs/data/cookbook.json` by `scripts/build_cookbook.py`; the dashboard shows the same recipes._

## 1. Check it works on your machine (30 seconds)

**For:** Anyone. No data, no cluster, no downloads.

```bash
git clone https://github.com/xX-its-amit-Xx/mockcheck && cd mockcheck
python -m unittest discover -s tests
python -m mockcheck demo
```

**You should see:** Tests print OK. The demo scores 72 committed classifier outputs and writes one plain-language report per run to out/demo/reports/ and a table to out/demo/summary.tsv.

**If it goes wrong:** 'No module named mockcheck' means you are not in the repo folder. 'results/raw not found' means you ran it from a copy without the results folder.

## 2. Score one classifier run you already have

**For:** You ran Kraken2, Bracken or MetaPhlAn on a ZymoBIOMICS mock community and want to know whether to trust it.

```bash
python -m mockcheck report my_run.report --mock zymo-d6300-even --out my_run.md
python -m mockcheck report my_run.bracken --mock zymo-d6300-even
python -m mockcheck report my_profile.txt --mock zymo-d6310-log
```

**You should see:** A markdown report with a scorecard, then three lists: what the result supports, what it does not support, what to check next. A last section explains the microbiology behind each surprise.

**If it goes wrong:** The file type is detected automatically. If it guesses wrong, add --format kraken2 | bracken | metaphlan | generic. An empty or wrong file gives a one-line error, not a traceback.

## 3. Choose between tools or settings

**For:** You have several outputs from the same mock and need one ranking.

```bash
python -m mockcheck compare results/raw/kraken2/even.16gb.c0.0.report.gz results/raw/kraken2/even.16gb.c0.1.report.gz results/raw/bracken/even.16gb.c0.0.bracken results/raw/metaphlan/even.profile.txt
```

**You should see:** A table sorted by Bray-Curtis distance to the known answer (lower is better), with how many expected organisms were found and how many false species appeared.

**If it goes wrong:** Compare only runs of the same mock. A tool scored on a different basis (cells versus DNA) is handled for MetaPhlAn automatically.

## 4. Use your own mock community

**For:** Your lab has a different standard, or an in-house spike-in.

```bash
printf 'organism\texpected_percent\nEscherichia coli\t50\nListeria monocytogenes\t50\n' > my_truth.tsv
python -m mockcheck report my_run.report --truth my_truth.tsv
```

**You should see:** Same report, scored against your table. Percentages summing to about 100 give abundance scoring; anything else gives presence scoring only. Add an expected_cell_percent column if the standard publishes cell shares.

**If it goes wrong:** Names must match the database naming after lower-casing. If a tool uses a renamed species, add one line to config/aliases.tsv (alias, canonical).

## 5. Measure a tool's false-call floor with a single-organism control

**For:** You want to know how many wrong species a tool invents when only one organism is present.

```bash
python -m mockcheck report results/raw/kraken2/iso_Listeria_monocytogenes.16gb.c0.0.report.gz --isolate "Listeria monocytogenes"
python -m mockcheck report results/raw/metaphlan/iso_Listeria_monocytogenes.profile.txt --isolate "Listeria monocytogenes"
```

**You should see:** How much of the profile is the true organism, how many other species passed the cutoff, and how much went to close relatives. Anything a tool reports above that level in a real sample is more likely a real signal.

**If it goes wrong:** Use a sequenced pure isolate, not a mix.

## 6. Make the classifier outputs yourself (fresh data, fresh databases)

**For:** You want to reproduce the benchmark from public reads, or run the same tools on your own reads.

```bash
bash run_all.sh                       # on a Slurm cluster: installs tools and databases, streams the reads, runs everything
bash scripts/collect_raw.sh           # shrinks outputs into results/raw (small enough to commit)
bash run_all.sh --local               # re-scores the benchmark, rebuilds the tables and dashboard data
kraken2 --db DB --paired r1.fq.gz r2.fq.gz --report run.report --output -
bracken -d DB -i run.report -o run.bracken -r 150 -l S -t 0
metaphlan r1.fq.gz,r2.fq.gz --input_type fastq --mapout run.bt2out --nproc 8 --db_dir MPA_DB -o run.profile.txt
```

**You should see:** The first three lines reproduce this repository. The last three are the exact tool calls to run on your own reads; their outputs go straight into recipe 2.

**If it goes wrong:** Tool versions and database dates are pinned in env/environment.lock.txt and the README. A different database changes the results: say which one you used.

## 7. A pathogen-sounding name appeared. What now?

**For:** Biosurveillance and biodefense readers, for example a founder or analyst at a company such as Valthos, who need to know whether a species call is a signal or noise.

```bash
python -m mockcheck report my_run.report --mock zymo-d6300-even
python -m mockcheck report results/raw/kraken2/iso_Bacillus_subtilis.16gb.c0.0.report.gz --isolate "Bacillus subtilis"       # a pure Bacillus isolate: what do the tools say?
```

**You should see:** The report lists species from pathogen-containing genera that are not in the tube, and says why a harmless mix can produce them. In the committed Kraken2 run on the harmless even mock, Mycobacterium canetti reached 0.37% and Xanthomonas euvesicatoria 0.12%, so a name alone is not a detection.

**If it goes wrong:** This tool only tells you how reliable the classifier was on a known mix. It is not a threat assessment: a species label from short reads does not give strain, virulence or intent. Corroborate with reads spread across the genome, a second method and a negative control before escalating.

## 8. How rare can a microbe be before my pipeline misses it? (blood-like, host-heavy samples)

**For:** Anyone who classifies samples that are mostly human (blood, plasma, tissue) and needs a detection limit, not just a species list.

```bash
python -m mockcheck spike results/raw/phase2/kraken2/spike_P0.16gb.c0.1.report.gz --sample spike_P0
python -m mockcheck spike results/raw/phase2/kraken2/spike_H0.16gb.c0.1.report.gz --sample spike_H0 --min-reads 3
python scripts/make_spikein.py READS_DIR OUT_DIR 4000000      # build your own: your human reads + your isolate reads at 1e-2 to 1e-5
```

**You should see:** A report in absolute read pairs: which added organisms were detected at which fraction, how many wrong species appeared (the human-only sample is your false-alarm floor), and what each result does not support. Truth tier T2: real reads, constructed mixture.

**If it goes wrong:** The spike-in truth lives in results/spike_manifest.tsv. Species are called at a minimum number of read pairs (default 10); change it with --min-reads. Detection limits here are optimistic for real samples: reads come from clean isolates, with no kit contamination.

## 9. Add another classifier

**For:** You use Centrifuge, Kaiju, mOTUs or similar.

```bash
# 1. write a reader in mockcheck/parsers.py that returns Profile(tool, {species: percent}, path, meta)
# 2. register it in READERS
python -m mockcheck report my_output.tsv --format generic      # works today with any 'organism<TAB>percent' file
```

**You should see:** The generic format needs two columns, organism and percent, and already gets the full report. A proper reader adds the tool-specific warnings (such as Kraken2's reads stuck at genus level).

**If it goes wrong:** Check tests/test_mockcheck.py for how a reader is tested against the committed tables.

## 10. Check how an AI agent handles these reports (with traces)

**For:** Anyone who wants to see what a real agent does with a classifier report, and score its steps and its answer, not just the final number.

```bash
python agent_eval/run_agent_evals.py build                      # makes two throwaway sandboxes (with and without the mockcheck tool) under C:/Temp
python agent_eval/run_agent_evals.py run --models haiku --reps 1   # real Claude Code agents, one per task and condition
python agent_eval/run_agent_evals.py score                        # writes results/agent_eval_runs.tsv and one trace file per run
```

**You should see:** One row per run: did it open the data and run code (L1), is every number in its answer traceable to tool output, does the answer equal what mockcheck's own code computes (L2), and does it decline when the data cannot answer. Each run's steps are in results/agent_traces/.

**If it goes wrong:** Needs the Claude Code CLI installed and logged in (it calls the models and costs a few cents per run). Set MC_AGENT_WORK to change where sandboxes go. If every agent disagrees with gold in the same way, check the task wording and the gold first: that is what happened in our round 1.
