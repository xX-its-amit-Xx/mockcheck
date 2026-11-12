"""Tie the tool to the benchmark tables and check the format readers. Run: python -m unittest -v"""
import csv, os, subprocess, sys, tempfile, unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from mockcheck.parsers import read_profile, detect_format
from mockcheck.score import analyze
from mockcheck.truth import load_mock

RAW = f"{ROOT}/results/raw"
PATHS = {"kraken2_08gb_c0.0": f"{RAW}/kraken2/even.08gb.c0.0.report.gz", "kraken2_16gb_c0.1": f"{RAW}/kraken2/even.16gb.c0.1.report.gz",
         "bracken_16gb_c0.0": f"{RAW}/bracken/even.16gb.c0.0.bracken", "metaphlan4": f"{RAW}/metaphlan/even.profile.txt"}


class T(unittest.TestCase):
    def test_formats_detected(self):
        self.assertEqual(detect_format(PATHS["kraken2_08gb_c0.0"]), "kraken2")
        self.assertEqual(detect_format(PATHS["bracken_16gb_c0.0"]), "bracken")
        self.assertEqual(detect_format(PATHS["metaphlan4"]), "metaphlan")

    def test_matches_benchmark_tables(self):
        rows = list(csv.DictReader(open(f"{ROOT}/results/even_summary.tsv"), delimiter="\t"))
        mock = load_mock("zymo-d6300-even")
        for method, path in PATHS.items():
            a = analyze(read_profile(path), mock, 0.1)
            r = next(x for x in rows if x["method"] == method and float(x["threshold_pct"]) == 0.1)
            self.assertEqual((a["at_cutoff"]["tp"], a["at_cutoff"]["fp"], a["at_cutoff"]["fn"]), (int(r["tp"]), int(r["fp"]), int(r["fn"])), method)
            r0 = next(x for x in rows if x["method"] == method and float(x["threshold_pct"]) == 0.0)
            self.assertAlmostEqual(a["bray_curtis"], float(r0["bray_curtis_native_basis"]), places=3, msg=method)

    def test_custom_truth_and_generic(self):
        with tempfile.TemporaryDirectory() as d:
            t, p = f"{d}/t.tsv", f"{d}/p.tsv"
            open(t, "w").write("organism\texpected_percent\nEscherichia coli\t50\nListeria monocytogenes\t50\n")
            open(p, "w").write("organism\tpercent\nEscherichia coli\t48\nListeria monocytogenes\t47\nBacillus cereus\t5\n")
            a = analyze(read_profile(p), load_mock(t), 0.1)
            self.assertEqual(a["at_cutoff"]["fp"], 1)
            self.assertLess(a["bray_curtis"], 0.1)

    def test_cli_report_and_demo(self):
        with tempfile.TemporaryDirectory() as d:
            for args in (["report", PATHS["metaphlan4"], "--out", f"{d}/r.md"], ["demo", "--out", f"{d}/demo"]):
                r = subprocess.run([sys.executable, "-m", "mockcheck"] + args, cwd=ROOT, capture_output=True, text=True)
                self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIn("What it does not support", open(f"{d}/r.md", encoding="utf-8").read())

    def test_bad_input_is_a_clear_error(self):
        with tempfile.TemporaryDirectory() as d:
            open(f"{d}/e.txt", "w").write("")
            r = subprocess.run([sys.executable, "-m", "mockcheck", "report", f"{d}/e.txt"], cwd=ROOT, capture_output=True, text=True)
            self.assertNotEqual(r.returncode, 0)
            self.assertNotIn("Traceback", r.stderr)

    def test_spike_report_on_real_phase2_output(self):
        from mockcheck.spike import analyze_spike, render_spike
        rep = os.path.join(ROOT, "results", "raw", "phase2", "kraken2", "spike_P0.16gb.c0.0.report.gz")
        r = analyze_spike(rep, "spike_P0")
        text = render_spike(r, "spike_P0")
        self.assertIn("pseudomonas aeruginosa", text.lower())
        self.assertNotIn("Traceback", text)

    def test_agent_eval_gold_and_grounding(self):
        sys.path.insert(0, os.path.join(ROOT, "agent_eval"))
        import tasks, run_agent_evals as rae
        g = {k: t["gold"]() for k, t in tasks.TASKS.items()}
        self.assertEqual(g["t1_even"], {"found": 7, "false": 1})
        self.assertEqual(g["t3_negctl"]["max_species_reads"], 3)
        self.assertEqual(g["t4_plasma"]["verdict"], "cannot_determine")
        # the renamed genus is the trap: literal reading differs from rename-aware reading
        self.assertNotIn("cutibacterium acnes", g["t4_plasma"]["reagent_genus_flagged"])
        self.assertIn("cutibacterium acnes", g["t4_plasma"]["reagent_flagged_rename_aware"])
        # scorer: right answer passes, wrong answer fails
        self.assertTrue(all(tasks.TASKS["t1_even"]["score"]({"found": 7, "false": 1}, g["t1_even"]).values()))
        self.assertFalse(all(tasks.TASKS["t1_even"]["score"]({"found": 6, "false": 2}, g["t1_even"]).values()))
        # grounding: a number the tools never printed is not grounded
        self.assertEqual(rae.grounded({"n": 4242}, "tool said 17"), 0.0)
        self.assertEqual(rae.grounded({"n": 17}, "tool said 17"), 1.0)
        self.assertEqual(rae.parse_result("x\nRESULT: {\"a\": 1}"), {"a": 1})


if __name__ == "__main__":
    unittest.main()
