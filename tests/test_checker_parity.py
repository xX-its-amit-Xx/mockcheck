"""The in-browser checker (docs/checker.js) must give the same numbers as the Python scorer."""
import glob
import json
import math
import os
import shutil
import subprocess
import unittest

from mockcheck import truth
from mockcheck.parsers import read_profile
from mockcheck.score import analyze

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@unittest.skipUnless(shutil.which("node"), "node not installed")
class CheckerParity(unittest.TestCase):
    def test_matches_python(self):
        files = sorted(glob.glob(f"{ROOT}/results/raw/*/even.*") + glob.glob(f"{ROOT}/results/raw/*/log.*"))
        self.assertTrue(files)
        res = subprocess.run(["node", f"{ROOT}/tests/checker_parity.js", f"{ROOT}/docs/data/checker_truth.json", *files],
                             capture_output=True, text=True, check=True)
        js = json.loads(res.stdout)
        for f in files:
            key = "zymo-d6300-even" if os.path.basename(f).startswith("even") else "zymo-d6310-log"
            prof = read_profile(f)
            py = analyze(prof, truth.load_mock(key), 0.1)
            j = js[f]
            with self.subTest(file=os.path.relpath(f, ROOT)):
                self.assertEqual(j["tool"], prof.tool)
                for k in ("tp", "fp", "fn"):
                    self.assertEqual(j[k], py["at_cutoff"][k])
                self.assertAlmostEqual(j["fpmass"], py["false_positive_mass_pct"], places=6)
                self.assertEqual(j["nrep"], py["n_species_reported"])
                if "bray_curtis" in py:
                    self.assertTrue(math.isclose(j["bray"], py["bray_curtis"], abs_tol=1e-9))


if __name__ == "__main__":
    unittest.main()
