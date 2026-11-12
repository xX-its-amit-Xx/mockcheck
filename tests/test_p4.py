"""Tests for the Phase 4 ground-truth set: the evidence rule, the strain-name parsers, and the holdout invariants."""
import collections, csv, gzip, os, sys, unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
sys.argv = [sys.argv[0]]          # the build script reads sys.argv[1] at import time
import p4_build_groundtruth as G  # noqa: E402

GT = os.path.join(ROOT, "data", "groundtruth", "p4_groundtruth.tsv.gz")


def rows():
    with gzip.open(GT, "rt") as f:
        return list(csv.DictReader(f, delimiter="\t"))


class EvidenceRule(unittest.TestCase):
    def test_tier_logic(self):
        self.assertEqual(G.tier("swine", {"swine"}), "A")       # independent evidence agrees
        self.assertEqual(G.tier("swine", set()), "B")           # no independent evidence: weak
        self.assertEqual(G.tier("swine", {"human"}), "X")       # conflict: label is dropped
        self.assertEqual(G.tier("swine", {"swine", "human"}), "B")   # ambiguous evidence is not promoted

    def test_influenza_strain_name(self):
        rec = dict(virus="influenza_a", strain="A/swine/Iowa/A02524810/2020(H1N1)", isolate="", definition="")
        self.assertEqual(G.evidence(rec), ({"swine"}, {"usa"}, {2020}))
        human = dict(virus="influenza_a", strain="A/Idaho/24/2018", isolate="", definition="")
        self.assertEqual(G.evidence(human), ({"human"}, {"usa"}, {2018}))

    def test_rotavirus_structured_name(self):
        rec = dict(virus="rotavirus_a", strain="RVA/Cow-tc/NGA/NGRBg8/1998/G8P[1]", isolate="", definition="")
        self.assertEqual(G.evidence(rec), ({"bovine"}, {"nigeria"}, {1998}))

    def test_host_genus_groups_common_and_latin_names(self):
        self.assertEqual(G.host_genus("canine"), G.host_genus("Canis lupus familiaris"))
        self.assertEqual(G.host_genus("porcine"), G.host_genus("Sus scrofa"))


@unittest.skipUnless(os.path.exists(GT), "ground-truth file not present")
class GroundTruthInvariants(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r = rows()

    def test_unique_accessions_and_known_splits(self):
        self.assertEqual(len({x["accession"] for x in self.r}), len(self.r))
        self.assertTrue({x["split"] for x in self.r} <= {"train", "dev", "test_geo", "test_host", "test_time", "blind"})

    def test_gold_labels_have_agreeing_evidence(self):
        for x in self.r:
            if x["host_tier"] == "A":
                self.assertEqual(x["host_evidence"], x["host_class"], x["accession"])
            if x["geo_tier"] == "A":
                self.assertEqual(x["geo_evidence"], x["country"], x["accession"])

    def test_no_country_or_host_genus_leaks_from_holdout_into_train(self):
        train = [x for x in self.r if x["split"] == "train"]
        geo = {x["country"] for x in self.r if x["split"] == "test_geo"}
        self.assertFalse(geo & {x["country"] for x in train}, "a held-out country appears in training")
        host = {G.host_genus(x["host_species_raw"]) for x in self.r if x["split"] == "test_host"}
        self.assertFalse(host & {G.host_genus(x["host_species_raw"]) for x in train}, "a held-out host genus appears in training")

    def test_blind_clusters_are_sealed(self):
        blind = {x["cluster"] for x in self.r if x["split"] == "blind"}
        self.assertFalse(blind & {x["cluster"] for x in self.r if x["split"] == "train"}, "a blind cluster has training members")

    def test_no_select_agent_viruses(self):
        self.assertTrue({x["virus"] for x in self.r} <= {"influenza_a", "rabies", "rotavirus_a"})


if __name__ == "__main__":
    unittest.main()
