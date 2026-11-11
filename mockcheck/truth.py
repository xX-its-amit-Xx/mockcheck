"""Mock community definitions: built-in (from config/expected.tsv) or a user TSV."""
import csv
import os
from dataclasses import dataclass, field

from .parsers import norm, load_aliases

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONFIG = os.path.join(ROOT, "config")

BUILTIN = {
    "zymo-d6300-even": ("even", "ZymoBIOMICS D6300 (even)", "quantitative"),
    "zymo-d6310-log": ("log", "ZymoBIOMICS D6310 (log distribution)", "presence"),
}


@dataclass
class Mock:
    key: str
    title: str
    kind: str                      # quantitative | presence
    dna: dict                      # organism -> expected percent (theoretical genomic DNA)
    cell: dict                     # organism -> expected percent of cells (may be empty)
    species: set                   # organisms known to be in the tube
    anchors: dict = field(default_factory=dict)   # presence mocks: organism -> reliably known percent
    source: str = ""
    organisms: dict = field(default_factory=dict) # microbiology notes by organism


def _tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def load_organisms():
    p = os.path.join(CONFIG, "organisms.tsv")
    return {norm(r["organism"]): r for r in _tsv(p)} if os.path.exists(p) else {}


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def load_contaminants():
    """Genera found in sequenced blank controls (Salter et al. 2014, BMC Biol 12:87, Table 1)."""
    p = os.path.join(CONFIG, "reagent_contaminant_genera.tsv")
    return {r["genus"] for r in _tsv(p)} if os.path.exists(p) else set()


def list_mocks():
    return [(k, v[1], v[2]) for k, v in BUILTIN.items()]


def load_mock(name):
    load_aliases(os.path.join(CONFIG, "aliases.tsv"))
    orgs = load_organisms()
    if name in BUILTIN:
        sample, title, kind = BUILTIN[name]
        exp = os.path.join(CONFIG, "expected.tsv")
        if not os.path.exists(exp):
            raise FileNotFoundError("config/expected.tsv not found; run from a repo checkout or pass --truth FILE")
        rows = _tsv(exp)
        even = [r for r in rows if r["sample"] == "even"]
        mine = [r for r in rows if r["sample"] == sample]
        species = {norm(r["organism"]) for r in even}  # log mock uses the same tube strains
        if kind == "quantitative":
            dna = {norm(r["organism"]): _num(r["expected_percent"]) for r in mine}
            cell = {norm(r["organism"]): _num(r["expected_cell_percent"]) for r in mine}
            return Mock(name, title, kind, dna, cell, species, {}, mine[0]["source"], orgs)
        anchors = {norm(r["organism"]): _num(r["expected_percent"]) for r in mine}
        return Mock(name, title, kind, {}, {}, species, anchors, mine[0]["source"], orgs)
    if os.path.exists(name):
        rows = _tsv(name)
        dna = {norm(r["organism"]): _num(r.get("expected_percent")) for r in rows}
        dna = {k: v for k, v in dna.items() if v is not None}
        cell = {norm(r["organism"]): _num(r.get("expected_cell_percent")) for r in rows}
        cell = {k: v for k, v in cell.items() if v is not None}
        quant = 90 <= sum(dna.values()) <= 110 and len(dna) >= 2
        return Mock(os.path.basename(name), f"custom mock ({os.path.basename(name)})",
                    "quantitative" if quant else "presence", dna, cell, set(dna) or {norm(r["organism"]) for r in rows},
                    {} if quant else dna, "user-supplied truth table", orgs)
    raise ValueError(f"unknown mock '{name}'. Built-in: {', '.join(BUILTIN)}; or pass a TSV with columns organism, expected_percent")
