"""Read classifier outputs into one shape: species -> percent (renormalised to 100)."""
import gzip
import os
from dataclasses import dataclass, field

# Names that differ between databases but mean the same thing for scoring.
BUILTIN_ALIASES = {
    "limosilactobacillus fermentum": "lactobacillus fermentum",
    "limosilactobacillus": "lactobacillus",
}
_aliases = dict(BUILTIN_ALIASES)


def load_aliases(path):
    """Add aliases from a two-column TSV (alias, canonical)."""
    if path and os.path.exists(path):
        with open(path) as f:
            next(f, None)
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 2:
                    _aliases[p[0].strip().lower()] = p[1].strip().lower()


def norm(name):
    n = name.replace("_", " ").strip().lower()
    return _aliases.get(n, n)


def genus(name):
    g = name.split()[0] if name.split() else name
    return _aliases.get(g, g)


@dataclass
class Profile:
    tool: str                      # kraken2 | bracken | metaphlan4 | generic
    species: dict                  # normalised name -> percent of species-level profile
    path: str = ""
    meta: dict = field(default_factory=dict)  # tool-specific facts used by the report


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def _renorm(d):
    s = sum(d.values())
    return {k: 100.0 * v / s for k, v in d.items()} if s > 0 else {}


def detect_format(path):
    with _open(path) as f:
        head = [next(f, "") for _ in range(40)]
    text = "".join(head)
    if "fraction_total_reads" in text:
        return "bracken"
    if "#mpa_" in text or "s__" in text or "relative_abundance" in text:
        return "metaphlan"
    for line in head:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 6:
            try:
                float(p[0]); int(p[1]); int(p[2])
                return "kraken2"
            except ValueError:
                continue
    if all(len(l.split("\t")) >= 2 for l in head if l.strip()):
        return "generic"
    raise ValueError(f"cannot recognise the format of {path}; pass --format")


def read_kraken(path):
    sp, total, unclass = {}, 0, 0
    with _open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            reads, rank, name = int(p[1]), p[3], p[5].strip()
            if rank == "U":
                unclass = reads
            elif rank == "R" and name == "root":
                total = reads
            elif rank == "S":
                sp[norm(name)] = sp.get(norm(name), 0) + reads
    allr = total + unclass
    s = sum(sp.values())
    meta = dict(total_reads=allr,
                unclassified_frac=(unclass / allr if allr else None),
                species_level_frac=(s / allr if allr else None))
    return Profile("kraken2", _renorm(sp), path, meta)


def read_kraken_counts(path):
    """Kraken2 report -> (species {name: clade read pairs}, genera {name: clade reads}, totals)."""
    sp, ge = {}, {}
    tot = dict(root=0, unclassified=0)
    with _open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            reads, rank, name = int(p[1]), p[3], p[5].strip()
            if rank == "U":
                tot["unclassified"] = reads
            elif rank == "R" and name == "root":
                tot["root"] = reads
            elif rank == "D":
                tot[name.lower()] = reads
            elif rank == "S":
                sp[norm(name)] = sp.get(norm(name), 0) + reads
            elif rank == "G":
                ge[norm(name)] = ge.get(norm(name), 0) + reads
    tot["total"] = tot["root"] + tot["unclassified"]
    return sp, ge, tot


def read_bracken(path):
    sp, frac = {}, 0.0
    with _open(path) as f:
        header = next(f).rstrip("\n").split("\t")
        ni, fi = header.index("name"), header.index("fraction_total_reads")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) <= max(ni, fi):
                continue
            n = norm(p[ni])
            sp[n] = sp.get(n, 0) + float(p[fi])
            frac += float(p[fi])
    return Profile("bracken", _renorm(sp), path, dict(species_level_frac=frac))


def read_metaphlan(path):
    sp = {}
    with _open(path) as f:
        for line in f:
            if line.startswith("#") or "s__" not in line:
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            last = p[0].split("|")[-1]
            if last.startswith("s__") and "t__" not in p[0]:
                sp[norm(last[3:])] = float(p[2])
    return Profile("metaphlan4", _renorm(sp), path, {})


def read_generic(path):
    """Two columns: organism <TAB> percent (header optional)."""
    sp = {}
    with _open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 2:
                continue
            try:
                sp[norm(p[0])] = sp.get(norm(p[0]), 0) + float(p[1])
            except ValueError:
                continue
    return Profile("generic", _renorm(sp), path, {})


READERS = {"kraken2": read_kraken, "bracken": read_bracken, "metaphlan": read_metaphlan,
           "metaphlan4": read_metaphlan, "generic": read_generic}


def read_profile(path, fmt="auto"):
    fmt = detect_format(path) if fmt == "auto" else fmt
    if fmt not in READERS:
        raise ValueError(f"unknown format {fmt}; choose from {sorted(READERS)}")
    prof = READERS[fmt](path)
    if not prof.species:
        raise ValueError(f"{path}: no species-level rows found (is this a {fmt} file, and is it empty?)")
    return prof
