"""Turn numbers into the report layer: what the run supports, what it does not, what to check next,
plus the microbiology behind each surprise. Rule-based and deterministic: no model is called."""
from .parsers import genus

EUK_GENERA = {"saccharomyces", "cryptococcus", "candida", "aspergillus", "malassezia", "penicillium", "fusarium"}
KRAKEN_LIKE = ("kraken2", "bracken")
# Genera that contain recognised human, animal or agricultural pathogens. A call here in a mock that
# contains none of them is a false alarm by construction; the report says how to treat such calls.
WATCH_GENERA = {"mycobacterium", "yersinia", "brucella", "francisella", "burkholderia", "vibrio", "shigella",
                "clostridium", "clostridioides", "bordetella", "coxiella", "legionella", "neisseria", "corynebacterium",
                "streptococcus", "campylobacter", "klebsiella", "acinetobacter", "xanthomonas"}


def pct(x, d=2):
    return f"{x:.{d}f}%"


def title_case(n):
    return n[:1].upper() + n[1:]


def _micro(mock, org):
    return mock.organisms.get(org)


def build_findings(a, mock):
    S, N, C, M = [], [], [], []   # supports, does not support, check next, microbiology
    tool, cut, c = a["tool"], a["cutoff"], a["at_cutoff"]
    truth_n = len(mock.species)
    meta = a["meta"] or {}

    # ---- presence
    if a["detected"]:
        S.append(f"{c['tp']} of {truth_n} organisms known to be in the tube were reported at or above {pct(cut)}"
                 f" ({', '.join(title_case(o) for o in a['detected'][:truth_n])}).")
    if c["fp"] <= 1:
        S.append(f"Few extra species: {c['fp']} organism(s) not in the tube passed the {pct(cut)} cutoff (precision {c['precision']:.2f}).")
    # ---- composition
    if "bray_curtis" in a:
        bc = a["bray_curtis"]
        basis = "cells" if a["basis"] == "cell" else "genomic DNA"
        if bc <= 0.25:
            S.append(f"The abundance profile is close to the known answer (Bray-Curtis {bc:.2f} against {basis}; 0 is perfect, 1 is no overlap).")
        elif bc <= 0.5:
            N.append(f"Do not read the percentages as quantities: the profile only roughly matches (Bray-Curtis {bc:.2f} against {basis}). Rank order may still be useful; check which organisms are far off below.")
        else:
            N.append(f"Do not use these percentages: the profile is far from the known answer (Bray-Curtis {bc:.2f} against {basis}).")
        if a["within_2x"]:
            S.append(f"{len(a['within_2x'])} organisms were within two-fold of the expected share ({', '.join(title_case(o) for o in a['within_2x'])}).")
        for r in a["off_by_more_than_2x"][:5]:
            N.append(f"{title_case(r['organism'])} was reported at {pct(r['observed'])} versus {pct(r['expected'])} expected ({r['fold']:.1f}x): "
                     f"do not trust this organism's abundance from this run.")
    # ---- anchors (presence mocks)
    for an in a.get("anchors", []):
        ratio = an["observed"] / an["expected"] if an["expected"] else None
        if an["expected"] >= 1:
            if ratio and 0.8 <= ratio <= 1.2:
                S.append(f"The dominant organism ({title_case(an['organism'])}) was quantified well: {pct(an['observed'])} against {pct(an['expected'])} expected.")
            else:
                N.append(f"The dominant organism ({title_case(an['organism'])}) came out at {pct(an['observed'])} against {pct(an['expected'])} expected.")
        else:
            tr = meta.get("total_reads")
            exp_reads = f" That is about {an['expected'] / 100 * tr:.1f} expected reads in this sample." if tr else ""
            N.append(f"{title_case(an['organism'])} is expected at {an['expected']}%, so the run cannot show whether the method sees organisms that rare.{exp_reads}")
    # ---- missed organisms
    for org in a["missed"]:
        info = _micro(mock, org)
        N.append(f"{title_case(org)} was not reported at {pct(cut)} or above. Its absence from the output is not evidence it is absent.")
        if info:
            M.append(f"**{title_case(org)}** ({info['gram_or_wall']}, about {info['approx_genome_mb']} Mb): {info['why_it_matters_for_detection']}")
        if genus(org) in EUK_GENERA and tool in KRAKEN_LIKE:
            C.append(f"{title_case(org)} is a eukaryote. If this run used the Kraken2 standard database, that database has no fungi, so check what the reads were forced onto (rerun with a database that includes fungi, or use a marker-gene tool).")
    if a["spizizenii_note"]:
        N.append("Strict scoring counts Bacillus subtilis as missed, but the tool reported Bacillus spizizenii, a former subspecies of B. subtilis that databases now treat as its own species. This is a naming difference, not a missed organism.")
        C.append("Check the species names your database uses (taxonomy renames are common) and add an alias line in config/aliases.tsv before judging recall.")
    # ---- false positives
    if a["false_positives"]:
        top = ", ".join(f"{title_case(k)} {pct(v)}" for k, v in a["false_positives"][:4])
        N.append(f"The extra species are not evidence those organisms are present. Largest: {top}. At cutoff 0% this run reported {a['n_species_reported']} species in total.")
        known_genera = {genus(o) for o in mock.species}
        for k, v in a["false_positives"][:8]:
            if genus(k) in known_genera:
                rel = next(o for o in mock.species if genus(o) == genus(k))
                M.append(f"**{title_case(k)}** shares a genus with {title_case(rel)}: a close-relative false call, typical when reads match conserved regions of related genomes.")
    watch = [(k, v) for k, v in a["false_positives"] if genus(k) in WATCH_GENERA and genus(k) not in {genus(o) for o in mock.species}]
    if watch:
        names = ", ".join(f"{title_case(k)} {pct(v)}" for k, v in watch[:4])
        N.append(f"Species from pathogen-containing genera were reported ({names}) although nothing like them is in the tube. A species name on its own is not a detection of a threat.")
        C.append("Before escalating any call from a pathogen-containing genus, require corroboration: many reads spread across the genome (not a few conserved regions), a second tool or a marker-gene check, and a negative control run through the same database.")
        M.append("**Why a harmless mix produces pathogen-genus calls:** short reads from conserved genes (ribosomal and housekeeping genes) match many genomes in a database, and species in large genera such as Mycobacterium or Burkholderia share long stretches of sequence. A likely cause (not verified for this specific run): low-confidence k-mer assignment lands a few reads on a better-known relative.")
    if a["false_positive_mass_pct"] > 5:
        N.append(f"{pct(a['false_positive_mass_pct'], 1)} of the profile sits on organisms that are not in the tube.")
    # ---- kraken-specific
    sl, un = meta.get("species_level_frac"), meta.get("unclassified_frac")
    if tool == "kraken2" and sl is not None and sl < 0.6:
        N.append(f"Only {sl * 100:.0f}% of reads were assigned at species level; the rest stopped at genus or above. Raw Kraken2 species percentages therefore do not add up to the sample, and read-share is not abundance.")
        C.append("Run Bracken on the same Kraken2 report to redistribute genus-level reads to species before quoting abundances.")
    if un is not None and un > 0.3:
        N.append(f"{un * 100:.0f}% of reads were unclassified, so these percentages describe only the classified remainder.")
        C.append("Compare against a larger database or a lower confidence setting; a high unclassified share usually means the index is too small or the confidence threshold too strict.")
    # ---- basis
    if "bray_curtis" in a:
        if a["basis"] == "dna":
            N.append("The mock's headline shares (for example 12% per bacterium) are genomic DNA, not cell counts. A tool that reports cells or taxa, not reads, should be compared against the cell column instead.")
        else:
            N.append("This tool reports a taxon-level (cell-like) abundance, so it was scored against the cell-number column, not the DNA column.")
    # ---- cutoff advice
    ok = [s for s in a["sweep"] if s["fp"] <= 1]
    if a["sweep"] and c["fp"] > 1:
        if ok:
            best = ok[0]
            C.append(f"False species fall as the cutoff rises. The lowest cutoff here with at most one false species is {pct(best['cutoff'])} (recall {best['recall']:.2f}). This was tuned on one mock, so confirm it with a negative control from your own sample type before relying on it.")
        else:
            C.append("No cutoff up to 2% removes all but one false species. Consider a stricter confidence setting or a different tool for presence calls.")
    if a["false_positives"]:
        top_fp = a["false_positives"][0][0]
        if not (a["spizizenii_note"] and top_fp == "bacillus spizizenii"):
            C.append(f"Run a single-organism control for the largest false call ({title_case(top_fp)}) to see where the reads really come from.")
    if tool == "kraken2" and mock.kind == "quantitative":
        pass
    C.append("Repeat on a second mock or a negative control from your own lab before choosing this tool for real samples.")
    # ---- breadth of what the mock can tell you
    if mock.kind == "quantitative":
        low = min(v for v in (mock.cell if a["basis"] == "cell" else mock.dna).values())
        N.append(f"Nothing below {low}% can be judged from this tube; sensitivity to rarer organisms is untested.")
    N.append("A species label from short reads does not tell you the strain, serovar or whether an organism is harmful. Treat it as a lead to confirm, not a final identification.")
    return S, N, C, M


def render_markdown(a, mock, label=""):
    S, N, C, M = build_findings(a, mock)
    c = a["at_cutoff"]
    h = [f"# mockcheck report{': ' + label if label else ''}", "",
         f"- Tool output type: **{a['tool']}**",
         f"- Mock: **{mock.title}** ({mock.kind}), truth basis for the score: **{'cell number' if a['basis'] == 'cell' else 'genomic DNA'}**",
         f"- Species-call cutoff: **{pct(a['cutoff'])}** of the species-level profile", "", "## Scorecard", ""]
    rows = [("Organisms in tube found", f"{c['tp']} of {len(mock.species)}"),
            ("False species at cutoff", str(c["fp"])), ("Precision / recall / F1", f"{c['precision']:.2f} / {c['recall']:.2f} / {c['f1']:.2f}"),
            ("Share of profile on organisms not in the tube", pct(a["false_positive_mass_pct"], 1))]
    if "bray_curtis" in a:
        rows.insert(3, ("Bray-Curtis distance to the known answer (0 best)", f"{a['bray_curtis']:.3f}"))
    m = a["meta"] or {}
    if m.get("species_level_frac") is not None:
        rows.append(("Reads placed at species level", f"{m['species_level_frac'] * 100:.0f}%"))
    if m.get("unclassified_frac") is not None:
        rows.append(("Reads unclassified", f"{m['unclassified_frac'] * 100:.0f}%"))
    h += ["| Measure | Value |", "|---|---|"] + [f"| {k} | {v} |" for k, v in rows]
    if a["spizizenii_note"]:
        h += ["", f"_Relaxed reading (B. spizizenii counted as B. subtilis): {c['tp'] + 1} of {len(mock.species)} found, {c['fp'] - 1} false species. The strict numbers above stay the headline._"]
    h += ["", "## What this result supports", ""] + ([f"- {s}" for s in S] or ["- Nothing strong: see the next section."])
    h += ["", "## What it does not support", ""] + [f"- {s}" for s in N]
    h += ["", "## What to check next", ""] + [f"- {s}" for s in C]
    if M:
        h += ["", "## The microbiology behind the surprises", ""] + [f"- {s}" for s in M]
    h += ["", "## Cutoff sweep", "", "| Cutoff | Found | False species | Precision | Recall |", "|---|---|---|---|---|"]
    h += [f"| {pct(s['cutoff'])} | {s['tp']} | {s['fp']} | {s['precision']:.2f} | {s['recall']:.2f} |" for s in a["sweep"]]
    if "per_organism" in a:
        h += ["", "## Per organism", "", "| Organism | Expected | Observed | Fold |", "|---|---|---|---|"]
        h += [f"| {title_case(r['organism'])} | {pct(r['expected'])} | {pct(r['observed'])} | {r['fold']:.2f}x |" for r in a["per_organism"]]
    h += ["", f"_Truth source: {mock.source}. Generated by mockcheck; rule-based, no model involved._", ""]
    return "\n".join(h)


def render_isolate(i, label=""):
    ok = i["truth_pct"] >= 50
    S = [f"The true organism holds {pct(i['truth_pct'], 1)} of the profile."] if i["truth_pct"] >= 50 else []
    N = []
    if i["other_ge_cutoff"]:
        N.append(f"{i['other_ge_cutoff']} other species passed the {pct(i['cutoff'])} cutoff although only one organism was sequenced. Largest: {title_case(i['top_other'])} {pct(i['top_other_pct'])}.")
    if i["same_genus_other_pct"] > 1:
        N.append(f"{pct(i['same_genus_other_pct'], 1)} went to other species of the same genus: the method cannot separate close relatives cleanly here.")
    if i["truth_pct"] < 50:
        N.append(f"The true organism is only {pct(i['truth_pct'], 1)} of the profile; check naming differences and whether the organism is in the database.")
    C = ["Use this as the false-call floor for the tool: anything it reports above this level in a real sample is more likely signal."] if ok else \
        ["Check whether the database contains this organism, and under which name."]
    return "\n".join([f"# mockcheck isolate control{': ' + label if label else ''}", "", f"Truth: **{title_case(i['truth'])}** (single isolate). Tool output type: {i['tool']}.", "",
                      "## What this supports", ""] + [f"- {s}" for s in S or ["Nothing strong."]] +
                     ["", "## What it does not support", ""] + [f"- {s}" for s in N or ["No failure seen at this cutoff."]] +
                     ["", "## What to check next", ""] + [f"- {s}" for s in C] + [""])
