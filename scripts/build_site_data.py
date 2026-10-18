#!/usr/bin/env python3
"""Aggregate the committed campaign tables into the JSON the project web pages read.

Inputs are tables already under results/ (nothing is recomputed):
  BFD Ascomycota      results/2026-10-03_ascomycota_v060/{genomes,loci}.tsv
  BFD Basidiomycota   results/2026-10-03_basidiomycota_v060/{genomes,loci}.tsv
  BFD Mucoromycotina  results/2026-10-03_mucoromycotina_mat/calls.tsv   (source BFD)
  ZyGoLife (LCG, Jena held-out sets, Mucoromycotina)  same calls.tsv, source LCG / JENA

Also writes one per-genome table per dataset, genomes_<dataset>.tsv (git-ignored; built with the site).
Also writes cocci_reads.json and cocci_reads.tsv from datasets/reads/coccidioides/results when present.
Writes web/intro/data/taxa.json (one row per dataset x lineage x species, with genome counts and allele
counts) and summary.json (per-dataset totals, locus classes, confidence). Python standard library only.

Usage: python3 scripts/build_site_data.py [--out web/intro/data]
"""
import argparse
import csv
import datetime
import json
import os
import re
from collections import Counter, defaultdict

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(REPO, "results")

NOT_A_GENUS = {"Ascomycota", "Basidiomycota", "Mucoromycota", "Fungal", "Fungi", "uncultured", "Candidatus", "unidentified"}
GENUS_RE = re.compile(r"^[A-Z][a-z-]+$")


def tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def clean(v):
    v = (v or "").strip()
    return "" if v in ("NA", "None") else v


def genus_of(species):
    tok = (species or "").replace("[", "").replace("]", "").split()
    if not tok or tok[0] in NOT_A_GENUS or not GENUS_RE.match(tok[0]):
        return "(unassigned)"
    return tok[0]


def species_of(species):
    tok = (species or "").replace("[", "").replace("]", "").split()
    if len(tok) < 2 or genus_of(species) == "(unassigned)":
        return "(unassigned)"
    return f"{tok[0]} {tok[1]}"


class Agg:
    """Counts genomes per (dataset, lineage, species). A genome is one of: uncalled, called with no
    determined allele (undet), called with one allele, called with two or more distinct alleles (two)."""

    def __init__(self):
        self.rows = defaultdict(lambda: {"n": 0, "called": 0, "undet": 0, "two": 0, "alleles": Counter(), "types": Counter()})

    def add(self, key, called, alleles, types):
        r = self.rows[key]
        r["n"] += 1
        if not called:
            return
        r["called"] += 1
        for t in sorted(set(types)):
            r["types"][t] += 1
        if not alleles:
            r["undet"] += 1
        elif len(alleles) > 1:
            r["two"] += 1
            r["alleles"]["two:" + "+".join(sorted(alleles))] += 1
        else:
            r["alleles"][next(iter(alleles))] += 1


def loci_by_genome(rows, genome_col="genome"):
    out = defaultdict(list)
    for r in rows:
        out[r[genome_col]].append(r)
    return out


GENOME_COLS = ["dataset", "genome", "species", "species_recorded", "phylum", "class", "order", "family", "genus", "status", "call", "alleles",
               "locus_types", "n_loci", "best_confidence", "loci"]
_CONF = {"high": 3, "medium": 2, "low": 1}


def genome_row(ds_id, genome, species, phylum, cls, order, family, status, called, alleles, types, loci):
    """loci: list of (contig, start, end, idiomorph, confidence, locus_class). One row per genome."""
    if not called:
        call = status if status not in ("called", "") else "uncalled"
    elif not alleles:
        call = "undetermined"
    else:
        call = "+".join(sorted(alleles))
    best = max((l[4] for l in loci), key=lambda c: _CONF.get(c, 0), default="")
    return [ds_id, genome, species_of(species), species, phylum, cls, order, family, genus_of(species), status, call, ",".join(sorted(alleles)),
            ",".join(sorted(set(types))), len(loci), best, ";".join(f"{c}:{a}-{b}|{i or 'undetermined'}|{cf}|{lc}" for c, a, b, i, cf, lc in loci)]


def strip_family(f):
    return f.split(":", 1)[-1]


COCCI_READS = os.path.join(REPO, "datasets/reads/coccidioides/results/cocci2_reads_vs_detect.tsv")


def cocci_reads(out):
    """Coccidioides read typing (datasets/reads/coccidioides): counts per species and call, the omitted and both strains, agreement with
    `detect`; plus a compact per-strain table. Nothing is recomputed."""
    if not os.path.exists(COCCI_READS):
        return
    rows = tsv(COCCI_READS)
    species = sorted({r["species"] for r in rows})
    calls = ["MAT1-1", "MAT1-2", "both", "low_depth"]
    counts = {s: {c: sum(1 for r in rows if r["species"] == s and r["reads_call"] == c) for c in calls} for s in species}
    typed_det = [r for r in rows if r["triage"] == "typed" and r["detect_call"]]
    agree = sum(1 for r in typed_det if r["reads_call"] == r["detect_call"])
    summ = {
        "panel": "cocci2", "strains": len(rows), "species": species, "counts": counts,
        "disagree_no_locus": sum(1 for r in typed_det if r["reads_call"] != r["detect_call"] and r["detect_call"] == "not_detected"),
        "disagree_other": sum(1 for r in typed_det if r["reads_call"] != r["detect_call"] and r["detect_call"] != "not_detected"),
        "with_assembly": sum(1 for r in rows if r["detect_call"]), "typed_with_detect": len(typed_det), "agree_with_detect": agree,
        "relabelled": [r["strain"] for r in rows if r["relabelled"] == "yes"],
        "omitted": [{"strain": r["strain"], "species": r["species"], "center": r["center"], "shared_depth": float(r["shared_depth"])} for r in rows if r["triage"] == "review"],
        "both": [{"strain": r["strain"], "species": r["species"], "MAT1-1_depth": float(r["MAT1-1_depth"]), "MAT1-2_depth": float(r["MAT1-2_depth"])} for r in rows if r["triage"] == "both"],
        "centers": dict(sorted(Counter(r["center"] or "unknown" for r in rows).items())),
    }
    with open(os.path.join(out, "cocci_reads.json"), "w") as fh:
        json.dump(summ, fh, indent=1)
    cols = ["strain", "species", "relabelled", "center", "reads_call", "triage", "MAT1-1_depth", "MAT1-2_depth", "shared_depth", "detect_call", "detect_confidence"]
    with open(os.path.join(out, "cocci_reads.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        w.writerows([r[c] for c in cols] for r in rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(REPO, "web/intro/data"))
    out = ap.parse_args().out
    os.makedirs(out, exist_ok=True)

    agg = Agg()
    grows = defaultdict(list)
    summary = {"datasets": {}, "generated": datetime.date.today().isoformat()}

    def dataset(ds_id, label, phylum, source, code, n_genomes):
        summary["datasets"][ds_id] = {
            "label": label, "phylum": phylum, "source": source, "code": code, "genomes": n_genomes,
            "loci": 0, "locus_class": Counter(), "confidence": Counter(), "status": Counter(), "routing": Counter(),
        }
        return summary["datasets"][ds_id]

    # ---- Ascomycota (BFD) ----
    d = os.path.join(R, "2026-10-03_ascomycota_v060")
    genomes, loci = tsv(os.path.join(d, "genomes.tsv")), loci_by_genome(tsv(os.path.join(d, "loci.tsv")))
    s = dataset("bfd_asco", "BFD Ascomycota", "Ascomycota", "BFD", "v0.6.0 (2026-10-03)", len(genomes))
    for g in genomes:
        ls = [l for l in loci.get(g["genome"], [])]
        called = g["status"] == "called"
        alleles = {l["idiomorph"] for l in ls if l["idiomorph"] not in ("", "undetermined")}
        key = ("BFD", "Ascomycota", clean(g["class_"]), clean(g["order"]), clean(g["family"]), genus_of(g["species"]), species_of(g["species"]))
        agg.add(key, called, alleles, ["MAT locus"] if called else [])
        grows["bfd_asco"].append(genome_row("bfd_asco", g["genome"], g["species"], "Ascomycota", clean(g["class_"]), clean(g["order"]), clean(g["family"]),
                                            g["status"], called, alleles, ["MAT locus"] if called else [],
                                            [(l["contig"], l["start"], l["end"], l["idiomorph"], l["confidence"], l["locus_class"]) for l in ls]))
        s["status"][g["status"]] += 1
        s["routing"][g["routing"] or "unrecorded"] += 1
        for l in ls:
            s["loci"] += 1
            s["locus_class"][l["locus_class"]] += 1
            s["confidence"][l["confidence"]] += 1

    # ---- Basidiomycota (BFD) ----
    d = os.path.join(R, "2026-10-03_basidiomycota_v060")
    genomes, loci = tsv(os.path.join(d, "genomes.tsv")), loci_by_genome(tsv(os.path.join(d, "loci.tsv")))
    s = dataset("bfd_basidio", "BFD Basidiomycota", "Basidiomycota", "BFD", "v0.6.0 (2026-10-03)", len(genomes))
    for g in genomes:
        ls = loci.get(g["genome"], [])
        called = g["status"] == "called"
        alleles = {l["idiomorph"] for l in ls if l["idiomorph"] not in ("", "undetermined")}
        types = [strip_family(l["family_called"]) for l in ls]
        key = ("BFD", "Basidiomycota", clean(g["class_"]), clean(g["order"]), clean(g["family"]), genus_of(g["species"]), species_of(g["species"]))
        agg.add(key, called, alleles, types)
        grows["bfd_basidio"].append(genome_row("bfd_basidio", g["genome"], g["species"], "Basidiomycota", clean(g["class_"]), clean(g["order"]), clean(g["family"]),
                                               g["status"], called, alleles, types,
                                               [(l["contig"], l["start"], l["end"], l["idiomorph"], l["confidence"], l["locus_class"]) for l in ls]))
        s["status"][g["status"]] += 1
        s["routing"][g["routing"] or "unrecorded"] += 1
        for l in ls:
            s["loci"] += 1
            s["locus_class"][l["locus_class"]] += 1
            s["confidence"][l["confidence"]] += 1

    # ---- Mucoromycotina: BFD rows and ZyGoLife (LCG + Jena) ----
    calls = tsv(os.path.join(R, "2026-10-03_mucoromycotina_mat/calls.tsv"))
    per = defaultdict(list)
    for r in calls:
        if r["source"] == "source":
            continue
        per[(r["source"], r["genome"])].append(r)
    meta = {"BFD": ("bfd_mucoro", "BFD Mucoromycotina", "BFD", "Mucoromycota"),
            "LCG": ("zygolife_lcg", "ZyGoLife LCG (held out)", "ZyGoLife", "Mucoromycota"),
            "JENA": ("zygolife_jena", "ZyGoLife Jena (held out)", "ZyGoLife", "Mucoromycota")}
    counts = Counter(src for src, _ in per)
    for src, (ds_id, label, source, phylum) in meta.items():
        dataset(ds_id, label, phylum, source, "v0.6.0 (2026-10-03)", counts[src])
    for (src, genome), ls in per.items():
        ds_id, _, source, phylum = meta[src]
        s = summary["datasets"][ds_id]
        called = any(l["status"] == "called" for l in ls)
        loc = [l for l in ls if l["status"] == "called"]
        alleles = {l["idiomorph"] for l in loc if l["idiomorph"] not in ("", "undetermined")}
        order = clean(ls[0]["bfd_order"])
        name = ls[0]["name"]
        key = (source, phylum, "Mucoromycotina", order, "", genus_of(name), species_of(name), src)
        agg.add(key, called, alleles, ["sexP/sexM MAT locus"] if called else [])
        grows[ds_id].append(genome_row(ds_id, genome, name, phylum, "Mucoromycotina", order, "", "called" if called else "uncalled", called, alleles,
                                       ["sexP/sexM MAT locus"] if called else [],
                                       [(l["contig"], l["locus_start"], l["locus_end"], l["idiomorph"], l["confidence"], l["locus_class"]) for l in loc]))
        s["status"]["called" if called else "uncalled"] += 1
        for l in loc:
            s["loci"] += 1
            s["locus_class"][l["locus_class"]] += 1
            s["confidence"][l["confidence"]] += 1

    # compact rows: dictionary-encoded strings keep the file small
    strings, idx = [], {}

    def enc(v):
        if v not in idx:
            idx[v] = len(strings)
            strings.append(v)
        return idx[v]

    rows = []
    for key in sorted(agg.rows, key=lambda k: tuple(str(x) for x in k)):
        r = agg.rows[key]
        cohort = key[7] if len(key) > 7 else ""
        rows.append([enc(x) for x in key[:7]] + [enc(cohort), r["n"], r["called"], r["undet"], r["two"],
                     {enc(a): n for a, n in r["alleles"].items() if n}, {enc(a): n for a, n in r["types"].items()}])
    taxa = {"fields": ["dataset", "phylum", "class", "order", "family", "genus", "species", "cohort", "n", "called", "undet", "two", "alleles", "types"],
            "strings": strings, "rows": rows}
    with open(os.path.join(out, "taxa.json"), "w") as fh:
        json.dump(taxa, fh, separators=(",", ":"))
    for ds_id, grs in grows.items():
        with open(os.path.join(out, f"genomes_{ds_id}.tsv"), "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(GENOME_COLS)
            w.writerows(sorted(grs, key=lambda r: r[1]))
        summary["datasets"][ds_id]["genome_table"] = f"data/genomes_{ds_id}.tsv"
    for ds in summary["datasets"].values():
        for k in ("locus_class", "confidence", "status", "routing"):
            ds[k] = dict(ds[k])
    with open(os.path.join(out, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    cocci_reads(out)
    print(f"{len(rows)} rows, {len(strings)} strings; {sum(r[8] for r in rows)} genomes")


if __name__ == "__main__":
    main()
