#!/usr/bin/env python3
"""Aggregate the committed campaign tables into the JSON the project web pages read.

Inputs are tables already under results/ (nothing is recomputed). Which folders are read is set in ONE place,
`RELEASES` below, and chosen with --release (default `cleanrun`); --asco / --basidio / --mucoro override one
folder. Per release:
  BFD Ascomycota      <asco>/{genomes,loci}.tsv
  BFD Basidiomycota   <basidio>/{genomes,loci}.tsv
  BFD Mucoromycotina  <mucoro>/{genomes,loci}.tsv with a `source` column (scripts/build_cleanrun_tables.py
                      --layout mucoro), or the v0.6.0 <mucoro>/calls.tsv                     (source BFD)
  ZyGoLife (LCG, Jena held-out sets, Mucoromycotina)  the same Mucoromycotina tables, source LCG / JENA
A missing folder is an error: the site is never built from a mix of releases by accident.

The clean-run tables (scripts/build_cleanrun_tables.py) carry the curator's post-run reporting columns; when a
table has them they are passed on, and they are simply absent for v0.6.0:
  genomes_<dataset>.tsv  three columns appended after the v0.6.0 ones (existing columns unchanged):
                         `label_unsupported` (number of the genome's loci marked, #138), `tier_notes` and
                         `high_basis` (distinct values over its loci), and `locus_notes`: one entry per locus,
                         `;`-separated in the order of the `loci` column, holding that locus' non-empty values of
                         tier_note, label_unsupported (its own `;` written as `+`) and high_basis joined by `,`,
                         or `-` (the three vocabularies do not overlap).
  summary.json           per dataset: `label_unsupported` (loci), `tier_note` and `high_basis` (locus counts);
                         top level: `release` (name, code label, folders).
  taxa.json              one field appended to each row: `unsupported`, the genomes with at least one marked
                         locus (readers that index the first 14 fields are unaffected).

Also writes one per-genome table per dataset, genomes_<dataset>.tsv (git-ignored; built with the site).
Also writes cocci_reads.json and cocci_reads.tsv from datasets/reads/coccidioides/results when present.
Writes web/intro/data/taxa.json (one row per dataset x lineage x species, with genome counts and allele
counts) and summary.json (per-dataset totals, locus classes, confidence). Python standard library only.

Usage: python3 scripts/build_site_data.py [--out web/intro/data] [--release cleanrun|v060]
                                          [--asco DIR] [--basidio DIR] [--mucoro DIR]
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

# The result folders behind each published release (relative to results/). Edit HERE when a campaign folder
# changes; build_locus_shards.py reads the same table.
RELEASES = {
    "cleanrun": {
        "code": "cleanrun-2026-10b (56e11cd, October 2026)",
        "asco": "2026-10-23_ascomycota_cleanrun",
        "basidio": "2026-10-23_basidiomycota_cleanrun",
        "mucoro": None,   # the clean-run Mucoromycotina campaign has not been run yet: pass --mucoro or set it here
    },
    "v060": {
        "code": "v0.6.0 (2026-10-03)",
        "asco": "2026-10-03_ascomycota_v060",
        "basidio": "2026-10-03_basidiomycota_v060",
        "mucoro": "2026-10-03_mucoromycotina_mat",
    },
}
DEFAULT_RELEASE = "cleanrun"

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
    determined allele (undet), called with one allele, called with two or more distinct alleles (two).
    `unsupported` counts the genomes with at least one locus whose label is marked unsupported (#138)."""

    def __init__(self):
        self.rows = defaultdict(lambda: {"n": 0, "called": 0, "undet": 0, "two": 0, "unsupported": 0, "alleles": Counter(), "types": Counter()})

    def add(self, key, called, alleles, types, unsupported=False):
        r = self.rows[key]
        r["n"] += 1
        if not called:
            return
        r["called"] += 1
        r["unsupported"] += bool(unsupported)
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
               "locus_types", "n_loci", "best_confidence", "loci", "label_unsupported", "tier_notes", "high_basis", "locus_notes"]
_CONF = {"high": 3, "medium": 2, "low": 1}
NOTE_COLS = ("tier_note", "label_unsupported", "high_basis")


def genome_row(ds_id, genome, species, phylum, cls, order, family, status, called, alleles, types, loci, notes=None):
    """loci: list of (contig, start, end, idiomorph, confidence, locus_class); notes: one
    (tier_note, label_unsupported, high_basis) per locus, or None when the tables have no such columns."""
    if not called:
        call = status if status not in ("called", "") else "uncalled"
    elif not alleles:
        call = "undetermined"
    else:
        call = "+".join(sorted(alleles))
    best = max((l[4] for l in loci), key=lambda c: _CONF.get(c, 0), default="")
    row = [ds_id, genome, species_of(species), species, phylum, cls, order, family, genus_of(species), status, call, ",".join(sorted(alleles)),
           ",".join(sorted(set(types))), len(loci), best, ";".join(f"{c}:{a}-{b}|{i or 'undetermined'}|{cf}|{lc}" for c, a, b, i, cf, lc in loci)]
    if notes is None:
        return row + ["", "", "", ""]
    return row + [sum(bool(n[1]) for n in notes), "|".join(sorted({n[0] for n in notes if n[0]})),
                  "|".join(sorted({n[2] for n in notes if n[2]})),
                  ";".join(",".join(v.replace(";", "+") for v in n if v) or "-" for n in notes)]


def locus_notes(ls):
    """(tier_note, label_unsupported, high_basis) per locus row, or None for tables without the columns."""
    if not ls or "label_unsupported" not in ls[0]:
        return None
    return [tuple(l.get(c, "") or "" for c in NOTE_COLS) for l in ls]


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


def resolve_release(release, asco=None, basidio=None, mucoro=None):
    """{"name", "code", "asco", "basidio", "mucoro"} with absolute folders; raises when one is unset or missing."""
    rel = dict(RELEASES[release], name=release)
    for key, override in (("asco", asco), ("basidio", basidio), ("mucoro", mucoro)):
        if override:
            rel[key] = override
        if not rel[key]:
            raise SystemExit(f"release {release!r}: no {key} folder is set; pass --{key} or set it in RELEASES")
        # an override is a path (absolute, or relative to the working directory); a RELEASES entry is under results/
        path = os.path.abspath(rel[key]) if override else os.path.join(R, rel[key])
        if not os.path.isdir(path):
            raise SystemExit(f"release {release!r}: {key} folder {path} does not exist")
        rel[key] = path
    return rel


def mucoro_rows(folder):
    """Per-locus rows keyed like the v0.6.0 calls.tsv (source, genome, name, bfd_order, status, idiomorph, confidence,
    locus_class, contig, locus_start, locus_end, plus the note columns when present), from either table format."""
    if not os.path.exists(os.path.join(folder, "genomes.tsv")):
        return tsv(os.path.join(folder, "calls.tsv"))
    out = []
    loci = defaultdict(list)
    for l in tsv(os.path.join(folder, "loci.tsv")):
        loci[(l["source"], l["genome"])].append(l)
    for g in tsv(os.path.join(folder, "genomes.tsv")):
        base = dict(source=g["source"], genome=g["genome"], name=g["species"], bfd_order=g["order"])
        ls = loci.get((g["source"], g["genome"]), [])
        if g["status"] != "called" or not ls:
            out.append(dict(base, status=g["status"] if g["status"] != "called" else "uncalled"))
            continue
        for l in ls:
            out.append(dict(base, status="called", idiomorph=l["idiomorph"], confidence=l["confidence"], locus_class=l["locus_class"],
                            contig=l["contig"], locus_start=l["start"], locus_end=l["end"],
                            **{c: l.get(c, "") for c in NOTE_COLS}))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(REPO, "web/intro/data"))
    ap.add_argument("--release", choices=sorted(RELEASES), default=DEFAULT_RELEASE,
                    help="which campaign folders to read (default %(default)s; see RELEASES)")
    ap.add_argument("--asco", help="Ascomycota results folder (overrides the release's)")
    ap.add_argument("--basidio", help="Basidiomycota results folder (overrides the release's)")
    ap.add_argument("--mucoro", help="Mucoromycotina results folder (overrides the release's)")
    a = ap.parse_args(argv)
    out = a.out
    rel = resolve_release(a.release, a.asco, a.basidio, a.mucoro)
    os.makedirs(out, exist_ok=True)

    agg = Agg()
    grows = defaultdict(list)
    summary = {"datasets": {}, "generated": datetime.date.today().isoformat(),
               "release": {"name": rel["name"], "code": rel["code"],
                           "folders": {k: os.path.relpath(rel[k], REPO) for k in ("asco", "basidio", "mucoro")}}}
    code = rel["code"]

    def dataset(ds_id, label, phylum, source, code, n_genomes):
        summary["datasets"][ds_id] = {
            "label": label, "phylum": phylum, "source": source, "code": code, "genomes": n_genomes,
            "loci": 0, "locus_class": Counter(), "confidence": Counter(), "status": Counter(), "routing": Counter(),
            "label_unsupported": 0, "tier_note": Counter(), "high_basis": Counter(),
        }
        return summary["datasets"][ds_id]

    def count_locus(s, l):
        s["loci"] += 1
        s["locus_class"][l["locus_class"]] += 1
        s["confidence"][l["confidence"]] += 1
        s["label_unsupported"] += bool(l.get("label_unsupported"))
        if l.get("tier_note"):
            s["tier_note"][l["tier_note"]] += 1
        if l.get("high_basis"):
            s["high_basis"][l["high_basis"]] += 1

    # ---- BFD Ascomycota and Basidiomycota: genomes.tsv + loci.tsv ----
    for ds_id, label, phylum, folder, types_of in (
            ("bfd_asco", "BFD Ascomycota", "Ascomycota", rel["asco"], lambda ls, called: ["MAT locus"] if called else []),
            ("bfd_basidio", "BFD Basidiomycota", "Basidiomycota", rel["basidio"], lambda ls, called: [strip_family(l["family_called"]) for l in ls])):
        genomes, loci = tsv(os.path.join(folder, "genomes.tsv")), loci_by_genome(tsv(os.path.join(folder, "loci.tsv")))
        s = dataset(ds_id, label, phylum, "BFD", code, len(genomes))
        for g in genomes:
            ls = loci.get(g["genome"], [])
            called = g["status"] == "called"
            alleles = {l["idiomorph"] for l in ls if l["idiomorph"] not in ("", "undetermined")}
            types = types_of(ls, called)
            notes = locus_notes(ls)
            key = ("BFD", phylum, clean(g["class_"]), clean(g["order"]), clean(g["family"]), genus_of(g["species"]), species_of(g["species"]))
            agg.add(key, called, alleles, types, any(n[1] for n in notes or ()))
            grows[ds_id].append(genome_row(ds_id, g["genome"], g["species"], phylum, clean(g["class_"]), clean(g["order"]), clean(g["family"]),
                                           g["status"], called, alleles, types,
                                           [(l["contig"], l["start"], l["end"], l["idiomorph"], l["confidence"], l["locus_class"]) for l in ls], notes))
            s["status"][g["status"]] += 1
            s["routing"][g["routing"] or "unrecorded"] += 1
            for l in ls:
                count_locus(s, l)

    # ---- Mucoromycotina: BFD rows and ZyGoLife (LCG + Jena) ----
    calls = mucoro_rows(rel["mucoro"])
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
        dataset(ds_id, label, phylum, source, code, counts[src])
    for (src, genome), ls in per.items():
        ds_id, _, source, phylum = meta[src]
        s = summary["datasets"][ds_id]
        called = any(l["status"] == "called" for l in ls)
        loc = [l for l in ls if l["status"] == "called"]
        alleles = {l["idiomorph"] for l in loc if l["idiomorph"] not in ("", "undetermined")}
        order = clean(ls[0]["bfd_order"])
        name = ls[0]["name"]
        notes = locus_notes(loc)
        key = (source, phylum, "Mucoromycotina", order, "", genus_of(name), species_of(name), src)
        agg.add(key, called, alleles, ["sexP/sexM MAT locus"] if called else [], any(n[1] for n in notes or ()))
        grows[ds_id].append(genome_row(ds_id, genome, name, phylum, "Mucoromycotina", order, "", "called" if called else "uncalled", called, alleles,
                                       ["sexP/sexM MAT locus"] if called else [],
                                       [(l["contig"], l["locus_start"], l["locus_end"], l["idiomorph"], l["confidence"], l["locus_class"]) for l in loc],
                                       notes))
        s["status"]["called" if called else "uncalled"] += 1
        for l in loc:
            count_locus(s, l)

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
                     {enc(a): n for a, n in r["alleles"].items() if n}, {enc(a): n for a, n in r["types"].items()}, r["unsupported"]])
    taxa = {"fields": ["dataset", "phylum", "class", "order", "family", "genus", "species", "cohort", "n", "called", "undet", "two", "alleles", "types",
                       "unsupported"],
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
        for k in ("locus_class", "confidence", "status", "routing", "tier_note", "high_basis"):
            ds[k] = dict(sorted(ds[k].items())) if k in ("tier_note", "high_basis") else dict(ds[k])
    with open(os.path.join(out, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    cocci_reads(out)
    print(f"{rel['name']}: {len(rows)} rows, {len(strings)} strings; {sum(r[8] for r in rows)} genomes")


if __name__ == "__main__":
    main()
