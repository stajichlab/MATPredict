#!/usr/bin/env python
"""Per-family performance table from a directory of `detection_report.yaml` files, and optional held-out scores.

    family_performance_table.py REPORTS_DIR [--out table.tsv] [--md table.md] [--family Basidiomycota:redHD ...]
        [--holdout-rows ROW.json|DIR ...] [--holdout-manifest manifest.tsv]

REPORTS_DIR is searched recursively for `detection_report.yaml`, so every layout the campaigns used works: one
directory per genome (`<taxid>_<accession>/`, as `batch_runner.run_batch` writes it), batches of those
(`batch_003/<taxid>_<accession>/`), or `<genome>/<mode>/` pairs. A genome is the directory that holds the report;
if the same genome name appears twice the later path wins and the count of such repeats is printed. A genome
directory with no report (a failed run) is counted once, under `no_report`, and belongs to no family.

For each family (`Phylum:locus`), every number is a count of GENOMES unless it says `loci`:

    in_scope        reports whose `families_attempted` lists the family. Scope-only families (order.yml
                    `fallback_searchable: false`) are attempted only for genomes inside their taxonomic scope, so
                    this is the genomes the family was searched in, not the whole run
    called          in-scope genomes with at least one locus of the family
    loci            locus calls of the family (a genome can have several)
    conf_*          loci by confidence, and best_* the best confidence of each called genome
    idiomorph       loci by idiomorph/allele label (`undetermined` when the call carries no allele)
    undetermined    called genomes whose every locus of the family is `undetermined`
    partial         genomes with a `partial_locus` call; gene_only: genomes with an `idiomorph_gene_only` call
    multi_locus     genomes with two or more loci of the family
    both_alleles    genomes whose loci carry two or more different labels other than `undetermined`
    arrangement     families with a known two-gene pair (redHD: HD1/HD2; rustHD: bE/bW): loci by orientation of
                    the pair on one contig (divergent, convergent, tandem), the intergenic bp (min/median/max),
                    and `pair_incomplete` for a locus that does not carry both genes on one contig

Held-out scores (`--holdout-rows`): the JSON rows `run_asco_holdout.py run` writes. With `--holdout-manifest`
the references-left counts are recomputed from that manifest, as `run_asco_holdout.py summarize --manifest` does.
For each family and radius: `correct` (enum families: the record's idiomorph is among the called labels; pattern
families: a locus of the family is found, because allele labels are not callable by homology), `wrong`,
`undetermined`, `no_locus`, `error`, and `no_reference` (nothing of the family, or of the record's own idiomorph,
is left to search with: a database gap, not a miss). `tested` = rows that had a reference; `correct/tested` is the
recovery. The verdict rules are those of `run_asco_holdout.summarize` (a test pins them to it).
"""
from __future__ import annotations

import argparse
import csv
import json
import statistics
import sys
from collections import Counter
from pathlib import Path

import yaml

try:
    _Loader = yaml.CSafeLoader
except AttributeError:  # pragma: no cover
    _Loader = yaml.SafeLoader

REPORT = "detection_report.yaml"
CONF_ORDER = ["high", "medium", "low"]
PAIRS = {"Basidiomycota:redHD": ("HD1", "HD2"), "Basidiomycota:rustHD": ("bE", "bW")}
RADII = ["record", "species", "genus", "family", "order"]
VERDICTS = ["correct", "wrong", "undetermined", "no_locus", "error", "no_reference"]


# --------------------------------------------------------------------------- reports
def find_reports(root: Path):
    """({genome: report path}, repeats, [genome dirs without a report])."""
    found, repeats = {}, 0
    for p in sorted(root.rglob(REPORT)):
        g = p.parent.name
        repeats += g in found
        found[g] = p
    reported = {p.parent for p in found.values()}
    holders = {d.parent for d in reported}  # the directories whose children are genome directories
    empty = [d.name for h in sorted(holders) for d in sorted(h.iterdir())
             if d.is_dir() and d not in reported and not any(x.is_dir() for x in d.iterdir())]
    return found, repeats, sorted(set(empty))


def orientation(a, b):
    """Orientation of two gene evidence entries on one contig: divergent, convergent or tandem; and the gap in bp."""
    left, right = sorted((a, b), key=lambda g: g["start"])
    gap = max(0, right["start"] - left["end"] - 1)
    if left["strand"] == "-" and right["strand"] == "+":
        kind = "divergent"
    elif left["strand"] == "+" and right["strand"] == "-":
        kind = "convergent"
    else:
        kind = "tandem"
    return kind, gap


def pair_arrangement(locus, pair):
    best = {}
    for g in locus.get("gene_evidence") or []:
        if g.get("gene") in pair and not g.get("lost_to") and g.get("contig") == locus["contig"]:
            if g["gene"] not in best or (g.get("bitscore") or 0) > (best[g["gene"]].get("bitscore") or 0):
                best[g["gene"]] = g
    if len(best) < 2:
        return "pair_incomplete", None
    return orientation(best[pair[0]], best[pair[1]])


class FamilyStats:
    def __init__(self, family):
        self.family = family
        self.in_scope = 0
        self.called = 0
        self.loci = 0
        self.conf = Counter()
        self.best = Counter()
        self.idiomorph = Counter()
        self.locus_class = Counter()
        self.undetermined = 0
        self.partial = 0
        self.gene_only = 0
        self.multi_locus = 0
        self.both_alleles = 0
        self.arr = Counter()
        self.gaps = []
        self.genomes_called = []
        self.genomes_multi = []

    def add_genome(self, genome, in_scope, loci):
        self.in_scope += in_scope
        if not loci:
            return
        self.called += 1
        self.loci += len(loci)
        self.genomes_called.append(genome)
        for c in loci:
            self.conf[c.get("confidence") or "none"] += 1
            self.idiomorph[str(c.get("idiomorph"))] += 1
            self.locus_class[c.get("locus_class") or "none"] += 1
        confs = [c.get("confidence") for c in loci]
        self.best[next((k for k in CONF_ORDER if k in confs), "none")] += 1
        labels = {str(c.get("idiomorph")) for c in loci}
        self.undetermined += labels <= {"undetermined", "None"}
        self.partial += any(c.get("locus_class") == "partial_locus" for c in loci)
        self.gene_only += any(c.get("locus_class") == "idiomorph_gene_only" for c in loci)
        if len(loci) > 1:
            self.multi_locus += 1
            self.genomes_multi.append(genome)
        self.both_alleles += len(labels - {"undetermined", "None"}) >= 2
        pair = PAIRS.get(self.family)
        if pair:
            for c in loci:
                kind, gap = pair_arrangement(c, pair)
                self.arr[kind] += 1
                if gap is not None:
                    self.gaps.append(gap)

    def row(self):
        r = dict(family=self.family, in_scope=self.in_scope, called=self.called, loci=self.loci,
                 conf_high=self.conf["high"], conf_medium=self.conf["medium"], conf_low=self.conf["low"],
                 best_high=self.best["high"], best_medium=self.best["medium"], best_low=self.best["low"],
                 undetermined=self.undetermined, partial=self.partial, gene_only=self.gene_only,
                 multi_locus=self.multi_locus, both_alleles=self.both_alleles,
                 idiomorph=";".join(f"{k}={v}" for k, v in sorted(self.idiomorph.items())),
                 locus_class=";".join(f"{k}={v}" for k, v in sorted(self.locus_class.items())))
        if self.family in PAIRS:
            g = self.gaps
            r["arrangement"] = ";".join(f"{k}={v}" for k, v in sorted(self.arr.items())) + (
                f";gap_bp_min/median/max={min(g)}/{statistics.median(g):g}/{max(g)}" if g else "")
        else:
            r["arrangement"] = ""
        return r


def summarize_reports(root: Path, only: set[str] | None = None):
    found, repeats, empty = find_reports(root)
    stats: dict[str, FamilyStats] = {}
    not_searched = 0
    for genome, path in found.items():
        doc = yaml.load(path.read_text(), Loader=_Loader) or {}
        not_searched += bool(doc.get("not_searched_reason"))
        by_family: dict[str, list] = {}
        for d in doc.get("detected") or []:
            by_family.setdefault(d["family"], []).append(d)
        for fam in set(doc.get("families_attempted") or []) | set(by_family):
            if only and fam not in only:
                continue
            attempted = fam in (doc.get("families_attempted") or [])
            stats.setdefault(fam, FamilyStats(fam)).add_genome(genome, attempted or fam in by_family, by_family.get(fam, []))
    return dict(stats=stats, genomes=len(found), repeats=repeats, no_report=empty, not_searched=not_searched)


# --------------------------------------------------------------------------- held out
def holdout_verdict(r: dict) -> str:
    """The rules of run_asco_holdout.summarize, on one row."""
    exp = set(r["expected"].split("|")) - {""}
    labels = r.get("pipeline_labels") or []
    found = bool(r.get("loci"))
    if r["status"] in ("detect_error", "no_genome_in_library", "no_reference_family"):
        v = r["status"]
    elif not found:
        v = "no_reference_family" if r.get("ref_family_left") == 0 else "no_locus"
    elif not r.get("enum_family"):
        v = "locus_found"
    elif exp & set(labels):
        v = "correct"
    elif labels and set(labels) <= {"undetermined"}:
        v = "undetermined"
    else:
        v = "wrong"
    if v in ("no_locus", "wrong") and r.get("enum_family") and r.get("ref_idiomorph_left") == 0:
        v = "no_reference_idiomorph"
    return v


def holdout_bucket(v: str) -> str:
    return {"locus_found": "correct", "no_reference_family": "no_reference", "no_reference_idiomorph": "no_reference",
            "detect_error": "error", "no_genome_in_library": "error"}.get(v, v)


def load_holdout(paths, manifest: Path | None):
    files = []
    for p in map(Path, paths):
        files += sorted(p.rglob("*.json")) if p.is_dir() else [p]
    rows = [json.loads(f.read_text()) for f in files]
    rows = [r for r in rows if "record_id" in r and "radius" in r]
    if manifest:
        man = {m["record_id"]: m for m in csv.DictReader(open(manifest), delimiter="\t")}
        for r in rows:
            m = man.get(r["record_id"])
            if m and f"{r['radius']}_same_family_left" in m:
                r["ref_family_left"] = int(m[f"{r['radius']}_same_family_left"])
                r["ref_idiomorph_left"] = int(m[f"{r['radius']}_same_idiomorph_left"])
    return rows


def summarize_holdout(rows):
    out: dict[tuple, Counter] = {}
    for r in rows:
        r["verdict"] = holdout_verdict(r)
        c = out.setdefault((r["family"], r["radius"]), Counter())
        c[holdout_bucket(r["verdict"])] += 1
    result = []
    for (fam, rad), c in sorted(out.items(), key=lambda kv: (kv[0][0], RADII.index(kv[0][1]) if kv[0][1] in RADII else 9)):
        tested = sum(c[k] for k in ("correct", "wrong", "undetermined", "no_locus"))
        result.append(dict(family=fam, radius=rad, rows=sum(c.values()), tested=tested, correct=c["correct"],
                           wrong=c["wrong"], undetermined=c["undetermined"], no_locus=c["no_locus"],
                           error=c["error"], no_reference=c["no_reference"]))
    return result


# --------------------------------------------------------------------------- output
REPORT_COLS = ["family", "in_scope", "called", "loci", "conf_high", "conf_medium", "conf_low", "best_high", "best_medium",
               "best_low", "undetermined", "partial", "gene_only", "multi_locus", "both_alleles", "idiomorph",
               "locus_class", "arrangement"]
HOLD_COLS = ["family", "radius", "rows", "tested", "correct", "wrong", "undetermined", "no_locus", "error", "no_reference"]


def write_tsv(path, cols, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def markdown(cols, rows):
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    lines += ["| " + " | ".join(str(r.get(c, "")) for c in cols) + " |" for r in rows]
    return "\n".join(lines)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reports", type=Path)
    ap.add_argument("--out", type=Path, help="TSV of the per-family table (the holdout table goes to <out>.holdout.tsv)")
    ap.add_argument("--md", type=Path, help="markdown of both tables")
    ap.add_argument("--family", action="append", help="only this Phylum:locus (repeatable)")
    ap.add_argument("--holdout-rows", nargs="*", default=[], help="row JSON files or directories from run_asco_holdout.py run")
    ap.add_argument("--holdout-manifest", type=Path)
    a = ap.parse_args(argv)
    res = summarize_reports(a.reports, set(a.family) if a.family else None)
    rows = [s.row() for _, s in sorted(res["stats"].items())]
    hold = summarize_holdout(load_holdout(a.holdout_rows, a.holdout_manifest)) if a.holdout_rows else []
    if a.family:
        hold = [h for h in hold if h["family"] in set(a.family)]
    head = (f"{res['genomes']} reports, {len(res['no_report'])} genome directories without a report, "
            f"{res['repeats']} repeated genome names, {res['not_searched']} reports with no family searched")
    print(head)
    text = markdown(REPORT_COLS, rows)
    if hold:
        text += "\n\nheld out\n\n" + markdown(HOLD_COLS, hold)
    print(text)
    if a.out:
        write_tsv(a.out, REPORT_COLS, rows)
        if hold:
            write_tsv(a.out.with_suffix(".holdout.tsv"), HOLD_COLS, hold)
    if a.md:
        a.md.write_text(head + "\n\n" + text + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
