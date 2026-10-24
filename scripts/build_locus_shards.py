#!/usr/bin/env python3
"""Compact per-locus extracts from the campaign reports, for the website's drill-down view.

Reads each campaign's reports_all.tar.zst (one detection_report.yaml per genome, at any depth: v0.6.0
`wave_N/runs/<genome>/`, clean-run `wave_N/small_NNN/runs/<genome>/` or `runs/small_NNN/runs/<genome>/`,
Mucoromycotina `<SOURCE>__<id>/`) and keeps, for every detected locus, the contig, span, idiomorph,
confidence, locus class and the gene models (name, span, strand, role, status). When the campaign folder also
holds the clean-run `loci.tsv` (scripts/build_cleanrun_tables.py), each locus gets three more fields after the
gene list: tier_note, label_unsupported and high_basis (the viewer reads fields 0-6 and ignores the rest).
Writes web/intro/data/loci_<dataset>.json.gz (committed; about 1.9 MB for all 27.8k loci). Run this when a
campaign changes; build_site_data.py and the site build do not need the tarballs.

The campaign folders are the release's, from RELEASES in build_site_data.py (--release, default as there);
--asco / --basidio / --mucoro override one, as in build_site_data.py. If a genome has several run folders in
a tarball, the lexicographically last path with loci wins (the rule of build_cleanrun_tables.py).

Needs PyYAML and the zstd command. Usage: python3 scripts/build_locus_shards.py [--out web/intro/data]
    [--release cleanrun|v060] [--asco DIR] [--basidio DIR] [--mucoro DIR]
"""
import argparse
import csv
import gzip
import json
import os
import subprocess
import sys
import tarfile
from collections import defaultdict

import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_site_data import DEFAULT_RELEASE, RELEASES, REPO, resolve_release  # noqa: E402

MUCORO = {"BFD": "bfd_mucoro", "LCG": "zygolife_lcg", "JENA": "zygolife_jena"}
FIELDS = ["contig", "start", "end", "idiomorph", "confidence", "locus_class", "genes[name,start,end,strand,role,status]"]
NOTE_FIELDS = ["tier_note", "label_unsupported", "high_basis"]


def reports(tar_path):
    """Yield (member path, genome directory name, parsed detection_report.yaml) from a .tar.zst without unpacking it."""
    proc = subprocess.Popen(["zstd", "-dc", tar_path], stdout=subprocess.PIPE)
    with tarfile.open(fileobj=proc.stdout, mode="r|") as tf:
        for m in tf:
            if m.isfile() and m.name.endswith("/detection_report.yaml"):
                yield os.path.normpath(m.name), os.path.basename(os.path.dirname(m.name)), yaml.load(tf.extractfile(m), Loader=yaml.CSafeLoader) or {}
    if proc.wait():
        raise RuntimeError(f"zstd failed on {tar_path}")


def locus_row(l):
    genes = [[g["gene"], g["start"], g["end"], g["strand"], g["role"], g["status"]] for g in (l.get("gene_evidence") or [])]
    return [l["contig"], l["start"], l["end"], l.get("idiomorph") or "undetermined", l.get("confidence") or "", l.get("locus_class") or "", genes]


def load_notes(folder):
    """{genome (SOURCE__id for Mucoromycotina): [(tier_note, label_unsupported, high_basis) per locus]} from the
    campaign's clean-run loci.tsv, or None when it has no such columns."""
    path = os.path.join(folder, "loci.tsv")
    if not os.path.exists(path):
        return None
    out = defaultdict(list)
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        if "label_unsupported" not in (rd.fieldnames or []):
            return None
        for r in rd:
            g = f"{r['source']}__{r['genome']}" if r.get("source") else r["genome"]
            out[g].append([r.get(c, "") or "" for c in NOTE_FIELDS])
    return out


def campaign_shards(key, folder):
    """{dataset: {genome: [locus rows]}} for one campaign folder; the second value says whether notes were added."""
    notes = load_notes(folder)
    best = {}
    for path, name, rep in reports(os.path.join(folder, "reports_all.tar.zst")):
        loci = [locus_row(l) for l in rep.get("detected") or []]
        if not loci:
            continue
        if notes is not None and len(notes.get(name, [])) == len(loci):   # loci.tsv rows are in report order
            loci = [l + n for l, n in zip(loci, notes[name])]
        if name not in best or path > best[name][0]:
            best[name] = (path, loci)
    shards = {}
    for name, (_, loci) in best.items():
        if key == "mucoromycotina":
            src, _, genome = name.partition("__")
            shards.setdefault(MUCORO[src], {})[genome] = loci
        else:
            shards.setdefault(key, {})[name] = loci
    return shards, notes is not None


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(REPO, "web/intro/data"))
    ap.add_argument("--release", choices=sorted(RELEASES), default=DEFAULT_RELEASE,
                    help="which campaign folders to read (default %(default)s; RELEASES in build_site_data.py)")
    ap.add_argument("--asco")
    ap.add_argument("--basidio")
    ap.add_argument("--mucoro")
    a = ap.parse_args(argv)
    rel = resolve_release(a.release, a.asco, a.basidio, a.mucoro)
    os.makedirs(a.out, exist_ok=True)
    shards, with_notes = {}, {}
    for key, folder in (("bfd_asco", rel["asco"]), ("bfd_basidio", rel["basidio"]), ("mucoromycotina", rel["mucoro"])):
        part, noted = campaign_shards(key, folder)
        for ds, per in part.items():
            shards[ds] = per
            with_notes[ds] = noted
    for ds, per in sorted(shards.items()):
        fields = FIELDS + (NOTE_FIELDS if with_notes[ds] else [])
        data = json.dumps({"fields": fields, "loci": dict(sorted(per.items()))}, separators=(",", ":")).encode()
        path = os.path.join(a.out, f"loci_{ds}.json.gz")
        with open(path, "wb") as fh:  # mtime=0 keeps the file byte-identical between runs
            fh.write(gzip.compress(data, 9, mtime=0))
        print(f"{ds}: {len(per)} genomes, {sum(len(v) for v in per.values())} loci, {os.path.getsize(path) / 1e6:.2f} MB")


if __name__ == "__main__":
    main()
