#!/usr/bin/env python3
"""Compact per-locus extracts from the campaign reports, for the website's drill-down view.

Reads each campaign's reports_all.tar.zst (one detection_report.yaml per genome) and keeps, for every detected
locus, the contig, span, idiomorph, confidence, locus class and the gene models (name, span, strand, role, status).
Writes web/intro/data/loci_<dataset>.json.gz (committed; about 1.9 MB for all 27.8k loci). Run this when a
campaign changes; build_site_data.py and the site build do not need the tarballs.

Needs PyYAML and the zstd command. Usage: python3 scripts/build_locus_shards.py [--out web/intro/data]
"""
import argparse
import gzip
import json
import os
import subprocess
import tarfile

import yaml

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(REPO, "results")
CAMPAIGNS = [("bfd_asco", "2026-10-03_ascomycota_v060"), ("bfd_basidio", "2026-10-03_basidiomycota_v060"), ("mucoromycotina", "2026-10-03_mucoromycotina_mat")]
MUCORO = {"BFD": "bfd_mucoro", "LCG": "zygolife_lcg", "JENA": "zygolife_jena"}
FIELDS = ["contig", "start", "end", "idiomorph", "confidence", "locus_class", "genes[name,start,end,strand,role,status]"]


def reports(tar_path):
    """Yield (genome directory name, parsed detection_report.yaml) from a .tar.zst without unpacking it to disk."""
    proc = subprocess.Popen(["zstd", "-dc", tar_path], stdout=subprocess.PIPE)
    with tarfile.open(fileobj=proc.stdout, mode="r|") as tf:
        for m in tf:
            if m.isfile() and m.name.endswith("/detection_report.yaml"):
                yield os.path.basename(os.path.dirname(m.name)), yaml.load(tf.extractfile(m), Loader=yaml.CSafeLoader) or {}
    if proc.wait():
        raise RuntimeError(f"zstd failed on {tar_path}")


def locus_row(l):
    genes = [[g["gene"], g["start"], g["end"], g["strand"], g["role"], g["status"]] for g in (l.get("gene_evidence") or [])]
    return [l["contig"], l["start"], l["end"], l.get("idiomorph") or "undetermined", l.get("confidence") or "", l.get("locus_class") or "", genes]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(REPO, "web/intro/data"))
    out = ap.parse_args().out
    os.makedirs(out, exist_ok=True)
    shards = {}
    for key, folder in CAMPAIGNS:
        for name, rep in reports(os.path.join(R, folder, "reports_all.tar.zst")):
            loci = [locus_row(l) for l in rep.get("detected") or []]
            if not loci:
                continue
            if key == "mucoromycotina":
                src, _, genome = name.partition("__")
                shards.setdefault(MUCORO[src], {})[genome] = loci
            else:
                shards.setdefault(key, {})[name] = loci
    for ds, per in sorted(shards.items()):
        data = json.dumps({"fields": FIELDS, "loci": dict(sorted(per.items()))}, separators=(",", ":")).encode()
        path = os.path.join(out, f"loci_{ds}.json.gz")
        with open(path, "wb") as fh:  # mtime=0 keeps the file byte-identical between runs
            fh.write(gzip.compress(data, 9, mtime=0))
        print(f"{ds}: {len(per)} genomes, {sum(len(v) for v in per.values())} loci, {os.path.getsize(path) / 1e6:.2f} MB")


if __name__ == "__main__":
    main()
