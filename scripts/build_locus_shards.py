#!/usr/bin/env python3
"""Compact per-locus extracts from the campaign reports, for the website's drill-down view.

Reads each campaign's reports_all.tar.zst (one detection_report.yaml per genome, at any depth: v0.6.0
`wave_N/runs/<genome>/`, clean-run `wave_N/small_NNN/runs/<genome>/` or `runs/small_NNN/runs/<genome>/`,
Mucoromycotina `<SOURCE>__<id>/`) and keeps, for every detected locus, the contig, span, idiomorph,
confidence, locus class and the gene models (name, span, strand, role, status; and, for a gene that took part in a resolved
idiomorph pair, a seventh item [W|L, the other gene, the basis, the scores]: W won the call, L was set aside). When the campaign folder also
holds the clean-run `loci.tsv` (scripts/build_cleanrun_tables.py), each locus gets five more fields after the
gene list: baseline_confidence, tier_note, label_unsupported_reason, high_basis and defective_genes (the viewer reads fields 0-6
and ignores the rest). A genome whose loci.tsv rows do not match its report's loci in number stops the build.
Writes web/intro/data/loci_<dataset>.json.gz (committed; about 1.9 MB for all 27.8k loci). Run this when a
campaign changes; build_site_data.py and the site build do not need the tarballs.

The campaign folders are the release's, from RELEASES in build_site_data.py (--release, default as there);
--asco / --basidio / --mucoro override one, as in build_site_data.py. If a genome has several run folders in
a tarball, the run is chosen by `choose_run` from build_cleanrun_tables.py (a non-empty report, then the last
path), so the maps and the tables describe the same run.

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
from build_cleanrun_tables import choose_run  # noqa: E402
from build_site_data import DEFAULT_RELEASE, RELEASES, REPO, resolve_release  # noqa: E402

MUCORO = {"BFD": "bfd_mucoro", "LCG": "zygolife_lcg", "JENA": "zygolife_jena"}
FIELDS = ["contig", "start", "end", "idiomorph", "confidence", "locus_class",
          "genes[name,start,end,strand,role,status,(verdict: W|L, other gene, basis, scores)]"]
NOTE_FIELDS = ["baseline_confidence", "tier_note", "label_unsupported_reason", "high_basis", "defective_genes"]


def reports(tar_path):
    """Yield (member path, genome directory name, report bytes non-empty, parsed detection_report.yaml) from a
    .tar.zst without unpacking it."""
    proc = subprocess.Popen(["zstd", "-dc", tar_path], stdout=subprocess.PIPE)
    with tarfile.open(fileobj=proc.stdout, mode="r|") as tf:
        for m in tf:
            if m.isfile() and m.name.endswith("/detection_report.yaml"):
                data = tf.extractfile(m).read()
                yield (os.path.normpath(os.path.dirname(m.name)), os.path.basename(os.path.dirname(m.name)), bool(data),
                       (yaml.load(data, Loader=getattr(yaml, "CSafeLoader", yaml.SafeLoader)) or {}) if data else {})
    if proc.wait():
        raise RuntimeError(f"zstd failed on {tar_path}")


def _score(r):
    """The verdict's scores as text, winner first: the HMM bitscores when it used them, else the alignment identities; or ""."""
    w, lo = r.get("winner_hmm_score"), r.get("loser_hmm_score")
    if w is not None and lo is not None:
        return f"{w:.0f} vs {lo:.0f}"
    w, lo = r.get("winner_identity"), r.get("loser_identity")
    return f"{w:.1f}% vs {lo:.1f}%" if w is not None and lo is not None else ""


def verdicts(l):
    """{gene entry index: ["W"|"L", other gene, basis, "winner vs loser score"]} for the genes of a locus that took part in a
    resolved idiomorph pair: W won, L was set aside (the report's gene `lost_to`, `idiomorph_resolutions`). The pair is matched by
    gene names and, when the report gives them, by the gene's own span. A gene with no verdict gets no entry."""
    genes, out = l.get("gene_evidence") or [], {}
    for r in l.get("idiomorph_resolutions") or []:
        win, lose = r.get("winner"), r.get("loser")
        basis, score = r.get("basis") or "", _score(r)
        for role, name, other, span in (("W", win, lose, r.get("winner_span")), ("L", lose, win, r.get("loser_span"))):
            for i, g in enumerate(genes):
                if g.get("gene") != name or i in out or (span and [g.get("start"), g.get("end")] != list(span)):
                    continue
                if role == "L" and g.get("lost_to") not in (None, other):
                    continue
                out[i] = [role, other, basis, score]
                break
    return out


def locus_row(l):
    verdict = verdicts(l)
    genes = [[g["gene"], g["start"], g["end"], g["strand"], g["role"], g["status"]] + ([verdict[i]] if i in verdict else [])
             for i, g in enumerate(l.get("gene_evidence") or [])]
    return [l["contig"], l["start"], l["end"], l.get("idiomorph") or "undetermined", l.get("confidence") or "", l.get("locus_class") or "", genes]


def load_notes(folder):
    """{genome (SOURCE__id for Mucoromycotina): [NOTE_FIELDS values per locus]} from the campaign's clean-run loci.tsv,
    or None when it has no such columns. A column the table predates (defective_genes, before #151) reads as blank."""
    path = os.path.join(folder, "loci.tsv")
    if not os.path.exists(path):
        return None
    out = defaultdict(list)
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        if "label_unsupported_reason" not in (rd.fieldnames or []):
            return None
        for r in rd:
            g = f"{r['source']}__{r['genome']}" if r.get("source") else r["genome"]
            out[g].append([r.get(c, "") or "" for c in NOTE_FIELDS])
    return out


def campaign_shards(key, folder):
    """{dataset: {genome: [locus rows]}} for one campaign folder; the second value says whether notes were added."""
    notes = load_notes(folder)
    runs = {}
    for path, name, nonempty, rep in reports(os.path.join(folder, "reports_all.tar.zst")):
        runs.setdefault(name, {})[path] = (nonempty, [locus_row(l) for l in rep.get("detected") or []])
    best, bad = {}, []
    for name, cands in runs.items():
        loci = cands[choose_run([(p, ne) for p, (ne, _) in cands.items()])][1]
        if not loci:
            continue
        if notes is not None:
            if len(notes.get(name, [])) != len(loci):
                bad.append(name)
                continue
            loci = [l + n for l, n in zip(loci, notes[name])]   # loci.tsv rows are in report order
        best[name] = loci
    if bad:
        raise SystemExit(f"{folder}: loci.tsv and the reports disagree on the number of loci for {len(bad)} genome(s), "
                         f"e.g. {sorted(bad)[:5]}; rebuild loci.tsv from the same reports")
    shards = {}
    for name, loci in best.items():
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
