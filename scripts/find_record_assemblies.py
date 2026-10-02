#!/usr/bin/env python3
"""Find a same-strain genome assembly for curated records that lack one (B10).

Curator ruling 2026-10-01: for each record with no `locus.assembly_accession`
(most are locus-specific GenBank deposits that belong to no assembly), search
NCBI for an assembly of the SAME STRAIN; where none is found the self-check
reports "not possible", not a failure.

A strain-name match alone is not trusted. Three steps:

  search  NCBI Assembly by the record's taxid; keep assemblies whose strain or
          isolate equals the record's strain or one of its culture-collection
          IDs (compared without spaces, dashes or case). Writes candidates.tsv.
  verify  For candidates in the BFD library, blastn the record's locus
          sequence (locus.gbk) against the genome. A candidate passes only at
          >= 99% identity over >= 90% of every core segment, and the best
          hit's contig and span become the record's location in that
          assembly. Writes verified.tsv. Run it under SLURM (it decompresses
          genomes into $SCRATCH).
  write   Insert `assembly_accession`, `assembly_accession_basis:
          strain_match_verified` and `assembly_location` into each verified
          record. Nothing else in the record changes.

Usage:
  find_record_assemblies.py search --db DB --out candidates.tsv
  find_record_assemblies.py verify --db DB --candidates candidates.tsv --out verified.tsv
  find_record_assemblies.py write  --db DB --verified verified.tsv
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import re
import subprocess
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

import yaml

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
LIBRARY = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes"
MIN_IDENTITY = 99.0
MIN_SEGMENT_COVERAGE = 90.0


def _get(url, params, tries=4):
    q = urllib.parse.urlencode({**params, "tool": "MATPredict"})
    for i in range(tries):
        try:
            with urllib.request.urlopen(f"{EUTILS}/{url}?{q}", timeout=60) as r:
                return r.read()
        except Exception:  # noqa: BLE001 - retried, then raised
            if i == tries - 1:
                raise
            time.sleep(2 * (i + 1))
    return b""


def norm(s):
    return re.sub(r"[^A-Z0-9]", "", (s or "").upper())


def record_strains(doc):
    o = doc.get("organism") or {}
    st = o.get("strain")
    names = []
    if isinstance(st, dict):
        names.append(st.get("name"))
        names += st.get("culture_collection_ids") or []
    else:
        names.append(st)
    return sorted({norm(n) for n in names if n and norm(n) not in ("", "UNKNOWN")})


def records_without_assembly(db):
    """Records the self-check cannot place: no `locus.assembly_accession` and
    no assembly found by record_selfcall's own fallbacks (an assembly-type
    segment, or an assembly accession in the record text). Run
    scripts/backfill_record_assembly.py (sequence links) first."""
    from MATPredict.detect.record_selfcall import record_location
    for md in sorted(Path(db).glob("*/*/*/metadata.yaml")):
        if not (md.parent / "locus.gbk").exists():
            continue
        doc = yaml.safe_load(md.read_text())
        if record_location(md).assembly:
            continue
        yield md, doc


def assemblies_for_taxid(taxid):
    """[(accession, strain/isolate strings, organism)] for every assembly under taxid."""
    root = ET.fromstring(_get("esearch.fcgi", {"db": "assembly", "term": f"txid{taxid}[Organism:exp]",
                                               "retmax": 2000}))
    ids = [e.text for e in root.findall(".//IdList/Id")]
    time.sleep(0.4)
    out = []
    for i in range(0, len(ids), 200):
        root = ET.fromstring(_get("esummary.fcgi", {"db": "assembly", "id": ",".join(ids[i:i + 200])}))
        time.sleep(0.4)
        for d in root.iter("DocumentSummary"):
            acc = d.findtext("AssemblyAccession") or ""
            vals = [e.text for e in d.iter("Sub_value") if e.text]
            vals += [d.findtext("Biosource/Isolate") or ""]
            out.append((acc, vals, d.findtext("Organism") or ""))
    return out


def cmd_search(a):
    rows = []
    cache = {}
    for md, doc in records_without_assembly(a.db):
        rid = doc["record_id"]
        taxid = (doc.get("taxonomy") or {}).get("taxid")
        want = record_strains(doc)
        if not taxid or not want:
            rows.append(dict(record_id=rid, status="no_taxid_or_strain", strains=";".join(want)))
            continue
        if taxid not in cache:
            try:
                cache[taxid] = assemblies_for_taxid(taxid)
            except Exception as e:  # noqa: BLE001
                rows.append(dict(record_id=rid, status=f"ncbi_error:{e}", strains=";".join(want)))
                continue
        hits = [(acc, org) for acc, vals, org in cache[taxid] if any(norm(v) in want for v in vals)]
        if not hits:
            rows.append(dict(record_id=rid, status="no_same_strain_assembly", strains=";".join(want),
                             n_taxid_assemblies=len(cache[taxid])))
            continue
        for acc, org in hits:
            lib = sorted(glob.glob(f"{LIBRARY}/{acc}_*.fa.gz")) or sorted(glob.glob(f"{LIBRARY}/{acc}.fa.gz"))
            rows.append(dict(record_id=rid, status="candidate" if lib else "candidate_not_in_bfd",
                             strains=";".join(want), assembly=acc, organism=org,
                             bfd_file=lib[0] if lib else "", n_taxid_assemblies=len(cache[taxid])))
        print(f"{rid}: {len(hits)} same-strain assemblies", file=sys.stderr)
    fields = ["record_id", "status", "strains", "assembly", "organism", "bfd_file", "n_taxid_assemblies"]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def segment_lengths(gbk):
    """{record id: length} for the segment records in locus.gbk. Each record
    holds only its segment's sequence, so coverage is over its full length."""
    from Bio import SeqIO
    return {r.id: len(r) for r in SeqIO.parse(str(gbk), "genbank")}


def gbk_to_fasta(gbk, out):
    from Bio import SeqIO
    recs = list(SeqIO.parse(str(gbk), "genbank"))
    with open(out, "w") as fh:
        for r in recs:
            fh.write(f">{r.id}\n{str(r.seq)}\n")
    return [r.id for r in recs]


def cmd_verify(a):
    scratch = os.environ.get("SCRATCH")
    if not scratch:
        sys.exit("SCRATCH is not set -- run verify under SLURM")
    by_rid = {}
    for md in Path(a.db).glob("*/*/*/metadata.yaml"):
        by_rid[md.parent.name] = md
    rows = []
    for c in csv.DictReader(open(a.candidates), delimiter="\t"):
        if c["status"] != "candidate":
            continue
        md = by_rid.get(c["record_id"])
        if md is None:
            continue
        q = Path(scratch) / f"{c['record_id']}.locus.fa"
        gbk_to_fasta(md.parent / "locus.gbk", q)
        g = Path(scratch) / (Path(c["bfd_file"]).name[:-3])
        if not g.exists():
            subprocess.run(f"zcat {c['bfd_file']} > {g}", shell=True, check=True)
        res = subprocess.run(["blastn", "-task", "megablast", "-query", str(q), "-subject", str(g),
                              "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore",
                              "-evalue", "1e-20"], capture_output=True, text=True, check=True).stdout
        hsps = [l.split("\t") for l in res.splitlines()]
        qlen = segment_lengths(md.parent / "locus.gbk")
        seg_cov, best = [], None
        for qid, L in qlen.items():
            # Coverage of this segment by HSPs at >= MIN_IDENTITY on one contig.
            per_contig = {}
            for h in hsps:
                if h[0] == qid and float(h[2]) >= MIN_IDENTITY:
                    per_contig.setdefault(h[1], []).append(h)
            # Group each contig's HSPs into runs no further apart than the
            # segment length, so repeat copies elsewhere on a chromosome do not
            # stretch the span (a 1.4 Mb "locus" in Cryptococcus JEC20 before
            # this). Keep the run that covers most of the segment.
            covs, groups = {}, {}
            for contig, hs in per_contig.items():
                hs = sorted(hs, key=lambda h: min(int(h[6]), int(h[7])))
                runs, cur, cur_end = [], [], 0
                for h in hs:
                    a0, a1 = sorted((int(h[6]), int(h[7])))
                    if cur and a0 - cur_end > L:
                        runs.append(cur)
                        cur, cur_end = [], 0
                    cur.append(h)
                    cur_end = max(cur_end, a1)
                runs.append(cur)
                for run in runs:
                    cov = set()
                    for h in run:
                        qs, qe = sorted((int(h[4]), int(h[5])))
                        cov.update(range(qs, qe + 1))
                    run_cov = 100.0 * len(cov) / L
                    if run_cov > covs.get(contig, -1):
                        covs[contig], groups[contig] = run_cov, run
            if not covs:
                seg_cov.append(0.0)
                continue
            contig = max(covs, key=covs.get)
            seg_cov.append(covs[contig])
            hs = groups[contig]
            ss = [int(v) for h in hs for v in (h[6], h[7])]
            pid = sum(float(h[2]) * int(h[3]) for h in hs) / sum(int(h[3]) for h in hs)
            best = best or (contig, min(ss), max(ss), round(pid, 2))
        segs = list(qlen)
        span_ok = bool(best) and (best[2] - best[1] + 1) <= 1.5 * sum(qlen.values())
        ok = bool(segs) and all(v >= MIN_SEGMENT_COVERAGE for v in seg_cov) and span_ok
        rows.append(dict(record_id=c["record_id"], assembly=c["assembly"],
                         verdict="verified" if ok else ("span_too_long" if best and not span_ok else "not_verified"),
                         segment_coverage=";".join(f"{v:.1f}" for v in seg_cov),
                         contig=best[0] if best else "", start=best[1] if best else "",
                         end=best[2] if best else "", identity=best[3] if best else ""))
        print(f"{c['record_id']} {c['assembly']}: {rows[-1]['verdict']} {rows[-1]['segment_coverage']}",
              file=sys.stderr)
        if not a.keep_genomes:
            g.unlink(missing_ok=True)
    fields = ["record_id", "assembly", "verdict", "segment_coverage", "contig", "start", "end", "identity"]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def cmd_write(a):
    by_rid = {md.parent.name: md for md in Path(a.db).glob("*/*/*/metadata.yaml")}
    ver = {}
    for r in csv.DictReader(open(a.verified), delimiter="\t"):
        if r["verdict"] == "verified":
            ver.setdefault(r["record_id"], []).append(r)
    for rid, rs in sorted(ver.items()):
        if len({r["assembly"].split(".")[0] for r in rs}) != 1:
            print(f"{rid}: {len(rs)} different verified assemblies; not written", file=sys.stderr)
            continue
        r = max(rs, key=lambda x: int(x["assembly"].split(".")[1]))
        md = by_rid[rid]
        lines = md.read_text().splitlines(keepends=True)
        if any(l.startswith("  assembly_accession:") for l in lines):
            continue
        idx = next(i for i, ln in enumerate(lines) if ln.rstrip("\n") == "locus:")
        lines[idx + 1:idx + 1] = [
            f"  assembly_accession: {r['assembly']}\n",
            "  assembly_accession_basis: strain_match_verified\n",
            "  assembly_location:\n",
            f"    contig: {r['contig']}\n",
            f"    start: {r['start']}\n",
            f"    end: {r['end']}\n",
            f"    identity: {r['identity']}\n",
        ]
        md.write_text("".join(lines))
        print(f"{rid}: wrote {r['assembly']} {r['contig']}:{r['start']}-{r['end']}", file=sys.stderr)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("search"); s.add_argument("--db", required=True); s.add_argument("--out", required=True)
    v = sub.add_parser("verify"); v.add_argument("--db", required=True); v.add_argument("--candidates", required=True)
    v.add_argument("--out", required=True); v.add_argument("--keep-genomes", action="store_true")
    w = sub.add_parser("write"); w.add_argument("--db", required=True); w.add_argument("--verified", required=True)
    a = ap.parse_args(argv)
    {"search": cmd_search, "verify": cmd_verify, "write": cmd_write}[a.cmd](a)


if __name__ == "__main__":
    main()
