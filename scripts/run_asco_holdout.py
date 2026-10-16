#!/usr/bin/env python
"""Held-out test of the Ascomycota idiomorph call and the alpha-box / HMG-box profiles, one record and radius.

    run_asco_holdout.py run --manifest testset/holdout/ascomycota_holdout.tsv --db db --record RID --radius order \
        --library DIR --work DIR --out row.json
    run_asco_holdout.py summarize ROW.json [...]

`run`: withholds the record's taxon at the radius (detect.holdout.records_to_withhold), runs `detect --taxid` on the
record's own genome with those records excluded, and scores the call against the record's mating type: the
pipeline's idiomorph for a `MAT` locus, and the profiles' call. The profiles (MAT1-1-1 alpha box, MAT1-2-1 HMG box)
are rebuilt WITHOUT the withheld records, so the profile is as blind as the search. Each core-gene protein of the
locus is scored on both profiles; the locus is called by the profiles when the alpha protein's alpha score and the
HMG protein's HMG score differ by at least MIN_MARGIN bits (one core gene: which profile it scores higher on).
"""
import argparse
import csv
import json
import os
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

import yaml

MIN_MARGIN = 25.0
GENES = {"MAT1-1-1": "MAT1-1", "MAT1-2-1": "MAT1-2"}
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def coding_sequence(exons, strand, seq):
    dna = "".join(seq[a - 1:b] for a, b in sorted(exons))
    return dna.translate(_COMP)[::-1] if strand == "-" else dna


def protein_candidates(entry, seq):
    from Bio.Seq import Seq
    if entry.get("exons"):
        ex = tuple((int(x["start"]), int(x["end"])) for x in entry["exons"])
        dna = coding_sequence(ex, entry["strand"], seq)
        return [str(Seq(dna[: len(dna) // 3 * 3]).translate()).replace("*", "")]
    dna = coding_sequence(((int(entry["start"]), int(entry["end"])),), entry["strand"], seq)
    return [str(Seq(dna[f:][: (len(dna) - f) // 3 * 3]).translate()).replace("*", "X") for f in range(3)]


def build_profiles(db: Path, drop: set[str], work: Path):
    """Alpha-box and HMG-box HMMs from the curated records not in `drop` (same builder as the classifiers)."""
    from MATPredict.detect import classifier_build as cb
    from MATPredict.detect.family_registry import FamilyKey, load_all_families

    family = next(f for f in load_all_families(db) if f.key == FamilyKey("Ascomycota", "MAT"))
    rows = [r for r in cb.training_set(db, family, work) if r["gene"] in GENES
            and r["id"].split("|")[1] not in drop]
    hmms, n = {}, {}
    for gene in GENES:
        seqs = {r["id"]: r["sequence"] for r in rows if r["gene"] == gene}
        n[gene] = len(seqs)
        if not seqs:
            return None, n
        hmms[gene], _ = cb.build_hmm(seqs, gene, work)
    return hmms, n


def call(a, h):
    if a - h >= MIN_MARGIN: return "MAT1-1"
    if h - a >= MIN_MARGIN: return "MAT1-2"
    return "undetermined"


def run(a) -> int:
    from Bio import SeqIO
    from MATPredict.detect.classifier_build import score
    from MATPredict.detect.holdout import Radius, records_to_withhold

    rows = {r["record_id"]: r for r in csv.DictReader(open(a.manifest), delimiter="\t")}
    r = rows[a.record]
    drop = set(records_to_withhold(a.db, a.record, Radius(a.radius)))
    work = Path(a.work) / a.record / a.radius
    work.mkdir(parents=True, exist_ok=True)
    asm = next(Path(a.library).glob(f"{r['assembly']}_*.fa.gz"), None)
    out = dict(record_id=a.record, radius=a.radius, subphylum=r["subphylum"], family=r["family"], system=r["system"],
               expected=r["idiomorphs"], species=r["species"], withheld=len(drop))
    if asm is None:
        out["status"] = "no_genome_in_library"
        Path(a.out).write_text(json.dumps(out)); return 0
    fna = work / f"{r['assembly']}.fna"
    if not fna.exists():
        subprocess.run(f"zcat {asm} > {fna}", shell=True, check=True)
    det = work / "detect"
    cmd = [sys.executable, "-m", "MATPredict", "detect", "--genome", str(fna), "--out-dir", str(det),
           "--exclude-records", ",".join(sorted(drop)), "--taxid", str(r["taxid"])]
    env = dict(os.environ, MATPREDICT_DB_ROOT=str(a.db), MATPREDICT_HTML="0")
    p = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=7200)
    if p.returncode != 0:
        out.update(status="detect_error", stderr=p.stderr[-300:]); Path(a.out).write_text(json.dumps(out)); return 0
    doc = yaml.safe_load((det / "detection_report.yaml").read_text())
    loci = [c for c in doc.get("detected") or [] if c.get("family") == "Ascomycota:MAT"]
    out["loci"] = [dict(contig=c["contig"], start=c["start"], end=c["end"], idiomorph=c["idiomorph"],
                        confidence=c["confidence"], locus_class=c["locus_class"], genes=c["genes_found"])
                   for c in loci]
    out["pipeline_labels"] = sorted({c["idiomorph"] for c in loci})
    hmms, n = build_profiles(a.db, drop, work)
    out["profile_training"] = n
    if hmms is None:
        out["status"] = "no_profile_training"
    else:
        wanted = {c["contig"] for c in loci}
        seqs = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(str(fna), "fasta") if rec.id in wanted}
        prof = []
        for c in loci:
            ent = {}
            for g in c["gene_evidence"]:
                if g["gene"] in GENES and not g.get("lost_to"):
                    if g["gene"] not in ent or (g.get("bitscore") or 0) > (ent[g["gene"]].get("bitscore") or 0):
                        ent[g["gene"]] = g
            sc = {}
            for gene, e in ent.items():
                cands = protein_candidates(e, seqs[c["contig"]])
                sc[gene] = {pg: max(score(h, {f"p{i}": x for i, x in enumerate(cands) if x}).values(), default=0.0)
                            for pg, h in hmms.items()}
            a_on_a = sc.get("MAT1-1-1", {}).get("MAT1-1-1")
            h_on_h = sc.get("MAT1-2-1", {}).get("MAT1-2-1")
            if a_on_a is not None and h_on_h is not None:
                pc = call(a_on_a, h_on_h)
            elif sc:
                only = next(iter(sc.values())); pc = call(only["MAT1-1-1"], only["MAT1-2-1"])
            else:
                pc = "no_core_gene"
            prof.append(dict(idiomorph_pipeline=c["idiomorph"], alpha_on_alpha=a_on_a, hmg_on_hmg=h_on_h,
                             alpha_on_hmg=sc.get("MAT1-1-1", {}).get("MAT1-2-1"),
                             hmg_on_alpha=sc.get("MAT1-2-1", {}).get("MAT1-1-1"), profile_call=pc,
                             core_genes_with_protein=sorted(sc)))
        out["profile"] = prof
        out["status"] = "ok"
    Path(a.out).write_text(json.dumps(out, default=str, indent=1))
    return 0


def summarize(a) -> int:
    rows = [json.loads(Path(p).read_text()) for p in a.rows]
    print(f"{'record':36s} {'radius':7s} {'system':13s} {'expected':14s} {'wh':>3s} loci  pipeline        profile")
    for r in sorted(rows, key=lambda x: (x["record_id"], ["genus", "family", "order"].index(x["radius"]))):
        exp = set(r["expected"].split("|"))
        pl = ",".join(r.get("pipeline_labels") or []) or "-"
        pf = ",".join(p["profile_call"] for p in r.get("profile") or []) or r["status"]
        ok_p = (bool(exp & set(r.get("pipeline_labels") or [])) if r["system"] != "homothallic"
                else bool(exp & set(r.get("pipeline_labels") or [])))
        print(f"{r['record_id']:36s} {r['radius']:7s} {r['system']:13s} {r['expected']:14s} {r['withheld']:3d} "
              f"{len(r.get('loci') or []):4d}  {pl:15s}  {pf}  {'' if r['status']=='ok' else r['status']}")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--manifest", required=True); r.add_argument("--db", type=Path, required=True)
    r.add_argument("--record", required=True); r.add_argument("--radius", required=True)
    r.add_argument("--library", required=True); r.add_argument("--work", required=True); r.add_argument("--out", required=True)
    s = sub.add_parser("summarize"); s.add_argument("rows", nargs="+")
    a = ap.parse_args(argv)
    return run(a) if a.cmd == "run" else summarize(a)


if __name__ == "__main__":
    sys.exit(main())
