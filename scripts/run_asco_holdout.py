#!/usr/bin/env python
"""Held-out test of the locus and idiomorph call (any phylum) and, for Ascomycota MAT, the alpha/HMG profiles; one record and radius.

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
    pre = r.get(f"{a.radius}_drop")  # a manifest row may carry its own drop list (genomes that are not records)
    drop = set(pre.split(",")) - {""} if pre is not None and f"{a.radius}_drop" in r else \
        set(records_to_withhold(a.db, a.record, Radius(a.radius)))
    work = Path(a.work) / a.record / a.radius
    work.mkdir(parents=True, exist_ok=True)
    asm = next(Path(a.library).glob(f"{r['assembly']}_*.fa.gz"), None)
    fam_key = f"{r['phylum']}:{r['family']}"
    out = dict(record_id=a.record, radius=a.radius, phylum=r["phylum"], subphylum=r["subphylum"], family=fam_key,
               system=r["system"], expected=r["idiomorphs"], species=r["species"], withheld=len(drop))
    out["ref_family_left"] = int(r[f"{a.radius}_same_family_left"])
    out["ref_idiomorph_left"] = int(r[f"{a.radius}_same_idiomorph_left"])
    out["classifier_leak"] = fam_key == "Mucoromycota:MAT"  # prebuilt idiomorph classifier saw every record
    if out["ref_family_left"] == 0:  # nothing of this family is left to search with: a database gap, not a miss
        out["status"] = "no_reference_family"
        Path(a.out).write_text(json.dumps(out)); return 0
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
    from MATPredict.detect.family_registry import FamilyKey, load_all_families
    fam = next((f for f in load_all_families(a.db) if f.key == FamilyKey(r["phylum"], r["family"])), None)
    out["enum_family"] = bool(fam and fam.vocabulary_type == "enum")
    out["routing_mode"] = doc.get("routing_mode")
    out["all_loci_families"] = sorted({c.get("family") for c in doc.get("detected") or []})
    loci = [c for c in doc.get("detected") or [] if c.get("family") == fam_key]
    spans = json.loads(r.get("truth_spans") or "{}")
    contigs = {rec.id for rec in SeqIO.parse(str(fna), "fasta")} if spans else set()
    comparable = {k: v for k, v in spans.items() if k in contigs}
    out["truth_comparable"] = bool(comparable)
    out["found_at_truth"] = (any(c["contig"] in comparable and c["start"] <= comparable[c["contig"]][1]
                                 and c["end"] >= comparable[c["contig"]][0] for c in loci) if comparable else None)
    out["loci"] = [dict(contig=c["contig"], start=c["start"], end=c["end"], idiomorph=c["idiomorph"],
                        confidence=c["confidence"], locus_class=c["locus_class"], genes=c["genes_found"])
                   for c in loci]
    out["pipeline_labels"] = sorted({c["idiomorph"] for c in loci})
    if fam_key != "Ascomycota:MAT":
        out["status"] = "ok"
        Path(a.out).write_text(json.dumps(out, default=str, indent=1)); return 0
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
    if a.manifest:  # recompute the reference counts from the current manifest (same phylum and family only)
        man = {m["record_id"]: m for m in csv.DictReader(open(a.manifest), delimiter="\t")}
        for r in rows:
            m = man[r["record_id"]]
            r["ref_family_left"] = int(m[f"{r['radius']}_same_family_left"])
            r["ref_idiomorph_left"] = int(m[f"{r['radius']}_same_idiomorph_left"])
    order = ["genus", "family", "order"]
    print(f"{'record':34s} {'rad':6s} {'system':13s} {'expected':12s} {'wh':>3s} loci at-truth pipeline            profile")
    tally = {}
    for r in sorted(rows, key=lambda x: (x.get("phylum", ""), x["family"], x["record_id"], order.index(x["radius"]))):
        exp = set(r["expected"].split("|")) - {""}
        labels = r.get("pipeline_labels") or []
        found = bool(r.get("loci"))
        if r["status"] in ("detect_error", "no_genome_in_library", "no_reference_family"):
            verdict = r["status"]
        elif not found:
            verdict = "no_reference_family" if r.get("ref_family_left") == 0 else "no_locus"
        elif not r.get("enum_family"):
            verdict = "locus_found"
        elif exp & set(labels):
            verdict = "correct"
        elif labels and set(labels) <= {"undetermined"}:
            verdict = "undetermined"
        else:
            verdict = "wrong"
        if verdict in ("no_locus", "wrong") and r.get("enum_family") and r.get("ref_idiomorph_left") == 0:
            verdict = "no_reference_idiomorph"
        r["verdict"] = verdict
        pf = ",".join(p["profile_call"] for p in r.get("profile") or []) or "-"
        at = {None: "n/a", True: "yes", False: "NO"}[r.get("found_at_truth")]
        print(f"{r['record_id'][:34]:34s} {r['radius']:6s} {r['system'][:13]:13s} {r['expected'][:12]:12s} {r['withheld']:3d} "
              f"{len(r.get('loci') or []):4d} {at:5s}   {','.join(labels)[:16] or '-':16s} {verdict:14s} {pf}")
    print()
    for key in sorted({(r.get("phylum", ""), r.get("subphylum", ""), r["family"], r["radius"]) for r in rows}):
        sel = [r for r in rows if (r.get("phylum", ""), r.get("subphylum", ""), r["family"], r["radius"]) == key]
        from collections import Counter as C
        c = C(r["verdict"] for r in sel)
        print(f"{key[2]:24s} {key[1]:18s} {key[3]:6s} n={len(sel):3d} {dict(c)}")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--manifest", required=True); r.add_argument("--db", type=Path, required=True)
    r.add_argument("--record", required=True); r.add_argument("--radius", required=True)
    r.add_argument("--library", required=True); r.add_argument("--work", required=True); r.add_argument("--out", required=True)
    s = sub.add_parser("summarize"); s.add_argument("rows", nargs="+")
    s.add_argument("--manifest", help="recompute the references-left counts from this manifest")
    a = ap.parse_args(argv)
    return run(a) if a.cmd == "run" else summarize(a)


if __name__ == "__main__":
    sys.exit(main())
