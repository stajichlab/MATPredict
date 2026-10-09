#!/usr/bin/env python3
"""Validation of the report-only supported_span over the 34-genome panel at three bitscore floors.
usage: validate.py   (run in this folder; needs floor30/ floor39/ floor50/ and ../2026-10-08_core_span_panel/newdb)"""
import csv, glob, re
import yaml

FLOORS = (30, 33, 36, 39, 50)
import sys
sys.path.insert(0, "/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-supported2/src")
from pathlib import Path
from MATPredict.detect.family_registry import load_all_families
FAM_GENES = {f"{f.key.phylum}:{f.key.locus_name}": {g["name"] for g in f.genes} | {a for g in f.genes for a in (g.get("aliases") or [])}
             for f in load_all_families(Path("/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-supported2/db"))}
BASE = "../2026-10-08_core_span_panel/newdb"      # earlier run of the core-span code (no supported_span), same database

def load(d):
    out = {}
    for f in glob.glob(f"{d}/runs/*/detection_report.yaml"):
        out[f.split("/")[-2]] = yaml.safe_load(open(f))
    return out

def callkey(x):
    return (x["family"], x["contig"], x["start"], x["end"], x["confidence"], x["idiomorph"], x.get("locus_class"),
            tuple(sorted(x.get("genes_found") or [])), tuple(sorted(str(m) for m in (x.get("merged_from") or []))))

def sup_key(r):   # everything in a report except the new field
    return ([callkey(x) for x in r.get("detected") or []],
            sorted((s["family"], s["contig"], s["start"], s["end"], s.get("withheld_reason")) for s in r.get("suppressed_loci") or []))

R = {f: load(f"floor{f}") for f in FLOORS}
base = load(BASE)
print("genomes per floor:", {f: len(R[f]) for f in FLOORS}, "| earlier core-span run:", len(base))

print("\n1. CALL IDENTITY (cluster span, calls, genes, merged_from, withheld loci)")
for f in FLOORS:
    common = [g for g in R[f] if g in base]
    diff = [g for g in common if sup_key(R[f][g]) != sup_key(base[g])]
    print(f"   floor {f} vs earlier run: {len(common)} common genomes, {len(diff)} with any difference {diff[:3]}")
for a, b in ((30, 33), (33, 36), (36, 39), (39, 50)):
    diff = [g for g in R[a] if g in R[b] and sup_key(R[a][g]) != sup_key(R[b][g])]
    print(f"   floor {a} vs floor {b}: {len(diff)} genomes differ in anything but supported_span")

print("\n2. GENE CONTAINMENT (every reported gene on the locus contig lies inside core, supported and cluster span)")
for f in FLOORS:
    n = bad = 0
    for g, r in R[f].items():
        for x in r.get("detected") or []:
            sp = x.get("supported_span")
            if not sp: continue
            for e in x.get("gene_evidence") or []:
                if e.get("contig") != x["contig"]: continue
                n += 1
                if not (sp["start"] <= e["start"] and e["end"] <= sp["end"] and x["start"] <= e["start"] and e["end"] <= x["end"]): bad += 1
    print(f"   floor {f}: {n} genes checked, {bad} outside")

print("\n3. SWEEP over called loci")
rows = {}
for f in FLOORS:
    loci = [(g, x) for g, r in R[f].items() for x in r.get("detected") or []]
    sb = [x["supported_span"]["beyond_supported_bp"] for _, x in loci if x.get("supported_span")]
    cb = [x["core_span"]["beyond_core_bp"] for _, x in loci if x.get("core_span")]
    rows[f] = loci
    print(f"   floor {f}: loci {len(loci)}; beyond_supported > 0: {sum(1 for v in sb if v>0)}, >= 1 kb: {sum(1 for v in sb if v>=1000)}, >= 10 kb: {sum(1 for v in sb if v>=10000)}, "
          f"total {sum(sb)} bp, max {max(sb) if sb else 0} | for comparison beyond_core: > 0 {sum(1 for v in cb if v>0)}, >= 10 kb {sum(1 for v in cb if v>=10000)}, total {sum(cb)} bp")

print("\n4. THE 14 HAND-TALLIED STRETCHES (noise-only should be dropped; strong OWN-family hits kept; other-family hits dropped by design)")
tally = [r for r in csv.DictReader(open("../2026-10-08_core_span_panel/tally_newdb.tsv"), delimiter="\t") if int(r["stretch_bp"]) >= 1000]
def find(f, genome_prefix, contig, family):
    for g, r in R[f].items():
        if g.startswith(genome_prefix):
            for x in r.get("detected") or []:
                if x["contig"] == contig and x["family"].split(":")[1] == family: return x
for f in FLOORS:
    c = dict(own_kept=0, own_lost=0, other_kept=0, other_dropped=0, noise_dropped=0, noise_kept=0); detail = []
    for t in tally:
        x = find(f, t["genome"][:30], t["contig"], t["family"])
        if not x or not x.get("supported_span"): continue
        sp = x["supported_span"]; stretch_bp = int(t["stretch_bp"]); strong = int(t["strong"])
        bh = t["best_hit_in_stretch(query;pos;e;aa;pid)"].split(";")
        pos = [int(v) for v in bh[1].split("-")] if len(bh) > 1 and "-" in bh[1] else None
        if strong > 0 and pos:
            gene = bh[0].split("|")[-1]
            own = gene in FAM_GENES.get(x["family"], set())
            inside = sp["start"] <= pos[0] and pos[1] <= sp["end"]
            key = ("own_" if own else "other_") + ("kept" if inside else ("lost" if own else "dropped"))
            c[key] += 1
            if own and not inside: detail.append(("OWN-family strong hit LOST", t["genome"][:18], t["family"], t["side"], bh[0][-30:], bh[2], f"aa {bh[3]}"))
        else:
            trimmed = (t["side"] == "left" and sp["start"] > int(t["core"].split("-")[0]) - stretch_bp * 0.5) or (t["side"] == "right" and sp["end"] < int(t["core"].split("-")[1]) + stretch_bp * 0.5)
            c["noise_dropped" if trimmed else "noise_kept"] += 1
            if not trimmed: detail.append(("noise stretch KEPT", t["genome"][:18], t["family"], t["side"], f"stretch {stretch_bp} bp"))
    print(f"   floor {f}: {c}")
    for d in detail: print("      ", *d)

print("\n5. NAMED CASES at floor 39")
for g, fam, contig in (("GCA_982397435.1_T48-F", "HD", "CEVXIV010000004.1"), ("GCA_002105055.1_Leucr1", "bLocus", "MCGR01000028.1")):
    for f in FLOORS:
        x = find(f, g[:28], contig, fam)
        if x: print(f"   {g[:22]} {fam} floor {f}: cluster {x['start']}-{x['end']}, core {x['core_span']['start']}-{x['core_span']['end']}, supported {x['supported_span']['start']}-{x['supported_span']['end']} (beyond supported {x['supported_span']['beyond_supported_bp']})")
