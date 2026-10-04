#!/usr/bin/env python3
"""Build the sexM/sexP tree input set.

Ingroup: called core proteins (../core_proteins.faa; Plus -> sexP, Minus -> sexM,
polished models only) plus the curated record core proteins (db/Mucoromycota/*/*/
proteins.faa, role=core_MAT). Genomes whose override `use` says "hold out" are
dropped. Confirmed overrides are named by their likely identity.
Proteins shorter than MIN_AA are dropped. Redundancy is removed later by cd-hit
within each idiomorph (run_tree.sh).

Outgroup (curator ruling 2026-10-03, non-MAT HMG): the classifier's HMG-paralog
negative set (paralog_negatives.faa: non-locus HMG copies outside the sexM and sexP
clades) plus the P1 sexM-like paralog.

Writes ingroup_sexP.faa, ingroup_sexM.faa, outgroup.faa, tips.tsv (tip id -> label).
"""
import csv, glob, re, sys
C = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-03_mucoromycotina_mat"
T = f"{C}/tree"
WT = "/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-7c7ed99"
MIN_AA = 60

def fasta(p):
    s, n = {}, None
    for l in open(p):
        l = l.rstrip()
        if l.startswith(">"):
            n = l[1:]; s[n] = []
        elif n:
            s[n].append(l)
    return {k: "".join(v).replace("*", "") for k, v in s.items()}

def safe(x):
    return re.sub(r"[^A-Za-z0-9.-]+", "_", x).strip("_")

ov = {}
for l in open(f"{WT}/db/taxon_overrides.tsv"):
    f = l.rstrip("\n").split("\t")
    if len(f) >= 9 and not l.startswith("#") and f[0] != "genome_id":
        ov[f[0]] = f
calls = {r["locus_id"]: r for r in csv.DictReader(open(f"{C}/calls.tsv"), delimiter="\t") if r.get("locus_id")}
tips, seqs = [], {"Plus": {}, "Minus": {}}
dropped = {"hold_out": 0, "short": 0}
for h, s in fasta(f"{C}/core_proteins.faa").items():
    loc, idio, gene = h.split("|")
    r = calls[loc]
    o = ov.get(r["genome"])
    if o and "hold out" in o[6]:
        dropped["hold_out"] += 1; continue
    if len(s) < MIN_AA:
        dropped["short"] += 1; continue
    name = r["name"]
    if o and o[5] == "confirmed":
        name = f"{o[2]} [{r['name_in_source']}]"
    tid = f"{r['source']}_{safe(r['genome'])}_{loc.rsplit('__', 1)[1]}"
    seqs[idio][tid] = s
    tips.append(dict(tip=tid, label=f"{name} ({r['source']})", idiomorph=idio, gene=gene,
                     source=r["source"], genome=r["genome"], name=name, confidence=r["confidence"],
                     locus_class=r["locus_class"], aa=len(s), override=(o[5] if o else ""), kind="call"))
for p in sorted(glob.glob(f"{WT}/db/Mucoromycota/*/*_MAT_*/proteins.faa")):
    order = p.split("/")[-3]
    for h, s in fasta(p).items():
        if "role=core_MAT" not in h:
            continue
        rec = h.split("|")[0]
        gene = re.search(r"name=([^|]+)", h).group(1)
        idio = {"sexP": "Plus", "sexM": "Minus"}.get(gene)
        if not idio:
            continue
        tid = f"REC_{safe(rec)}"
        seqs[idio][tid] = s
        tips.append(dict(tip=tid, label=f"{rec} (record, {order})", idiomorph=idio, gene=gene, source="record",
                         genome=rec, name=rec, confidence="curated", locus_class="record", aa=len(s),
                         override="", kind="record"))
out = []
neg = fasta(f"{WT}/db/Mucoromycota/classifiers/MAT/paralog_negatives.faa")
for h, s in neg.items():
    tid = "OUT_" + safe(h)
    out.append((tid, s))
    tips.append(dict(tip=tid, label=f"{h} (non-MAT HMG)", idiomorph="outgroup", gene="HMG", source="BFD",
                     genome=h.split("|")[0], name="", confidence="", locus_class="", aa=len(s), override="", kind="outgroup"))
for h, s in fasta(f"{WT}/db/Mucoromycota/classifiers/MAT/paralogs/P1.faa").items():
    tid = "OUT_P1_" + safe(h)
    out.append((tid, s))
    tips.append(dict(tip=tid, label=f"P1 sexM-like paralog {h}", idiomorph="outgroup", gene="P1", source="BFD",
                     genome=h, name="", confidence="", locus_class="", aa=len(s), override="", kind="outgroup"))
for idio, g in (("Plus", "sexP"), ("Minus", "sexM")):
    with open(f"{T}/ingroup_{g}.faa", "w") as fo:
        for k, v in seqs[idio].items():
            fo.write(f">{k}\n{v}\n")
with open(f"{T}/outgroup.faa", "w") as fo:
    for k, v in out:
        fo.write(f">{k}\n{v}\n")
with open(f"{T}/tips.tsv", "w") as fo:
    w = csv.DictWriter(fo, fieldnames=list(tips[0]), delimiter="\t", lineterminator="\n")
    w.writeheader(); w.writerows(tips)
print({k: len(v) for k, v in seqs.items()}, "outgroup", len(out), "dropped", dropped)
