"""Flank-ortholog synteny at sexM/sexP loci, early-diverging scan.

Evidence per genome (from run_synteny.slurm):
  hits/<g>.tbn.zst      forward tblastn, curated Mucoromycota proteins, e<=10
  hits/<g>.back.zst     reverse blastp of each flank region vs 3 Mucorales proteomes

Definitions
  flank region     HSPs of one flank gene on one contig within 3 kb (regions.py)
  confirmed        the region's top reverse hit (bitscore) is that gene's ortholog
                   in a Mucorales proteome (ref/flank_orthologs.tsv)
  genome ortholog  the confirmed region with the best forward bitscore; RBH-strict
                   additionally requires it to be the gene's best forward region
  HMG hit          a sexM or sexP HSP with e <= HMG_E
Locus tests (loci from detection reports carrying sexM or sexP):
  flank_conf       >=1 flank region within the locus span (+-PAD) is confirmed
  hmg_between      an HMG hit lies between two confirmed regions of DIFFERENT
                   flank genes on the same contig, no more than SPAN apart
Genome test (independent of the pipeline):
  pair             the genome orthologs of tptA and rnhA share a contig within SPAN
  pair_hmg         an HMG hit lies between them
  near_hmg         an HMG hit within NEAR bp of any genome ortholog
"""
import collections, csv, glob, io, os, subprocess, sys
import yaml
from yaml import CSafeLoader as L

D = os.path.dirname(os.path.abspath(__file__))
ED = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-26_early_diverging"
import os as _o
PAD, SPAN, NEAR, HMG_E = 5000, int(_o.environ.get("SPAN",100000)), int(_o.environ.get("NEAR",20000)), float(_o.environ.get("HMG_E",1e-3))
meta = {r["ASMID"]: r for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"))}
ORTH = {}
for line in open(f"{D}/ref/flank_orthologs.tsv"):
    g, ids = line.rstrip("\n").split("\t")
    ORTH[g] = set(ids.split(","))


def zread(p):
    if not os.path.exists(p):
        return []
    return subprocess.run(["zstd", "-dc", p], capture_output=True, text=True).stdout.splitlines()


def load(g):
    regions = {}
    for q in zread(f"{D}/hits/{g}.regions.faa.zst"):
        if q.startswith(">"):
            _, gene, contig, s, e, bits = q[1:].split("|")
            regions[q[1:]] = dict(gene=gene, contig=contig, s=int(s), e=int(e), bits=float(bits), conf=False, top=None)
    best = {}
    for line in zread(f"{D}/hits/{g}.back.zst"):
        q, s, pid, ln, ev, bits = line.split("\t")
        if q not in best or float(bits) > best[q][1]:
            best[q] = (s, float(bits))
    for k, r in regions.items():
        if k in best:
            r["top"] = best[k][0]
            r["conf"] = best[k][0] in ORTH[r["gene"]]
    hmg = []
    for line in zread(f"{D}/hits/{g}.tbn.zst"):
        f = line.split("\t")
        gene = f[0].split("|")[2]
        if gene in ("sexM", "sexP") and float(f[8]) <= HMG_E:
            a, b = sorted((int(f[6]), int(f[7])))
            hmg.append((f[1], a, b))
    return list(regions.values()), hmg


def genome_orthologs(regions):
    out = {}
    for gene in ORTH:
        rs = [r for r in regions if r["gene"] == gene]
        if not rs:
            continue
        top = max(rs, key=lambda r: r["bits"])
        conf = [r for r in rs if r["conf"]]
        if conf:
            o = max(conf, key=lambda r: r["bits"])
            out[gene] = dict(o, strict=(o is top))
    return out


def between(hmg, contig, a, b):
    lo, hi = min(a["s"], b["s"]), max(a["e"], b["e"])
    return any(c == contig and s >= lo and e <= hi for c, s, e in hmg)


def near(hmg, r, dist):
    return any(c == r["contig"] and max(0, max(s, r["s"]) - min(e, r["e"])) <= dist for c, s, e in hmg)


rows = []
for grp in ["Mucoromycota", "Mortierellomycota", "Kickxellomycota", "chytrid_control"]:
    for rp in sorted(glob.glob(f"{ED}/{grp}/runs/*/detection_report.yaml")):
        g = rp.split("/")[-2]
        if not os.path.exists(f"{D}/hits/{g}.done"):
            continue
        regions, hmg = load(g)
        orth = genome_orthologs(regions)
        m = meta.get(g, {})
        pair = pair_hmg = False
        if "tptA" in orth and "rnhA" in orth:
            t, r = orth["tptA"], orth["rnhA"]
            if t["contig"] == r["contig"] and max(0, max(t["s"], r["s"]) - min(t["e"], r["e"])) <= SPAN:
                pair = True
                pair_hmg = between(hmg, t["contig"], t, r)
        base = dict(genome=g, group=grp, family=m.get("FAMILY", "?"), species=m.get("SPECIES", "?"),
                    orth_tptA=("tptA" in orth), orth_rnhA=("rnhA" in orth),
                    strict_tptA=orth.get("tptA", {}).get("strict", False),
                    strict_rnhA=orth.get("rnhA", {}).get("strict", False),
                    pair=pair, pair_hmg=pair_hmg,
                    near_hmg=any(near(hmg, o, NEAR) for o in orth.values()))
        doc = yaml.load(open(rp), Loader=L) or {}
        loci = [("called", x) for x in doc.get("detected") or []] + \
               [("withheld", x) for x in doc.get("suppressed_loci") or []]
        loci = [(k, x) for k, x in loci if {"sexM", "sexP"} & set(x.get("genes_found") or [])]
        base["n_loci"] = len(loci)
        base["called"] = any(k == "called" for k, _ in loci)
        lrows = []
        for kind, x in loci:
            at = [r for r in regions if r["contig"] == x["contig"]
                  and r["e"] >= x["start"] - PAD and r["s"] <= x["end"] + PAD]
            conf = [r for r in at if r["conf"]]
            hb = any(a["gene"] != b["gene"] and max(a["e"], b["e"]) - min(a["s"], b["s"]) <= SPAN
                     and between(hmg, x["contig"], a, b) for a in conf for b in conf)
            is_gorth = any(orth.get(r["gene"]) is not None and orth[r["gene"]]["s"] == r["s"]
                           and orth[r["gene"]]["contig"] == r["contig"] for r in conf)
            lrows.append(dict(kind=kind, flank_at=len(at), flank_conf=bool(conf),
                              conf_genes=",".join(sorted({r["gene"] for r in conf})),
                              is_genome_ortholog=is_gorth, hmg_between=hb))
        base["loci"] = lrows
        rows.append(base)

# ---------------------------------------------------------------- tables
def unit(r):
    return r["group"] if r["group"] != "Mucoromycota" else f"Mucoromycota/{r['family']}"

out = io.StringIO()
p = lambda *a: print(*a, file=out)
p(f"# PAD={PAD} SPAN={SPAN} NEAR={NEAR} HMG_E={HMG_E}\n")
p("## Genome level (independent of the pipeline)")
p(f"{'unit':38s} {'genomes':>7s} {'tptA_orth':>9s} {'rnhA_orth':>9s} {'strict_both':>11s} {'pair':>5s} {'pair+HMG':>8s} {'HMG<=20kb':>9s}")
for grp in ["Mucoromycota", "Mortierellomycota", "Kickxellomycota", "chytrid_control"]:
    units = sorted({unit(r) for r in rows if r["group"] == grp}, key=lambda u: -sum(unit(r) == u for r in rows))
    for u in ([grp] + units if grp == "Mucoromycota" else [grp]):
        rs = [r for r in rows if (unit(r) == u or (u == grp and r["group"] == grp))]
        c = lambda k: sum(bool(r[k]) for r in rs)
        p(f"{u:38s} {len(rs):7d} {c('orth_tptA'):9d} {c('orth_rnhA'):9d} "
          f"{sum(r['strict_tptA'] and r['strict_rnhA'] for r in rs):11d} {c('pair'):5d} {c('pair_hmg'):8d} {c('near_hmg'):9d}")

p("\n## Locus level (sexM/sexP loci from the detection reports)")
p(f"{'unit':38s} {'kind':8s} {'loci':>5s} {'flank@':>6s} {'conf':>5s} {'=gOrth':>6s} {'HMGbtw':>6s}")
for grp in ["Mucoromycota", "Mortierellomycota", "Kickxellomycota", "chytrid_control"]:
    for kind in ["called", "withheld"]:
        ls = [l for r in rows if r["group"] == grp for l in r["loci"] if l["kind"] == kind]
        c = lambda k: sum(bool(l[k]) for l in ls)
        p(f"{grp:38s} {kind:8s} {len(ls):5d} {sum(l['flank_at'] > 0 for l in ls):6d} {c('flank_conf'):5d} {c('is_genome_ortholog'):6d} {c('hmg_between'):6d}")

p("\n## Genomes with >=1 locus where HMG lies between confirmed flanks, per species")
for grp in ["Mortierellomycota", "Kickxellomycota", "Mucoromycota"]:
    sp = collections.defaultdict(lambda: [0, 0])
    for r in rows:
        if r["group"] != grp:
            continue
        sp[r["species"]][0] += 1
        sp[r["species"]][1] += any(l["hmg_between"] for l in r["loci"])
    pos = {k: v for k, v in sp.items() if v[1]}
    p(f"{grp}: species {len(sp)}, species with >=1 such genome {len(pos)}")
    if grp != "Mucoromycota":
        for k, v in sorted(pos.items(), key=lambda kv: -kv[1][1]):
            p(f"    {k:40s} {v[1]}/{v[0]}")

p("\n## Confirmed flank genes at loci, by group")
for grp in ["Mucoromycota", "Mortierellomycota", "Kickxellomycota", "chytrid_control"]:
    cc = collections.Counter(l["conf_genes"] or "-" for r in rows if r["group"] == grp for l in r["loci"])
    p(grp, dict(cc.most_common(8)))

p(f"\ngenomes analysed: {len(rows)}")
open(f"{D}/summary.txt", "w").write(out.getvalue())
print(out.getvalue())
with open(f"{D}/genomes_table.tsv", "w") as fo:
    keys = [k for k in rows[0] if k != "loci"] if rows else []
    fo.write("\t".join(keys + ["loci_hmg_between", "loci_flank_conf"]) + "\n")
    for r in rows:
        fo.write("\t".join(str(r[k]) for k in keys) + f"\t{sum(l['hmg_between'] for l in r['loci'])}\t{sum(l['flank_conf'] for l in r['loci'])}\n")
