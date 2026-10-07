"""Serinales-wide scan: new roster (PAP1/OBP1/PIK1 flanks, optional MTLalpha2) vs baseline.
usage: serinales_scan_analysis.py BASE_DIR NEW_DIR"""
import collections, csv, glob, os, sys, statistics as st
import yaml
from yaml import CSafeLoader as L
base_d, new_d = sys.argv[1], sys.argv[2]
meta = {r["ASMID"]: r for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"))}
FLANKS = {"PAP1", "OBP1", "PIK1"}
CORE = {"MTLA1", "MTLA2", "MTLalpha1", "MTLalpha2"}
def load(d):
    return {rp.split("/")[-2]: (yaml.load(open(rp), Loader=L) or {}) for rp in glob.glob(f"{d}/runs/*/detection_report.yaml")}
def genotype(det):
    s = set()
    for x in det:
        if x["locus_class"] == "homothallic_candidate": s.add("homothallic")
        elif x["idiomorph"] == "A": s.add("a")
        elif x["idiomorph"] == "alpha": s.add("alpha")
        else: s.add("undetermined")
    return "+".join(sorted(s)) or "none"
def modelled(e): return (e.get("status") or "").startswith("polished") or e.get("method") == "diamond_proteome"
base, new = load(base_d), load(new_d)
print(f"genomes: base {len(base)}  new {len(new)}  both {len(set(base)&set(new))}")
by = collections.defaultdict(lambda: collections.defaultdict(collections.Counter))
extra = collections.defaultdict(collections.Counter)
for g in sorted(set(base) & set(new)):
    sp = meta.get(g, {}).get("SPECIES") or "?"
    for label, doc in (("base", base[g]), ("new", new[g])):
        det = doc.get("detected") or []
        by[sp][label][genotype(det)] += 1
        by[sp][label]["_loci"] += len(det)
        if label == "new":
            for x in det:
                genes = set(x.get("genes_found") or [])
                ev = {e["gene"]: e for e in x.get("gene_evidence") or []}
                mcore = [n for n in CORE if n in ev and modelled(ev[n])]
                mfl = [n for n in FLANKS if n in ev and modelled(ev[n])]
                c = extra[sp]
                c["loci"] += 1
                c["with_flank"] += bool(genes & FLANKS)
                c["with_2plus_flanks"] += len(genes & FLANKS) >= 2
                if x["idiomorph"] == "alpha":
                    c["alpha_loci"] += 1; c["alpha_with_alpha2"] += "MTLalpha2" in genes
                if not mcore and mfl:
                    c["flank_carried"] += 1   # passed the bar on flank models alone
                    c["flank_carried_single_flank"] += len(genes & FLANKS) == 1
    ob, on = genotype(base[g].get("detected") or []), genotype(new[g].get("detected") or [])
    if ob != on: extra[sp][f"change {ob} -> {on}"] += 1
print(f"\n{'species':32s} {'n':>5}  genotype base -> new (counts)")
tot = collections.Counter()
for sp, v in sorted(by.items(), key=lambda kv: -sum(kv[1]['new'][k] for k in kv[1]['new'] if not k.startswith('_'))):
    n = sum(c for k, c in v["new"].items() if not k.startswith("_"))
    if n < 5: 
        for lab in ("base","new"): tot.update({f"{lab}:{k}": c for k, c in v[lab].items()})
        continue
    b = {k: c for k, c in v["base"].items() if not k.startswith("_")}; nn = {k: c for k, c in v["new"].items() if not k.startswith("_")}
    print(f"{sp[:32]:32s} {n:5d}  base {dict(sorted(b.items()))}")
    print(f"{'':32s} {'':5s}  new  {dict(sorted(nn.items()))}  loci/genome {v['base']['_loci']/n:.2f} -> {v['new']['_loci']/n:.2f}")
    e = extra[sp]
    print(f"{'':32s} {'':5s}  new loci {e['loci']}: with a flank {e['with_flank']}, >=2 flanks {e['with_2plus_flanks']}, "
          f"alpha loci {e['alpha_loci']} (with alpha2 {e['alpha_with_alpha2']}), flank-carried {e['flank_carried']} (single-flank {e['flank_carried_single_flank']})")
    ch = {k: c for k, c in e.items() if k.startswith("change")}
    if ch: print(f"{'':32s} {'':5s}  changes: {dict(sorted(ch.items(), key=lambda x:-x[1])[:8])}")
print("\nspecies with <5 genomes, pooled:", {k: c for k, c in sorted(tot.items()) if not k.split(':')[1].startswith('_')})
