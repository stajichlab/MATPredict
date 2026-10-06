"""Flank ablation summary: which roster change produces which calls, at what cost."""
import collections, csv, glob, os, statistics as st
import yaml
from yaml import CSafeLoader as L
E = os.path.dirname(os.path.abspath(__file__))
sp = {l.split("\t")[0]: l.rstrip("\n").split("\t")[2] for l in open(f"{E}/genomes.tsv")}
grp = lambda s: "C. albicans" if s == "Candida albicans" else "C. auris" if s == "Candidozyma auris" else "other Serinales"
V = ["V0_base", "V1_alpha2only", "V2_flanksonly", "V3_both", "V3m_both_miniprotflanks", "V4_flanksminiprot_noalpha2", "V5_flanksminiprot_alpha2"]
CORE = {"MTLA1", "MTLA2", "MTLalpha1", "MTLalpha2"}; FL = {"PAP1", "OBP1", "PIK1"}
def modelled(e): return (e.get("status") or "").startswith("polished") or e.get("method") == "diamond_proteome"
def geno(det):
    s = {("homothallic" if x["locus_class"] == "homothallic_candidate" else x["idiomorph"]) for x in det}
    return "+".join(sorted(s)) or "none"
res = {}
for v in V:
    for rp in glob.glob(f"{E}/runs/{v}/*/detection_report.yaml"):
        g = rp.split("/")[-2]; d = yaml.load(open(rp), Loader=L) or {}
        w = open(os.path.dirname(rp) + "/wall_seconds").read().strip()
        res[(v, g)] = (d.get("detected") or [], int(w))
gs = sorted(g for g in sp if all((v, g) in res for v in V))
print(f"{len(gs)} genomes finished in all {len(V)} variants\n")
print(f"{'variant':26s} {'group':16s} {'called':>7s} {'loci/g':>7s} {'median s':>9s}  genotypes")
for v in V:
    for G in ("C. albicans", "C. auris", "other Serinales"):
        sub = [g for g in gs if grp(sp[g]) == G]
        if not sub: continue
        det = [res[(v, g)][0] for g in sub]; w = [res[(v, g)][1] for g in sub]
        gc = collections.Counter(geno(d) for d in det)
        print(f"{v:26s} {G:16s} {sum(bool(d) for d in det):3d}/{len(sub):<3d} {sum(map(len, det))/len(sub):7.2f} {st.median(w):9.0f}  {dict(sorted(gc.items()))}")
    print()
print("flank-carried calls (no modelled core gene) per variant:")
for v in V:
    n = sum(1 for g in gs for x in res[(v, g)][0]
            if not any(e["gene"] in CORE and modelled(e) for e in x["gene_evidence"]))
    print(f"  {v:26s} {n}")
print("\nagreement with V3_both (same genotype per genome):")
for v in V:
    same = sum(geno(res[(v, g)][0]) == geno(res[("V3_both", g)][0]) for g in gs)
    print(f"  {v:26s} {same}/{len(gs)}")
for x, y in (("V3_both", "V3m_both_miniprotflanks"), ("V3m_both_miniprotflanks", "V5_flanksminiprot_alpha2"), ("V5_flanksminiprot_alpha2", "V4_flanksminiprot_noalpha2")):
    diff = [g for g in gs if geno(res[(x, g)][0]) != geno(res[(y, g)][0])]
    print(f"\n{y} vs {x}: {len(gs) - len(diff)}/{len(gs)} same genotype; differ: {[(g, sp[g], geno(res[(x, g)][0]), geno(res[(y, g)][0])) for g in diff]}")
alpha2 = sum(1 for g in gs for x in res[("V5_flanksminiprot_alpha2", g)][0] if "MTLalpha2" in (x.get("genes_found") or []))
print(f"V5 calls that include MTLalpha2: {alpha2}")
