import csv, subprocess, sys, os, collections
sys.path.insert(0, "/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/flank-first-scan/scripts")
from homeodomain_tree import Prot, table
B = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-11-13_flank_first_scan"
D = B + "/stage1_all/per_genome"
P = Prot()
oeh = "Oehlia_diaphana_SD1"
# partner genomes: pairs with >=0.96 HD2 identity to Oehlia
partners = collections.Counter()
for cls in ("HD1", "HD2"):
    for r in csv.DictReader(open(f"{B}/hd_alleles/hd_pairs_AMF_{cls}.tsv"), delimiter="\t"):
        a, b = r["a"].split("|")[0], r["b"].split("|")[0]
        if float(r["aligned_identity"]) >= 0.96 and oeh in (a, b):
            partners[(cls, b if a == oeh else a)] += 1
print("partners at >=96% to Oehlia:", dict(partners))
def flank(g):
    cand = {c[0]: c for c in table(f"{D}/{g}.candidates.tsv")}
    out = {}
    for c in cand.values():
        if c[5] in ("HD1", "HD2", "choline", "PGAM"):
            out.setdefault(c[5], []).append(c[0])
    return out
def aln_id(s1, s2):
    import tempfile
    t = os.environ.get("SCRATCH", "/tmp")
    f = os.path.join(t, "p.faa")
    open(f, "w").write(f">a\n{s1}\n>b\n{s2}\n")
    o = subprocess.run(["mafft", "--quiet", "--auto", f], capture_output=True, text=True).stdout.split(">")[1:]
    A = ["".join(x.split("\n")[1:]).upper() for x in o]
    both = same = 0
    for x, y in zip(*A):
        if x != "-" and y != "-":
            both += 1; same += x == y
    return same / both if both else 0
F = flank(oeh)
print("Oehlia genes:", F)
for (cls, pg) in sorted(partners):
    G = flank(pg)
    print("== partner", pg, G)
    for k in ("HD1", "HD2", "choline", "PGAM"):
        best = 0
        for x in F.get(k, []):
            for y in G.get(k, []):
                best = max(best, aln_id(P.get(oeh)[x], P.get(pg)[y]))
        print(f"   {k}: best aligned identity {best:.3f}")
