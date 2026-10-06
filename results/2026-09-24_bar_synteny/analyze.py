"""Do withheld loci sit where MAT sits? Adjacency to each genome's SLA2/APN2 orthologs.

Ortholog = the top tblastn hit (bitscore) for that gene anywhere in the genome,
independent of the pipeline's own hits. A locus is ADJACENT when it is on the
same contig as either ortholog and the gap between the two intervals is
<= GAP bp. Populations, per pilot:
  CALLED    reported loci (validates the test: should be mostly adjacent)
  BG        withheld loci in called genomes (specificity: should be rare)
  WH-best   the best withheld locus in genomes with no call
  WH-any    any withheld locus in genomes with no call is adjacent
Also: what a bar of 1 or 0 would report in no-call genomes, split by adjacency.
"""
import collections, csv, glob, os, sys
import yaml
from yaml import CSafeLoader as L

GAP = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
D = os.path.dirname(os.path.abspath(__file__))
P = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-24_pilots"
meta = {r["ASMID"]: r for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"))}

def orthologs(g):
    p = f"{D}/hits/{g}.tsv"
    if not os.path.exists(p):
        return None
    best = {}
    for line in open(p):
        q, s, pid, ln, a, b, ev, bits = line.rstrip("\n").split("\t")
        if q == "none":
            continue
        gene = q.split("|")[0]
        bits = float(bits)
        if gene not in best or bits > best[gene][3]:
            best[gene] = (s, min(int(a), int(b)), max(int(a), int(b)), bits, float(pid))
    return best

import glob as _g
ROLES = {}
for oy in _g.glob("/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/polish-scope-cuts/db/*/order.yml"):
    doc = yaml.safe_load(open(oy))
    for loc in doc.get("loci") or []:
        for gene in loc.get("genes") or []:
            ROLES[(f"{doc['phylum']}:{loc['locus_name']}", gene["name"])] = gene["role"]

def has_core(locus):
    return any(ROLES.get((locus["family"], g)) == "core_MAT" for g in locus.get("genes_found") or [])

def adjacent(locus, orth):
    # A locus counts only if it carries a core MAT gene: SLA2/APN2 are roster
    # genes in some families, so a flank-only cluster is trivially "adjacent".
    if not has_core(locus):
        return False
    for s, a, b, _, _ in orth.values():
        if s == locus["contig"]:
            gap = max(0, max(a, locus["start"]) - min(b, locus["end"]))
            if gap <= GAP:
                return True
    return False

tot = collections.defaultdict(collections.Counter)
by_order = collections.defaultdict(collections.Counter)
for d in sorted(glob.glob(f"{P}/*/runs")):
    name = d.split("/")[-2]
    c = tot[name]
    for g in sorted(os.listdir(d)):
        rp = f"{d}/{g}/detection_report.yaml"
        orth = orthologs(g)
        if orth is None or not os.path.exists(rp):
            c["no_flank_data"] += 1; continue
        if not orth:
            c["no_ortholog_found"] += 1
        doc = yaml.load(open(rp), Loader=L) or {}
        det, sup = doc.get("detected") or [], doc.get("suppressed_loci") or []
        if det:
            c["called_genomes"] += 1
            c["called_adj"] += any(adjacent(x, orth) for x in det)
            for x in sup:
                c["bg_n"] += 1; c["bg_adj"] += adjacent(x, orth)
        elif sup:
            c["wh_genomes"] += 1
            best = max(sup, key=lambda x: (len(x["genes_found"]), x["polished_genes"]))
            c["whbest_adj"] += adjacent(best, orth)
            adj = [x for x in sup if adjacent(x, orth)]
            c["whany_adj"] += bool(adj)
            # What a lowered bar would report in this genome
            for bar in (1, 0):
                rep = [x for x in sup if x["polished_genes"] >= bar]
                if rep:
                    c[f"bar{bar}_genomes"] += 1
                    c[f"bar{bar}_loci"] += len(rep)
                    c[f"bar{bar}_adj_loci"] += sum(adjacent(x, orth) for x in rep)
                    c[f"bar{bar}_genome_has_adj"] += any(adjacent(x, orth) for x in rep)
            o = meta.get(g, {}).get("ORDER") or "?"
            by_order[(name, o)]["wh"] += 1
            by_order[(name, o)]["adj"] += bool(adj)

def pct(a, b):
    return f"{a}/{b} ({100*a/b:.0f}%)" if b else "0/0"

print(f"adjacency window: {GAP} bp")
for name, c in tot.items():
    print(f"\n## {name}  (no ortholog found in {c['no_ortholog_found']} genomes)")
    print(f"  CALLED  genomes with a called locus adjacent: {pct(c['called_adj'], c['called_genomes'])}")
    print(f"  BG      withheld loci adjacent in called genomes: {pct(c['bg_adj'], c['bg_n'])}")
    print(f"  WH-best best withheld locus adjacent: {pct(c['whbest_adj'], c['wh_genomes'])}")
    print(f"  WH-any  any withheld locus adjacent:  {pct(c['whany_adj'], c['wh_genomes'])}")
    for bar in (1, 0):
        print(f"  bar={bar}: would call {c[f'bar{bar}_genomes']} no-call genomes with {c[f'bar{bar}_loci']} loci; "
              f"adjacent loci {pct(c[f'bar{bar}_adj_loci'], c[f'bar{bar}_loci'])}; "
              f"genomes with an adjacent one {pct(c[f'bar{bar}_genome_has_adj'], c[f'bar{bar}_genomes'])}")
print("\nby order (withheld-only genomes with an adjacent withheld locus):")
for (name, o), v in sorted(by_order.items(), key=lambda x: -x[1]["wh"]):
    if v["wh"] >= 3:
        print(f"  {name:26s} {o:22s} {pct(v['adj'], v['wh'])}")
