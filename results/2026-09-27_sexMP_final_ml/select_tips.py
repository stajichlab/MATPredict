"""Select tips for the trimmed Mucoromycotina sexM/sexP final ML tree.

Source: the FastTree run (results/2026-09-27_sexMP_fasttree), whose loci match the
4457c3a scan (0 new / 0 lost calls vs 1a00b0a); labels updated to 4457c3a.
Scope: Mucoromycota orders Mucorales + Umbelopsidales only (curator ruling
2026-09-27: Mortierellomycota, Kickxellomycota, Endogonales etc. are discovery-only).

Included:
  * the 15 curated sexM/sexP references and the 2 Ascomycota MAT1-2-1
  * every tip standing for a CALLED Mucoromycotina locus
  * withheld / non-locus tips that fell in the sexP clade or sexM group
  * other-HMG outgroups: greedy identity clustering of the aligned 69-column
    HMG boxes of the remaining in-scope other-HMG tips (threshold OUT_ID),
    one representative per cluster (the longest ungapped), capped at OUT_CAP
    by keeping the largest clusters.
"""
import collections, csv

FT = "../2026-09-27_sexMP_fasttree"
OUT_ID, OUT_CAP = 0.55, 70
IN_ORDERS = {"Mucorales", "Umbelopsidales"}

def read_fa(f):
    d, n = {}, None
    for l in open(f):
        l = l.rstrip()
        if l.startswith(">"):
            n = l[1:].split()[0]; d[n] = []
        elif n:
            d[n].append(l)
    return {k: "".join(v) for k, v in d.items()}

samp = {r["ASMID"]: r for r in csv.DictReader(open(
    "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"))}
tips = list(csv.DictReader(open(f"{FT}/tip_names.tsv"), delimiter="\t"))
cand = {r["candidate_id"]: r for r in csv.DictReader(open(f"{FT}/candidates.tsv"), delimiter="\t")}
members = collections.defaultdict(list)
for r in csv.DictReader(open(f"{FT}/dedup_map.tsv"), delimiter="\t"):
    members[r["rep"]].append(r["member"])
aln = read_fa(f"{FT}/aln_hmm.afa")
undet = {}
for r in csv.DictReader(open("../2026-09-27_btbA_homothallic/compare_vs_1a00b0a.tsv"), delimiter="\t"):
    if r["kind"].startswith("changed:idiomorph"):
        undet[(r["genome"], r["contig"])] = r["new_idiomorph"]

def genome_of(cid):
    return cid.split("|")[0]

def in_scope(cid):
    s = samp.get(genome_of(cid), {})
    return s.get("PHYLUM") == "Mucoromycota" and s.get("ORDER") in IN_ORDERS

def label_of(cid):
    c = cand.get(cid, {})
    lab = c.get("idiomorph_label", "")
    return undet.get((c.get("genome"), c.get("contig")), lab)

rows, other = [], []
for t in tips:
    src = t["source"]
    if t["status"] == "reference":
        rows.append(dict(seqid=src, why="reference", kind=t["kind"], genome="", species=t["species"],
                         family="", status="reference", label="", ft_clade=t["tree_clade"], n_members=1))
        continue
    mem = members.get(src, [src])
    scoped = [m for m in mem if in_scope(m)]
    if not scoped:
        continue
    called = [m for m in scoped if cand.get(m, {}).get("status") == "called"]
    pick = called[0] if called else scoped[0]
    c = cand.get(pick, {})
    s = samp.get(genome_of(pick), {})
    base = dict(seqid=src, kind=t["kind"], genome=genome_of(pick), species=s.get("SPECIES", ""),
                family=s.get("FAMILY", ""), status=c.get("status", t["status"]), label=label_of(pick),
                ft_clade=t["tree_clade"], n_members=len(scoped), rep_member=pick)
    if called:
        rows.append(dict(base, why="called_locus"))
    elif t["tree_clade"] in ("sexP", "sexM"):
        rows.append(dict(base, why="withheld_or_nonlocus_in_clade"))
    else:
        other.append(base)

# outgroup subsample: greedy clustering on aligned HMG box identity
def ident(a, b):
    n = m = 0
    for x, y in zip(a, b):
        if x == "-" or y == "-":
            continue
        n += 1; m += x == y
    return m / n if n >= 30 else 0.0

other.sort(key=lambda r: -sum(ch != "-" for ch in aln[r["seqid"]]))
clusters = []
for r in other:
    s = aln[r["seqid"]]
    for cl in clusters:
        if ident(s, aln[cl[0]["seqid"]]) >= OUT_ID:
            cl.append(r); break
    else:
        clusters.append([r])
clusters.sort(key=lambda cl: -len(cl))
for cl in clusters[:OUT_CAP]:
    rows.append(dict(cl[0], why=f"outgroup_rep(cluster_n={len(cl)})"))

with open("selection.tsv", "w", newline="") as fo:
    keys = ["seqid", "why", "kind", "genome", "species", "family", "status", "label", "ft_clade",
            "n_members", "rep_member"]
    w = csv.DictWriter(fo, fieldnames=keys, delimiter="\t", extrasaction="ignore")
    w.writeheader(); w.writerows(rows)
with open("selected.afa", "w") as fo:
    for r in rows:
        fo.write(f">{r['seqid']}\n{aln[r['seqid']]}\n")
comp = collections.Counter((r["why"].split("(")[0], r.get("ft_clade", "")) for r in rows)
print(f"other-HMG in scope: {len(other)}; clusters at {OUT_ID}: {len(clusters)}; kept {min(len(clusters), OUT_CAP)}")
for k, v in sorted(comp.items()):
    print(k, v)
print("total", len(rows))
