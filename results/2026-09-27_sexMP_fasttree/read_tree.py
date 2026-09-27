"""Step 7: assign each candidate to the sexM clade, the sexP clade, or neither.

Unrooted-safe: for each type T, find the smallest tree side (bipartition) that holds
every curated reference of type T and no reference of the other type and no MAT1-2-1.
Its support (SH-aLRT/UFBoot) is the node label IQ-TREE writes. If no such side exists
(type not monophyletic), the pure side holding the most type-T references is used,
and the clades file says so. A candidate is placed
in T if it is on that side. Orthology is claimed only when UFBoot >= 95.
Writes candidates.tsv (one row per locus x extracted protein) and clade_summary.tsv.
Usage: read_tree.py PREFIX
"""
import collections, csv, sys
from Bio import Phylo
prefix = sys.argv[1] if len(sys.argv) > 1 else "tree_ft"
ids = dict(l.rstrip("\n").split("\t") for l in list(open("taxon_ids.tsv"))[1:])
tree = Phylo.read(f"{prefix}.treefile", "newick")
leaves = {t.name for t in tree.get_terminals()}
name = {t: ids[t] for t in leaves}
kind = {}
for t, n in name.items():
    if n.startswith("OUT|"): kind[t] = "MAT1-2-1"
    elif n.startswith("REF|"): kind[t] = n.split("|")[-1]
    else: kind[t] = "query"
refs = {k: {t for t in leaves if kind[t] == k} for k in ("sexM", "sexP", "MAT1-2-1")}
def label(cl):
    # FastTree writes one local support value (Shimodaira-Hasegawa-like, 0-1) per node;
    # it is reported in both slots, scaled to 0-100. It is NOT UFBoot.
    c = cl.confidence
    if c is None and cl.name:
        try: c = float(cl.name)
        except ValueError: c = None
    if c is None: return None, None
    c = c * 100 if c <= 1 else c
    return c, c
sides = []
for cl in tree.find_clades(order="level"):
    if cl.is_terminal() or cl == tree.root: continue
    S = {t.name for t in cl.get_terminals()}
    a, b = label(cl)
    sides.append((S, a, b)); sides.append((leaves - S, a, b))
result = {}
for T, other in (("sexM", "sexP"), ("sexP", "sexM")):
    X = refs[other] | refs["MAT1-2-1"]
    ok = [s for s in sides if refs[T] <= s[0] and not (s[0] & X)]
    if ok:
        S, a, b = min(ok, key=lambda s: len(s[0]))
        result[T] = (S, a, b, len(refs[T]))
    else:
        # not monophyletic: take the pure side (no other-type ref, no MAT1-2-1)
        # that holds the most refs of type T, smallest if tied
        pure = [s for s in sides if not (s[0] & X)]
        S, a, b = max(pure, key=lambda s: (len(s[0] & refs[T]), -len(s[0])))
        result[T] = (S, a, b, len(S & refs[T]))
with open(f"{prefix}.clades.txt", "w") as fo:
    for T, (S, a, b, k) in result.items():
        mono = "monophyletic" if k == len(refs[T]) else f"NOT monophyletic: best pure clade holds {k} of {len(refs[T])} refs"
        out = sorted(name[t] for t in refs[T] - S)
        fo.write(f"{T}: {mono}; clade size {len(S)}; FastTree local support {a}; refs outside: {out}\n")
print(open(f"{prefix}.clades.txt").read())
# map representative -> clade, then expand dedup members
rep_clade = {}
for t in leaves:
    c = "sexM" if t in result["sexM"][0] else "sexP" if t in result["sexP"][0] else "other_HMG"
    rep_clade[name[t]] = c
members = collections.defaultdict(list)
for r in csv.DictReader(open("dedup_map.tsv"), delimiter="\t"):
    members[r["rep"]].append(r["member"])
clade_of = {m: rep_clade[rep] for rep, ms in members.items() for m in ms if rep in rep_clade}
loci = {r["locus_id"]: r for r in csv.DictReader(open("loci.tsv"), delimiter="\t")}
clf = {}
for r in csv.DictReader(open("/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-26_hmm_classifier/classifier_calls.tsv"), delimiter="\t"):
    clf[(r["genome"], r["contig"], r["start"])] = (r["clf_verdict"], r["clf_margin"])
nodom = set(l.strip() for l in open("no_hmg_domain.txt") if l.strip())
rows = []
for l in open("candidates_raw.tsv"):
    f = l.rstrip("\n").split("\t")
    cid, g, knd, wid, method, q, score, cov, plen = f
    L = loci.get(wid, {})
    tree_clade = clade_of.get(cid, "no_HMG_box" if cid in nodom or method == "none" else "not_in_tree")
    lab = L.get("idiomorph", "")
    agree = ""
    if tree_clade in ("sexM", "sexP") and lab in ("Plus", "Minus"):
        agree = str((tree_clade == "sexM") == (lab == "Minus"))
    rows.append(dict(candidate_id=cid, genome=g, kind=knd, locus_or_copy=wid, method=method,
                     best_ref=q, score=score, ref_coverage_or_E=cov, protein_len=plen,
                     group=L.get("group", ""), family=L.get("family", ""), species=L.get("species", ""),
                     status=L.get("status", "nonlocus"), idiomorph_label=lab,
                     contig=L.get("contig", ""), start=L.get("start", ""), end=L.get("end", ""),
                     tree_clade=tree_clade, label_agrees=agree,
                     clf_verdict=clf.get((g, L.get("contig", ""), L.get("start", "")), ("", ""))[0],
                     clf_margin=clf.get((g, L.get("contig", ""), L.get("start", "")), ("", ""))[1]))
with open(f"candidates.tsv", "w", newline="") as fo:
    w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
# per locus: best clade present among its proteins
per = collections.defaultdict(set)
for r in rows:
    if r["kind"] == "cand": per[r["locus_or_copy"]].add(r["tree_clade"])
summ = collections.Counter()
for lid, cs in per.items():
    L = loci[lid]
    grp = L["group"] if L["group"] != "Mucoromycota" else f"Mucoromycota/{L['family']}"
    c = "sexM+sexP" if {"sexM", "sexP"} <= cs else "sexM" if "sexM" in cs else "sexP" if "sexP" in cs else \
        "other_HMG" if "other_HMG" in cs else "no_HMG_box"
    summ[(grp, L["status"], c)] += 1
    if c in ("sexM", "sexP") and L["idiomorph"] in ("Plus", "Minus"):
        summ[(grp, L["status"], "label_agree" if (c == "sexM") == (L["idiomorph"] == "Minus") else "label_disagree")] += 1
nonloc = collections.Counter(r["tree_clade"] for r in rows if r["kind"] == "nonlocus")
with open("clade_summary.tsv", "w") as fo:
    fo.write("group\tstatus\tclass\tloci\n")
    for k, v in sorted(summ.items()): fo.write("\t".join(k) + f"\t{v}\n")
    for k, v in sorted(nonloc.items()): fo.write(f"nonlocus_outgroup_copies\t-\t{k}\t{v}\n")
print(open("clade_summary.tsv").read())
