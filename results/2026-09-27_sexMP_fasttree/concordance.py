"""Step 6: label vs clade concordance, before (earlier labels, IQ-TREE clades) and
after (classifier-scan labels, FastTree clades). Unit = locus; a locus takes the
clade(s) of its extracted proteins (sexM, sexP, both, or neither). Called loci only
for label counts. Writes concordance.txt and disagreements.tsv."""
import collections, csv
def load(path):
    per = collections.defaultdict(set); info = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r["kind"] != "cand": continue
        per[r["locus_or_copy"]].add(r["tree_clade"]); info[r["locus_or_copy"]] = r
    return per, info
def clade(cs):
    return "both" if {"sexM", "sexP"} <= cs else "sexM" if "sexM" in cs else "sexP" if "sexP" in cs else "none"
out = []
for tag, path in (("before", "../2026-09-26_sexMP_phylogeny/candidates.tsv"), ("after", "candidates.tsv")):
    per, info = load(path)
    tab = collections.Counter(); fam = collections.defaultdict(collections.Counter)
    for lid, cs in per.items():
        r = info[lid]
        if r["status"] != "called": continue
        c = clade(cs); lab = r["idiomorph_label"]
        tab[(c, lab)] += 1
        grp = r["group"] if r["group"] != "Mucoromycota" else r["family"]
        if c in ("sexP", "sexM"):
            ok = (c == "sexP" and lab == "Plus") or (c == "sexM" and lab == "Minus")
            fam[grp]["agree" if ok else "disagree"] += 1
    out.append(f"== {tag} ({path})")
    for c in ("sexP", "sexM", "both"):
        out.append(f"  {c:5s} clade: " + ", ".join(f"{lab or '-'}={n}" for (cc, lab), n in sorted(tab.items()) if cc == c))
    tot = sum(v["agree"] for v in fam.values()), sum(v["disagree"] for v in fam.values())
    out.append(f"  label agrees with clade: {tot[0]} / {tot[0] + tot[1]} called loci placed in sexP or sexM")
    for g, v in sorted(fam.items()):
        out.append(f"    {g:28s} agree {v['agree']:3d}  disagree {v['disagree']:3d}")
open("concordance.txt", "w").write("\n".join(out) + "\n"); print("\n".join(out))
per, info = load("candidates.tsv")
with open("disagreements.tsv", "w") as fo:
    fo.write("locus\tgenome\tspecies\tfamily\tstatus\tlabel\tclf_verdict\tclf_margin\tclade\n")
    for lid, cs in sorted(per.items()):
        r = info[lid]; c = clade(cs); lab = r["idiomorph_label"]
        if c in ("sexP", "sexM") and lab in ("Plus", "Minus", "undetermined") and not ((c == "sexP" and lab == "Plus") or (c == "sexM" and lab == "Minus")):
            fo.write("\t".join([lid, r["genome"], r["species"], r["family"], r["status"], lab, r["clf_verdict"], r["clf_margin"], c]) + "\n")
