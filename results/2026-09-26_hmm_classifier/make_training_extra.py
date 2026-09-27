"""Select the extra sexP training proteins for the Mucoromycota:MAT classifier.

Rule (curator's option (a), 2026-09-26): members of the UFBoot-99 sexP clade
of the 2026-09-26 HMG-box tree (results/2026-09-26_sexMP_phylogeny), minus
  * curated references (already in the build from the database),
  * every protein in the disputed sets (kept as test material): Minus calls
    whose HMG box is in the sexP clade, and Lichtheimiaceae /
    Syncephalastraceae sexP-clade loci,
  * any genome whose strain is one of the Zygo 23 truth organisms (Zygo stays
    an independent held-out test),
  * proteins shorter than 50 aa, and exact duplicates.
sexM gets NO extra members: the sexM group is not supported (UFBoot 57).
Writes training_extra.faa (headers `>{id}|sexP|genus={genus}`) and
training_extra.tsv (why each candidate was kept or dropped).
"""
import csv, re, collections

PH = "../2026-09-26_sexMP_phylogeny"
SAMPLES = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"
ZYGO = "../2026-09-23_zygo_bar/zygo_truth.tsv"


def read_fa(path):
    seqs, cur = {}, None
    for l in open(path):
        if l.startswith(">"):
            cur = l[1:].split()[0]; seqs[cur] = []
        elif cur:
            seqs[cur].append(l.strip())
    return {k: "".join(v).replace("*", "") for k, v in seqs.items()}


def norm(s):
    return re.sub(r"[^a-z0-9]", "", s.lower())


allp = read_fa(f"{PH}/all_proteins.faa")
cand = list(csv.DictReader(open(f"{PH}/candidates.tsv"), delimiter="\t"))
samp = {r["ASMID"]: r for r in csv.DictReader(open(SAMPLES))}
zygo_orgs = {l.split("\t")[0] for l in open(ZYGO)}
# a Zygo organism key: genus+epithet + strain digits/letters, normalised
zkeys = set()
for o in zygo_orgs:
    o = re.sub(r"_(Plus|Minus|-)$", "", o)
    parts = o.split("_")
    strain = "".join(p for p in parts if re.search(r"\d", p))
    zkeys.add((norm(" ".join(parts[:2])), norm(strain)))

disputed_fams = {"Lichtheimiaceae", "Syncephalastraceae"}
rows, keep, seen = [], {}, set()
for r in cand:
    if r["tree_clade"] != "sexP":
        continue
    cid, g = r["candidate_id"], r["genome"]
    s = samp.get(g, {})
    sp = s.get("SPECIES") or r["species"] or g
    strain = s.get("STRAIN", "")
    genus = sp.split()[0]
    why = "kept"
    seq = allp.get(cid, "")
    if r["status"] == "called" and r["idiomorph_label"] == "Minus":
        why = "disputed:minus_in_sexP_clade"
    elif r["family"] in disputed_fams:
        why = "disputed:lichth_syncephal"
    elif (norm(" ".join(sp.split()[:2])), norm("".join(re.findall(r"[A-Za-z]*\d+[A-Za-z0-9-]*", strain)))) in zkeys:
        why = "zygo_strain"
    elif len(seq) < 50:
        why = "short_or_missing"
    elif seq in seen:
        why = "duplicate"
    rows.append(dict(candidate_id=cid, genome=g, species=sp, strain=strain, family=r["family"],
                     status=r["status"], label=r["idiomorph_label"], length=len(seq), decision=why))
    if why == "kept":
        seen.add(seq)
        keep[f"{cid.replace('|', ':')}|sexP|genus={genus}"] = seq

with open("training_extra.faa", "w") as fo:
    for k, v in keep.items():
        fo.write(f">{k}\n{v}\n")
with open("training_extra.tsv", "w", newline="") as fo:
    w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader(); w.writerows(rows)
print("sexP-clade candidates:", len(rows), collections.Counter(r["decision"] for r in rows))
print("kept:", len(keep), "genera:", len({k.split("genus=")[1] for k in keep}))
