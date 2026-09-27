"""Pick a diverse proteome sample: per-order caps, known model species forced in."""
import csv, random, collections
random.seed(20260927)
rows = [r for r in csv.DictReader(open("proteomes.tsv"), delimiter="\t") if r["proteome"]]
CAP = {"Agaricales": 25, "Boletales": 15, "Polyporales": 15, "Russulales": 12, "Hymenochaetales": 10,
       "Cantharellales": 10, "Tremellales": 15, "Sporidiobolales": 15, "Pucciniales": 12, "Ustilaginales": 12,
       "Wallemiales": 6, "Malasseziales": 6, "Microbotryales": 6, "Tilletiales": 6, "Trichosporonales": 6,
       "Auriculariales": 4, "Sebacinales": 3, "Filobasidiales": 4, "Cystofilobasidiales": 4}
FORCE = ["Coprinopsis cinerea", "Schizophyllum commune", "Ustilago maydis", "Cryptococcus neoformans",
         "Rhodotorula toruloides", "Puccinia graminis", "Wallemia mellicola", "Laccaria bicolor",
         "Microbotryum lychnidis-dioicae", "Sporisorium reilianum", "Rhodotorula mucilaginosa",
         "Phanerochaete chrysosporium", "Heterobasidion annosum", "Agaricus bisporus"]
by = collections.defaultdict(list)
for r in rows: by[r["ORDER"]].append(r)
pick = {}
for sp in FORCE:
    c = [r for r in rows if r["SPECIES"].strip() == sp]
    if c: pick[c[0]["ASMID"]] = c[0]
for o, lst in by.items():
    cap = CAP.get(o, 2)
    have = sum(1 for r in pick.values() if r["ORDER"] == o)
    random.shuffle(lst)
    seen_sp = {r["SPECIES"] for r in pick.values()}
    for r in lst:  # one per species first
        if have >= cap: break
        if r["SPECIES"] in seen_sp: continue
        pick[r["ASMID"]] = r; seen_sp.add(r["SPECIES"]); have += 1
with open("sample.tsv", "w") as fo:
    w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(pick.values())
print(len(pick), collections.Counter(r["ORDER"] for r in pick.values()).most_common(12))
print("forced found:", [sp for sp in FORCE if any(r["SPECIES"].strip()==sp for r in pick.values())])
