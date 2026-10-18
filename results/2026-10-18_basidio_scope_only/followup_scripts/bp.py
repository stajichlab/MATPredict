import csv, collections
B = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/ncbi_accessions.csv"
bp = {r["ASM_FOLDER"]: (r["BIOPROJECT"], r["N50"], r["ASM_LENGTH"]) for r in csv.DictReader(open(B))}
rows = list(csv.DictReader(open("/bigdata/stajichlab/jstajich/projects/mat-fam18/pr/base.tsv"), delimiter="\t"))
mu = [r for r in rows if r["species"] == "Rhodotorula mucilaginosa"]
byp = collections.defaultdict(collections.Counter)
for r in mu:
    folder = r["genome"].split("_", 1)[1]
    p = bp.get(folder, ("?",))[0]
    byp[p][r["label"] or "none"] += 1
print("R. mucilaginosa genomes", len({r['genome'] for r in mu}), "bioprojects", len(byp))
for p, c in sorted(byp.items(), key=lambda t: -sum(t[1].values()))[:8]: print("  ", p, dict(c))
# all genomes: bioprojects overall for sporidio
allp = collections.defaultdict(collections.Counter)
for r in rows:
    allp[bp.get(r["genome"].split("_", 1)[1], ("?",))[0]][r["label"] or "none"] += 1
print("all Sporidiobolales bioprojects", len(allp))
print("species-level unique bioproject-species pairs A1/A2:")
pairs = collections.defaultdict(set)
for r in rows:
    if r["label"]:
        pairs[(r["species"], bp.get(r["genome"].split("_", 1)[1], ("?",))[0])].add(r["label"])
c = collections.Counter("both" if len(v) > 1 else next(iter(v)) for v in pairs.values())
print("  (species,bioproject) units", len(pairs), dict(c))
