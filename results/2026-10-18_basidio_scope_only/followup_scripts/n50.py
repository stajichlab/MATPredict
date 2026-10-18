import csv, glob, yaml, collections, statistics
W = "/bigdata/stajichlab/jstajich/projects/mat-fam18/runs/"
B = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/ncbi_accessions.csv"
n50 = {r["ASM_FOLDER"]: int(r["N50"] or 0) for r in csv.DictReader(open(B))}
bins = [(0, 10_000, "<10 kb"), (10_000, 100_000, "10-100 kb"), (100_000, 1_000_000, "0.1-1 Mb"), (1_000_000, 10**12, ">=1 Mb")]
out = []
for run, fams in (("pucc", ["rustHD"]), ("sporidio", ["redHD", "redPR"]), ("walle", ["wallMAT"])):
    tab = {f: collections.defaultdict(lambda: [0, 0]) for f in fams}
    for p in glob.glob(W + run + "/*/detection_report.yaml"):
        d = yaml.load(open(p), Loader=yaml.CSafeLoader)
        folder = p.split("/")[-2].split("_", 1)[1]; n = n50.get(folder, 0)
        lab = next(b[2] for b in bins if b[0] <= n < b[1])
        for f in fams:
            c = any(x["family"].endswith(f) for x in d.get("detected") or [])
            tab[f][lab][0] += 1; tab[f][lab][1] += c
    for f in fams:
        out.append(f"{f}: " + "; ".join(f"{b[2]} {tab[f][b[2]][1]}/{tab[f][b[2]][0]}" for b in bins))
open("/bigdata/stajichlab/jstajich/projects/mat-fam18/rc/n50_bins.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
