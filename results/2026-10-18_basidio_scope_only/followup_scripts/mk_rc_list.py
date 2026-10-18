import csv, glob, random, yaml
W = "/bigdata/stajichlab/jstajich/projects/mat-fam18/"
J = {r["acc"]: r for r in csv.DictReader(open(W + "bfd_taxonomy_joined.tsv"), delimiter="\t")}
lst = {l.split("\t")[0]: l.rstrip("\n") for l in open(W + "lists/pucc_all.txt")}
called = {}
for p in glob.glob(W + "runs/pucc/*/detection_report.yaml"):
    d = yaml.load(open(p), Loader=yaml.CSafeLoader)
    acc = p.split("/")[-2].split("_", 1)[1]
    called[acc] = any(x["family"].endswith("rustHD") for x in d.get("detected") or [])
targets = {a for a in called if J[a]["tk_genus"] in ("Hemileia", "Endocronartium", "Uromyces", "Melampsora")}
unc = {a for a, c in called.items() if not c}
rng = random.Random(1)
ctrl = sorted(a for a, c in called.items() if c and a not in targets)
ctrl = set(rng.sample(ctrl, 12))
sel = sorted(targets | unc | ctrl)
print(len(targets), len(unc), len(ctrl), len(sel))
for i in range(10):
    open(W + f"lists/rc.{i:02d}", "w").write("".join(lst[a] + "\n" for a in sel[i::10]))
with open(W + "rc/groups.tsv", "w") as fh:
    for a in sel: fh.write(f"{a}\t{J[a]['tk_genus']}\t{'called' if called[a] else 'uncalled'}\t{'target' if a in targets else 'uncalled_other' if a in unc else 'control'}\n")
