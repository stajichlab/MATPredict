"""Per-genome table of miniprot hits of MAT-related reference proteins: best hit per gene class and its coordinates."""
import re, sys, glob, collections
for fn in sorted(glob.glob("data/*.miniprot.gff")):
    rows = []
    for l in open(fn):
        f = l.rstrip("\n").split("\t")
        if len(f) > 8 and f[2] == "mRNA":
            t = re.search(r"Target=(\S+)", f[8]).group(1)
            idn = float(re.search(r"Identity=([\d.]+)", f[8]).group(1)); pos = float(re.search(r"Positive=([\d.]+)", f[8]).group(1))
            m = re.search(r"name=([^|\s;]+)", t); gene = m.group(1) if m else t.split("|")[-1].split("_")[-1]
            rows.append((f[0], int(f[3]), int(f[4]), f[6], gene, idn, t))
    print("==", fn.split("/")[-1][:30])
    by = collections.defaultdict(list)
    for r in rows: by[(r[0], r[4])].append(r)
    for (c, g), rs in sorted(by.items(), key=lambda kv: (kv[0][0], min(x[1] for x in kv[1]))):
        b = max(rs, key=lambda x: x[5])
        if b[5] >= 0.3 or g in ("ONH67106", "XP_020071302"):
            print(f"  {c}\t{b[1]}-{b[2]}\t{b[3]}\t{g}\tid={b[5]:.2f}\tn={len(rs)}")
