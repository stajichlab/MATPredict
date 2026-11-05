"""Best miniprot hit per reference gene and 2 kb bin, for each window (refs.faa: curated Saccharomycotina set + the two MATog records)."""
import glob, re, sys
for f in sorted(glob.glob("data/*.miniprot.gff")):
    print("==", f.split("/")[1].split(".miniprot")[0])
    best = {}
    for l in open(f):
        if l.startswith("#") or "\tmRNA\t" not in l: continue
        c = l.rstrip().split("\t"); a = dict(kv.split("=", 1) for kv in c[8].split(";"))
        tgt = a["Target"].split()[0]
        m = re.search(r"name=([^|]+)", tgt); g = m.group(1) if m else tgt.split("_", 1)[-1]
        src = tgt.split("|")[0]
        r = (int(c[3]), int(c[4]), c[6], float(a.get("Identity", 0)), int(c[5]), src)
        k = (g, r[0] // 2000)
        if k not in best or r[4] > best[k][4]: best[k] = r
    for k, r in sorted(best.items(), key=lambda x: x[1][0]):
        print("  %-10s %6d-%6d %s id=%.2f score=%d  %s" % (k[0], r[0], r[1], r[2], r[3], r[4], r[5][:34]))
