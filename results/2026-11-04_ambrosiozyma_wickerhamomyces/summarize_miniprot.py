"""Best miniprot hit per gene name and 3 kb bin, for each window (default and relaxed runs, refs.faa)."""
import glob, re
def gname(tgt):
    m = re.search(r"name=([^|]+)", tgt)
    n = m.group(1) if m else re.split(r"[_|]", tgt)[-1]
    return n.upper()
for f in sorted(glob.glob("data/*.miniprot.gff")):
    print("==", f.split("/")[1].split(".miniprot")[0])
    best = {}
    for fn in (f, f.replace(".miniprot.gff", ".relaxed.gff")):
        for l in open(fn):
            if l.startswith("#") or "\tmRNA\t" not in l: continue
            c = l.rstrip().split("\t"); a = dict(kv.split("=", 1) for kv in c[8].split(";"))
            tgt = a["Target"].split()[0]; g = gname(tgt)
            k = (g, int(c[3]) // 3000, c[6])
            r = (int(c[3]), int(c[4]), c[6], float(a["Identity"]), int(c[5]), "Frameshift" in a, tgt[:30])
            if k not in best or r[4] > best[k][4]: best[k] = r
    for k, r in sorted(best.items(), key=lambda x: x[1][0]):
        if r[4] >= 100: print("  %-10s %6d-%6d %s id=%.2f score=%d fs=%s %s" % (k[0], r[0], r[1], r[2], r[3], r[4], r[5], r[6]))
