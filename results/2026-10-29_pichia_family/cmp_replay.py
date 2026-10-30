"""Before/after of adding the MATpk family: the 158 Pichia genomes of the v0.7.0 Ascomycota run.
Both runs use `--preset v0.8`; `before` is tree_v5 (no MATpk), `after` is tree_v6 (MATpk, scope genus Pichia 4919).
usage: python3 cmp_replay.py <yeast_routing dir> <list.tsv> [after-run directory name]   (run where /bigdata is mounted)"""
import collections, glob, os, sys

D, LIST = sys.argv[1], sys.argv[2]
AFTER = sys.argv[3] if len(sys.argv) > 3 else "out_pichia_after"
sp = {}
for l in open(LIST):
    g, t = l.rstrip("\n").split("\t")[:2]
    sp[g] = t

def loci(run, g):
    f = os.path.join(D, run, "runs", g, "detected_loci.gff3")
    if not os.path.exists(f):
        return None
    out = []
    for l in open(f):
        if l.startswith("#"):
            continue
        c = l.rstrip("\n").split("\t")
        if len(c) > 8 and c[2] == "MAT_locus":
            a = dict(x.split("=", 1) for x in c[8].split(";") if "=" in x)
            out.append((a.get("family", "?").split(":")[-1], c[0], int(c[3]), int(c[4]), a.get("idiomorph", ""), a.get("confidence", ""),
                        a.get("family_status", ""), a.get("cassette", "")))
    return out

names = {}
for l in open(os.path.join(D, "pichia_names.tsv")) if os.path.exists(os.path.join(D, "pichia_names.tsv")) else []:
    g, n = l.rstrip("\n").split("\t"); names[g] = n
rows = []
for g in sp:
    b, a = loci("out_pichia_before", g), loci(AFTER, g)
    rows.append((g, b, a))
def kind(x):
    if x is None: return "no_report"
    return "called" if x else "uncalled"
print("status before -> after:", collections.Counter((kind(b), kind(a)) for g, b, a in rows))
print("families before:", collections.Counter(tuple(sorted({x[0] for x in (b or [])})) for g, b, a in rows))
print("families after :", collections.Counter(tuple(sorted({x[0] for x in (a or [])})) for g, b, a in rows))
print("confidence after (MATpk):", collections.Counter(x[5] for g, b, a in rows for x in (a or []) if x[0] == "MATpk"))
print()
for g, b, a in rows:
    sb = ";".join(f"{x[0]}:{x[1]}:{x[2]}-{x[3]}:{x[4]}:{x[5]}" for x in (b or []))
    sa = ";".join(f"{x[0]}:{x[1]}:{x[2]}-{x[3]}:{x[4]}:{x[5]}{'/'+x[6] if x[6] else ''}" for x in (a or []))
    print(g, sp[g], "|", sb or kind(b), "=>", sa or kind(a))
