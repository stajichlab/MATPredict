"""Step 4: cut the HMG box from every candidate, outgroup and reference; deduplicate.

hmmsearch Pfam PF00505 (HMG_box, Pfam 38.2) over all proteins. Keep the best domain
per protein with i-Evalue <= 0.01, envelope extended 5 aa each side.
Identical domain sequences are collapsed; dedup_map.tsv lists every member.
"""
import collections, glob, subprocess
seqs = {}
def read(f):
    n = None
    for l in open(f):
        if l[0] == ">": n = l[1:].split()[0]; seqs[n] = ""
        else: seqs[n] += l.strip()
for f in sorted(glob.glob("extract2/*.faa")) + ["refs_sexMP.faa", "refs_MAT121.faa"]:
    read(f)
with open("all_proteins.faa", "w") as fo:
    for n, s in seqs.items():
        if s: fo.write(f">{n}\n{s}\n")
subprocess.run("hmmsearch --cpu 4 -E 10 --domE 10 --domtblout all_pf00505.domtbl PF00505.hmm all_proteins.faa > /dev/null", shell=True, check=True)
best = {}
for l in open("all_pf00505.domtbl"):
    if l[0] == "#": continue
    f = l.split()
    n, iE, a, b = f[0], float(f[12]), int(f[19]), int(f[20])
    if iE <= 0.01 and (n not in best or iE < best[n][0]):
        best[n] = (iE, a, b)
dom = {}
for n, (iE, a, b) in best.items():
    s = seqs[n]; dom[n] = s[max(0, a - 1 - 5): min(len(s), b + 5)]
groups = collections.defaultdict(list)
for n, d in dom.items():
    groups[d].append(n)
with open("domains.faa", "w") as fo, open("dedup_map.tsv", "w") as fm:
    fm.write("rep\tmember\n")
    for d, ms in groups.items():
        ms.sort(key=lambda x: (not x.startswith(("REF|", "OUT|")), x))
        rep = ms[0]
        fo.write(f">{rep}\n{d}\n")
        for m in ms: fm.write(f"{rep}\t{m}\n")
nodom = [n for n in seqs if n not in best]
open("no_hmg_domain.txt", "w").write("\n".join(nodom) + "\n")
print("proteins", len(seqs), "with HMG box", len(best), "unique domains", len(groups), "no domain", len(nodom))
