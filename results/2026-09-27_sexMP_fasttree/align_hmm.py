"""Step 5: hmmalign the HMG-box domains to PF00505 (match columns only), as in
results/2026-09-26_sexMP_phylogeny (alignment A), and write short ids for the tree."""
import subprocess
def read_afa(f):
    d, n = {}, None
    for l in open(f):
        l = l.rstrip()
        if l.startswith(">"): n = l[1:].split()[0]; d[n] = []
        elif n: d[n].append(l)
    return {k: "".join(v) for k, v in d.items()}
subprocess.run("hmmalign --trim --outformat afa PF00505.hmm domains.faa > aln_hmm_raw.afa", shell=True, check=True)
raw = read_afa("aln_hmm_raw.afa")
L = len(next(iter(raw.values())))
keep = [i for i in range(L) if all(not (s[i].islower() or s[i] == ".") for s in raw.values())]
A = {k: "".join(v[i] for i in keep).upper() for k, v in raw.items()}
with open("aln_hmm.afa", "w") as fa, open("aln_hmm.ids.fa", "w") as fi, open("taxon_ids.tsv", "w") as ft:
    ft.write("id\tname\n")
    for i, (k, v) in enumerate(A.items()):
        t = f"T{i:04d}"
        fa.write(f">{k}\n{v}\n"); fi.write(f">{t}\n{v}\n"); ft.write(f"{t}\t{k}\n")
print("sequences", len(A), "columns", len(keep))
