"""Flank genes around each called locus of the ten two-family genomes (run on the HPCC): miniprot of flank_refs_all.faa on locus +/- 50 kb.
Output flank10.tsv: genome, family, idiomorph, contig, locus, flank hits as gene:identity%@start-end(strand) relative to the locus start."""
import gzip, subprocess, re, sys, os
LIB = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes/"
MP = "/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin/miniprot"
HERE = os.path.dirname(os.path.abspath(__file__))
def seq_of(path, contig):
    n, s = None, []
    for l in gzip.open(path, "rt"):
        if l[0] == ">":
            if n == contig: return "".join(s)
            n = l[1:].split()[0]; s = []
        else: s.append(l.strip())
    return "".join(s) if n == contig else None
out = ["genome\tfamily\tidiomorph\tcontig\tlocus\tcontig_len\tflank_hits"]
for l in open(HERE + "/loci10.tsv"):
    g, fam, ctg, a, b, idio = l.rstrip().split("\t"); a, b = int(a), int(b)
    s = seq_of(LIB + g + ".fa.gz", ctg); lo = max(1, a - 50000); hi = min(len(s), b + 50000)
    open("/tmp/w189.fa", "w").write(">w\n" + s[lo - 1:hi] + "\n")
    res = subprocess.run([MP, "--gff", "-t2", "/tmp/w189.fa", HERE + "/flank_refs_all.faa"], capture_output=True, text=True).stdout
    best = {}
    for ln in res.splitlines():
        f = ln.split("\t")
        if len(f) > 8 and f[2] == "mRNA":
            at = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
            tgt = at["Target"].split()[0]; m = re.search(r"name=([^|]+)", tgt); nm = (m.group(1) if m else tgt.rsplit("_", 1)[-1]).lower()
            if nm in ("mata1", "mata2", "matalpha1", "matalpha2", "alpha1", "alpha2", "a1", "a2"): continue
            idn = float(at["Identity"]) * 100; st, en = int(f[3]) + lo - 1, int(f[4]) + lo - 1
            if nm not in best or int(f[5]) > best[nm][0]: best[nm] = (int(f[5]), idn, st, en, f[6])
    hits = ";".join(f"{k}:{v[1]:.0f}%@{v[2]-a:+d}({v[4]})" for k, v in sorted(best.items(), key=lambda kv: kv[1][2]) if v[0] > 150 and v[1] >= 35)
    out.append(f"{g}\t{fam}\t{idio}\t{ctg}\t{a}-{b}\t{len(s)}\t{hits}")
open(HERE + "/flank10.tsv", "w").write("\n".join(out) + "\n"); print(len(out) - 1)
