#!/usr/bin/env python3
"""Leave-own-genome-out gene models for the receptor proteins.

The first pass (build_genome.py) modelled each locus with the non-curated queries only. Many
models were fragments. Here each genome is modelled with every STE3 query (curated REF
receptors included) except the genome's own curated record(s) and the queries that came from the
same genome. A mating receptor is then modelled from its relatives, not from its own record, as
a pipeline would on a new genome. Updates prot_len, n_cds, aln_score in OUT/<asm>.loci.tsv and
rewrites OUT/<asm>.prot.faa.  Usage: reprotein.py ASM OUTDIR
"""
import csv
import importlib.util
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("build_genome_mod", os.path.join(HERE, "build_genome.py"))
src = open(os.path.join(HERE, "build_genome.py")).read().replace("\nmain()\n", "\n")
ns = {"__file__": os.path.join(HERE, "build_genome.py"), "__name__": "bg"}
exec(compile(src, "build_genome.py", "exec"), ns)
parse_gff, BIN, LIB = ns["parse_gff"], ns["BIN"], ns["LIB"]
EXCL = {
    "GCA_016772295.1": ["REF|5346_"], "GCF_000143185.2": ["REF|5334_"], "GCA_000988875.2": ["REF|5286_nbrc-0880"],
    "GCA_921037615.3": ["REF|5286_cbs-14_"], "GCF_000328475.2": ["REF|5270_"], "GCF_000091045.1": ["REF|40410_jec21"],
    "GCA_056621545.1": ["REF|40410_jec20"], "GCA_026119225.1": ["REF|29898_"], "GCA_920103745.3": ["REF|5535_"],
    "GCA_024748845.1": ["REF|5537_"], "GCA_023212685.1": ["REF|1652704_"], "GCA_023212835.1": ["REF|203535_"],
    "GCA_023212605.1": ["REF|203536_"], "GCA_023212725.1": ["REF|349360_"], "GCA_023212615.1": ["REF|49012_"],
    "GCA_023212695.1": ["REF|63387_"], "GCA_023212635.2": ["REF|84751_"], "GCA_056320075.1": ["REF|1708542_"],
    "GCF_000263375.1": ["REF|671144_"],
}
asm, out = sys.argv[1], sys.argv[2]
pre = asm[:15]
q = os.path.join(out, "tmp_" + asm + ".q.faa")
keep = False
with open(q, "w") as fo:
    for line in open(os.path.join(HERE, "ste3_all.faa")):
        if line.startswith(">"):
            h = line[1:]
            keep = not (h.startswith(pre) or any(h.startswith(x) for x in EXCL.get(pre, [])))
        if keep:
            fo.write(line)
gff = os.path.join(out, "tmp_" + asm + ".gff")
with open(gff, "w") as fo:
    subprocess.run([f"{BIN}/miniprot", "-t", os.environ.get("SLURM_CPUS_PER_TASK", "8"), "-I", "--gff", "--trans", "--outn=100",
                    f"{LIB}/{asm}.fa.gz", q], stdout=fo, stderr=subprocess.DEVNULL, check=True)
mr = parse_gff(gff)
rows = list(csv.DictReader(open(os.path.join(out, asm + ".loci.tsv")), delimiter="\t"))
fp = open(os.path.join(out, asm + ".prot.faa"), "w")
for r in rows:
    a, b = int(r["start"]), int(r["end"])
    ov = [m for m in mr if m["contig"] == r["contig"] and m["start"] <= b and m["end"] >= a and m["prot"]]
    m = max(ov, key=lambda m: m["score"]) if ov else None
    r["prot_len"], r["n_cds"], r["aln_score"] = (len(m["prot"]), m["ncds"], m["score"]) if m else (0, 0, 0)
    if m:
        fp.write(f">{asm}|{r['contig']}:{a}-{b}\n{m['prot']}\n")
fp.close()
with open(os.path.join(out, asm + ".loci.tsv"), "w") as fo:
    w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(rows)
os.remove(q)
os.remove(gff)
print(asm, "reprotein done", sum(1 for r in rows if int(r["prot_len"]) > 0), "/", len(rows))
