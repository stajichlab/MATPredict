"""Protein models for the two out-of-group Plus gains (part 3).

For each gain: (a) translate detect's polished sexP model and its alternate
model from the candidate report; (b) re-model the locus with miniprot (--trans)
on +-8 kb using the 18 curated sexM/sexP proteins, the Circinella-group queries
and the two new records' core proteins, keeping the best AS; (c) for LCG, the
funannotate protein overlapping the locus. Score each with the baseline
(polish-scope-cuts) and candidate (curation-circinella) classifiers.
Writes gains.faa (best full-length model per gain) and gains.tsv.
"""
import gzip, re, subprocess, sys
from pathlib import Path
import pyhmmer, yaml
from Bio.Seq import Seq

E = "/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin"
R = Path("/bigdata/stajichlab/jstajich/projects/MATPredict/results")
WT = Path("/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees")
LT = R / "2026-10-01_circinella_label_tree"
SCR = Path(sys.argv[1])
CLF = {"base": WT / "polish-scope-cuts/db/Mucoromycota/classifiers/MAT",
       "cand": WT / "curation-circinella/db/Mucoromycota/classifiers/MAT"}
LCG = "/bigdata/stajichlab/shared/projects/ZyGoLife/LCG/Annotation"
GAINS = [
    dict(gid="GCA_025677815.1_Phaart1", genome=SCR / "Phaart1.fna",
         report=SCR / "phaart_cand/GCA_025677815.1_Phaart1/detection_report.yaml", gff=None, prot=None),
    dict(gid="Absidia_sp._NRRL_3163", genome=Path(f"{LCG}/genomes/Absidia_sp._NRRL_3163.sorted.fasta"),
         report=R / "2026-10-01_circinella_curation/lcg/cand/Absidia_sp._NRRL_3163/detection_report.yaml",
         gff=f"{LCG}/annotate/Absidia_sp._NRRL_3163/annotate_results/Absidia_sp._NRRL_3163.gff3",
         prot=f"{LCG}/annotate/Absidia_sp._NRRL_3163/annotate_results/Absidia_sp._NRRL_3163.proteins.fa"),
]
alpha = pyhmmer.easel.Alphabet.amino()
H = {}
for side, d in CLF.items():
    for n, p in (("sexM", d / "sexM.hmm"), ("sexP", d / "sexP.hmm"), ("P1", d / "paralogs/P1.hmm")):
        with pyhmmer.plan7.HMMFile(str(p)) as fh:
            H[(side, n)] = fh.read()
pf = pyhmmer.plan7.HMMFile("/srv/projects/db/pfam/2024-06-04/Pfam-A.hmm") if False else None


def score(seq):
    s = seq.replace("*", "").replace("X", "")
    blk = pyhmmer.easel.DigitalSequenceBlock(alpha, [pyhmmer.easel.TextSequence(name=b"q", sequence=s).digitize(alpha)])
    out = {}
    for k, h in H.items():
        best = 0.0
        for top in pyhmmer.hmmsearch([h], blk, T=0.0, domT=0.0):
            for hit in top:
                best = max(best, hit.score)
        out[f"{k[0]}_{k[1]}"] = round(best, 1)
    return out


def fasta(path):
    op = gzip.open if str(path).endswith(".gz") else open
    seqs, nm = {}, None
    for line in op(path, "rt"):
        if line.startswith(">"):
            nm = line[1:].split()[0]; seqs[nm] = []
        else:
            seqs[nm].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def contig(path, name):
    op = gzip.open if str(path).endswith(".gz") else open
    on, buf = False, []
    for line in op(path, "rt"):
        if line.startswith(">"):
            if on: break
            on = line[1:].split()[0] == name
        elif on:
            buf.append(line.strip())
    return "".join(buf)


def translate(cs, exons, strand):
    nt = "".join(cs[s - 1:e] for s, e in sorted(exons))
    if strand == "-":
        nt = str(Seq(nt).reverse_complement())
    return str(Seq(nt[: len(nt) // 3 * 3]).translate())


queries = SCR / "gain_queries.faa"
with open(queries, "w") as q:
    for p in (R / "2026-10-01_lichtheimiaceae/sexMP_refs.faa", LT / "circ_queries2.faa"):
        q.write(open(p).read())
    for rec, gene in (("Mucorales/101103_nrrl1351_MAT_Plus", "sexP"), ("Mucorales/64656_rsa-1403_MAT_Minus", "sexM")):
        for k, v in fasta(WT / f"curation-circinella/db/Mucoromycota/{rec}/proteins.faa").items():
            if f"name={gene}|" in k + "|":
                q.write(f">REC_{rec.split('/')[1]}|{k}\n{v}\n")

fa = open("gains.faa", "w"); tab = open("gains.tsv", "w")
cols = ["gid", "model", "coords", "len", "base_sexM", "base_sexP", "base_P1", "cand_sexM", "cand_sexP", "cand_P1", "note"]
tab.write("\t".join(cols) + "\n")
for g in GAINS:
    rep = yaml.safe_load(open(g["report"]))
    loc = [d for d in rep["detected"] if d["family"] == "Mucoromycota:MAT"][0]
    sp = [e for e in loc["gene_evidence"] if e["gene"] == "sexP"][0]
    cs = contig(g["genome"], sp["contig"])
    models = {}
    models["detect_polished"] = (translate(cs, [(x["start"], x["end"]) for x in sp["exons"]], sp["strand"]),
                                 f"{sp['contig']}:{sp['start']}-{sp['end']}{sp['strand']}", sp["method"])
    alt = sp["alternate_model"]
    models["detect_alternate"] = (translate(cs, [(x["start"], x["end"]) for x in alt["exons"]], alt["strand"]),
                                  f"{alt['contig']}:{alt['start']}-{alt['end']}{alt['strand']}", alt["method"])
    lo, hi = max(0, sp["start"] - 8000), min(len(cs), sp["end"] + 8000)
    reg = SCR / f"{g['gid']}_reg.fa"; reg.write_text(f">r\n{cs[lo:hi]}\n")
    p = subprocess.run([f"{E}/miniprot", "--trans", "-p", "0.3", "-N", "5", str(reg), str(queries)],
                       capture_output=True, text=True, check=True)
    best, cur = None, None
    for line in p.stdout.splitlines():
        if not line.startswith("#") and "\t" in line:
            f = line.split("\t")
            ms = [x for x in f if x.startswith("AS:i:")]
            cur = dict(q=f[0], ts=lo + int(f[7]) + 1, te=lo + int(f[8]), st=f[4], AS=int(ms[0][5:]))
        elif line.startswith("##STA") and cur:
            cur["prot"] = line.split("\t")[1].strip()
            if best is None or cur["AS"] > best["AS"]:
                best = cur
            cur = None
    models["miniprot_best"] = (best["prot"], f"{sp['contig']}:{best['ts']}-{best['te']}{best['st']}", f"query {best['q']} AS {best['AS']}")
    if g["gff"]:
        ov = (0, None)
        for line in open(g["gff"]):
            f = line.rstrip("\n").split("\t")
            if len(f) == 9 and f[0] == sp["contig"] and f[2] == "mRNA":
                o = min(int(f[4]), sp["end"]) - max(int(f[3]), sp["start"])
                if o > ov[0]:
                    ov = (o, re.search(r"ID=([^;]+)", f[8]).group(1), f"{f[0]}:{f[3]}-{f[4]}{f[6]}")
        if ov[1]:
            models["funannotate"] = (fasta(g["prot"])[ov[1]], ov[2], ov[1])
    for name, (seq, coords, note) in models.items():
        s = score(seq)
        tab.write("\t".join(str(x) for x in [g["gid"], name, coords, len(seq.rstrip("*")), s["base_sexM"], s["base_sexP"], s["base_P1"],
                                              s["cand_sexM"], s["cand_sexP"], s["cand_P1"], note]) + "\n")
    # tree tip: detect's alternate model, the protein the classifier scored
    n, m = "detect_alternate", models["detect_alternate"]
    assert "*" not in m[0].rstrip("*"), g["gid"]
    fa.write(f">GAIN__{g['gid']}__{n}\n{m[0].rstrip('*')}\n")
    tab.write(f"# tree tip for {g['gid']}: {n}\n")
fa.close(); tab.close()
