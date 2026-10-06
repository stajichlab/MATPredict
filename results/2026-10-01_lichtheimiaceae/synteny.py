"""Neighbourhood synteny around candidate MAT genes in annotated Lichtheimiaceae genomes.

Inputs: genomes.tsv (family manifest), out/<gid>/{loci.tsv, clf_models.tsv, mp.tsv, clf_annot.tsv},
CDS03202.1.faa (L. ramosa SexM, Schulz 2016). For each annotated genome (gff3 + proteins):
  * anchor A = locus of the best tblastn hit of CDS03202.1 (literature sexM ortholog);
  * anchor B = best classifier-typed candidate (model or annotated protein, max MAT score).
Neighbour proteins within +-50 kb / +-100 kb of each anchor are written to neigh.faa; an
all-vs-all DIAMOND blastp (run outside) clusters them; conservation counts are computed per genus.
Also writes flank_positions.tsv: tptA/rnhA/glrA/algA/btbA best miniprot hit relative to anchor A.
"""
import csv, re, subprocess, sys
from pathlib import Path
from collections import defaultdict

E = "/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin"
HERE = Path(__file__).resolve().parent
SCR = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "tmp"
SCR.mkdir(parents=True, exist_ok=True)
GATE = 98.4


def read_tsv(p):
    return list(csv.DictReader(open(p), delimiter="\t")) if Path(p).exists() else []


def read_fasta(p):
    s, nm, b = {}, None, []
    for line in open(p):
        if line.startswith(">"):
            if nm: s[nm] = "".join(b)
            nm = line[1:].split()[0]; b = []
        else:
            b.append(line.strip())
    if nm: s[nm] = "".join(b)
    return s


def gff_genes(gff):
    """mRNA id -> (contig, start, end, strand, product)"""
    genes = {}
    for line in open(gff):
        if line.startswith("#") or "\tmRNA\t" not in line:
            continue
        f = line.rstrip("\n").split("\t")
        attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        genes[attrs.get("ID")] = (f[0], int(f[3]), int(f[4]), f[6], attrs.get("product", ""))
    return genes


rows = [r for r in read_tsv(HERE / "genomes.tsv") if r["set"] == "LCG" and r["gff"]]
neigh = open(HERE / "neigh.faa", "w")
anch_out = open(HERE / "anchors.tsv", "w")
anch_out.write("genome_id\tgenus\tspecies\tanchor\tcontig\tstart\tend\tdetail\n")
flank_out = open(HERE / "flank_positions.tsv", "w")
flank_out.write("genome_id\tgenus\tgene\tidentity\tcontig\tstart\tsame_contig_as_A\tdistance_to_A_kb\n")
for r in rows:
    gid = r["genome_id"]; od = HERE / "out" / gid
    if not (od / "done").exists():
        continue
    genes = gff_genes(r["gff"])
    prots = read_fasta(r["proteins"])
    # anchor A: tblastn CDS03202.1 vs genome
    db = SCR / gid
    subprocess.run([f"{E}/makeblastdb", "-in", r["genome"], "-dbtype", "nucl", "-out", str(db)], capture_output=True)
    t = subprocess.run([f"{E}/tblastn", "-query", str(HERE / "CDS03202.1.faa"), "-db", str(db), "-evalue", "1e-5", "-seg", "no",
                        "-outfmt", "6 sseqid pident length sstart send bitscore", "-max_target_seqs", "5"], capture_output=True, text=True)
    for f in db.parent.glob(db.name + ".*"):
        f.unlink()
    A = None
    for line in t.stdout.splitlines():
        s, pid, ln, ss, se, bs = line.split("\t")
        if A is None or float(bs) > A[3]:
            A = (s, min(int(ss), int(se)), max(int(ss), int(se)), float(bs), float(pid))
    anchors = {}
    if A:
        anchors["A_CDS03202"] = (A[0], A[1], A[2], f"bits={A[3]:.0f};pid={A[4]:.0f}")
    # anchor B: best typed candidate
    best = None
    for c in read_tsv(od / "clf_annot.tsv"):
        m, p, p1 = float(c["sexM"]), float(c["sexP"]), float(c["P1"])
        top = max(m, p)
        if top >= GATE and abs(m - p) >= 25 and not (p1 >= top + 25) and c["protein"] in genes:
            if best is None or top > best[0]:
                g = genes[c["protein"]]
                best = (top, g[0], g[1], g[2], f"annot {c['protein']} M={m} P={p} P1={p1}")
    loci = {x["mid"]: x for x in read_tsv(od / "loci.tsv")}
    for c in read_tsv(od / "clf_models.tsv"):
        m, p, p1 = float(c["sexM"]), float(c["sexP"]), float(c["P1"])
        top = max(m, p)
        L = loci.get(c["mid"])
        if L and top >= GATE and abs(m - p) >= 25 and not (p1 >= top + 25):
            if best is None or top > best[0]:
                best = (top, L["contig"], int(L["start"]), int(L["end"]), f"model {c['mid']} M={m} P={p} P1={p1}")
    if best:
        anchors["B_typed"] = best[1:]
    for an, (ctg, st, en, det) in anchors.items():
        anch_out.write(f"{gid}\t{r['genus']}\t{r['species']}\t{an}\t{ctg}\t{st}\t{en}\t{det}\n")
        for mid, (gc, gs, ge, strand, prod) in genes.items():
            if gc != ctg:
                continue
            d = 0 if ge >= st and gs <= en else min(abs(gs - en), abs(st - ge))
            if d <= 100000 and mid in prots:
                neigh.write(f">{gid}|{an}|{mid}|{d}|{r['genus']}\n{prots[mid]}\n")
    # flank genes relative to A
    best_f = {}
    for x in read_tsv(od / "mp.tsv"):
        gname = x["gene"]
        if gname in ("tptA", "rnhA", "glrA", "algA", "btbA"):
            idn = float(x["identity"] or 0)
            if gname not in best_f or idn > best_f[gname][0]:
                best_f[gname] = (idn, x["contig"], int(x["start"]))
    for gname, (idn, ctg, st) in best_f.items():
        same = bool(A and ctg == A[0])
        dist = f"{abs(st - A[1]) / 1000:.1f}" if same else ""
        flank_out.write(f"{gid}\t{r['genus']}\t{gname}\t{idn:.2f}\t{ctg}\t{st}\t{same}\t{dist}\n")
neigh.close(); anch_out.close(); flank_out.close()
print("done")
