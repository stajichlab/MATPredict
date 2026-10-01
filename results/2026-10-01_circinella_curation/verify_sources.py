"""Verify the source proteins of the two Circinella-group records against the
assembly: translate the annotated CDS exons from the NCBI contig and compare to
the deposited protein (100% identity required); check transl_table and phase;
find the PF00505 HMG box; score with the current classifier.

Usage: python verify_sources.py SCRATCH_DIR CLASSIFIER_DIR PF00505.hmm > verify_sources.tsv
"""
import re, sys
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
import pyhmmer

scr, clf_dir, pfam = map(Path, sys.argv[1:4])
SRC = {
    "GCF_025766255.1": ("NW_026516701.1", ["XP_052979470.1", "XP_052979471.1", "XP_052979472.1",
                                          "XP_052979473.1", "XP_052979474.1", "XP_052979475.1"]),
    "GCA_025093555.1": ("JAIWNF010000086.1", ["KAI7847719.1", "KAI7847720.1", "KAI7847721.1",
                                             "KAI7847722.1", "KAI7847723.1"]),
}
alpha = pyhmmer.easel.Alphabet.amino()
hmms = {n: pyhmmer.plan7.HMMFile(clf_dir / f"{n}.hmm").read() for n in ("sexM", "sexP")}
pf = pyhmmer.plan7.HMMFile(pfam).read()
print("assembly\tcontig\tprotein\tstrand\tstart\tend\tn_exons\texons\tphase0\ttransl_table\tlen_aa\tidentity\tPF00505\tsexM\tsexP")
for asm, (contig, prots) in SRC.items():
    d = scr / asm / "ncbi_dataset/data" / asm
    seq = str(next(SeqIO.parse(scr / f"{contig}.fa", "fasta")).seq).upper()
    faa = {r.id: str(r.seq) for r in SeqIO.parse(d / "protein.faa", "fasta")}
    cds = {}
    for line in open(d / "genomic.gff"):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[0] != contig or f[2] != "CDS":
            continue
        m = re.search(r"protein_id=([^;]+)", f[8])
        if m and m[1] in prots:
            tt = re.search(r"transl_table=(\d+)", f[8])
            cds.setdefault(m[1], []).append((int(f[3]), int(f[4]), f[6], int(f[7]), tt[1] if tt else "1"))
    for p in prots:
        ex = sorted(cds[p])
        strand = ex[0][2]
        nt = "".join(seq[s - 1:e] for s, e, *_ in ex)
        if strand == "-":
            nt = str(Seq(nt).reverse_complement())
        first = ex[-1] if strand == "-" else ex[0]
        phase = first[3]
        table = int(ex[0][4])
        aa = str(Seq(nt[phase:]).translate(table=table)).rstrip("*")
        ref = faa[p]
        ident = 100.0 if aa == ref else round(100 * sum(a == b for a, b in zip(aa, ref)) / max(len(aa), len(ref)), 1)
        ds = pyhmmer.easel.TextSequence(name=p.encode(), sequence=ref).digitize(alpha)
        hit = list(pyhmmer.hmmsearch([pf], [ds], E=1e-3))[0]
        box = ";".join(f"{dm.alignment.target_from}-{dm.alignment.target_to}(E={dm.i_evalue:.1e})" for h in hit for dm in h.domains if dm.i_evalue < 1e-3) or "-"
        sc = {}
        for n, h in hmms.items():
            r = list(pyhmmer.hmmsearch([h], [ds], E=1e9, domE=1e9))[0]
            sc[n] = round(max((x.score for x in r), default=0.0), 1)
        exs = ",".join(f"{s}-{e}" for s, e, *_ in ex)
        print(f"{asm}\t{contig}\t{p}\t{strand}\t{ex[0][0]}\t{ex[-1][1]}\t{len(ex)}\t{exs}\t{phase}\t{table}\t{len(ref)}\t{ident}\t{box}\t{sc['sexM']}\t{sc['sexP']}")
