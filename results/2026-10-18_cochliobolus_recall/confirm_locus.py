"""Independent check that a genome carries the curated Bipolaris maydis MAT gene, and what the gene model looks like.

usage: confirm_locus.py GENOME.fna OUTDIR PFAM_SUBSET.hmm
Runs, on the genome: tblastn of the two curated proteins (C5 MAT1-1-1, C4 MAT1-2-1), blastn of each record's locus.gbk,
miniprot (-G 5000 and default) and exonerate protein2genome (maxintron 5000) on the best tblastn window, then translates the
models and scores them with the Pfam alpha-box / HMG-box models (pyhmmer, GA cut-off) and reports the C-terminal/N-terminal
completeness. Everything is written as TSV lines to stdout and files in OUTDIR.
"""
import subprocess, sys, os, re, json
from pathlib import Path

DB = Path(os.environ.get("DBROOT", "db")) / "Ascomycota/Pleosporales"
RECS = {"C5_MAT1-1-1": "5016_c5_MAT_MAT1-1", "C4_MAT1-2-1": "5016_c4_MAT_MAT1-2"}
genome, out, pfam = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
out.mkdir(parents=True, exist_ok=True)

def sh(cmd, **kw):
    return subprocess.run(cmd, shell=True, check=False, capture_output=True, text=True, **kw)

CODON = {}
bases = "TCAG"
aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
i = 0
for a in bases:
    for b in bases:
        for c in bases:
            CODON[a + b + c] = aa[i]; i += 1
def rc(s):
    return s[::-1].translate(str.maketrans("ACGTacgtN", "TGCAtgcaN"))
def tr(s):
    return "".join(CODON.get(s[i:i + 3].upper(), "X") for i in range(0, len(s) - 2, 3))

def read_fasta(p):
    d = {}; k = None
    for l in open(p):
        if l.startswith(">"):
            k = l[1:].split()[0]; d[k] = []
        else:
            d[k].append(l.strip())
    return {k: "".join(v) for k, v in d.items()}

gen = read_fasta(genome)
sh(f"makeblastdb -in {genome} -dbtype nucl -out {out}/gdb")

# record proteins and locus nucleotide sequences
prots = {}
for name, rec in RECS.items():
    f = DB / rec / "proteins.faa"
    seq = "".join(l.strip() for l in open(f) if not l.startswith(">"))
    prots[name] = seq
    (out / f"{name}.faa").write_text(f">{name}\n{seq}\n")
    # locus.gbk -> fasta
    gbk = (DB / rec / "locus.gbk").read_text()
    m = re.search(r"ORIGIN(.*?)//", gbk, re.S)
    nt = re.sub(r"[^acgtnACGTN]", "", m.group(1)) if m else ""
    (out / f"{name}.locus.fna").write_text(f">{name}_locus\n{nt}\n")

rows = []
for name in RECS:
    r = sh(f"tblastn -query {out}/{name}.faa -db {out}/gdb -evalue 1e-5 -max_target_seqs 50 -outfmt '6 sseqid sstart send pident length qlen qstart qend evalue bitscore' ")
    hits = [l.split("\t") for l in r.stdout.strip().splitlines()]
    hits.sort(key=lambda h: -float(h[9]))
    with open(out / f"{name}.tblastn.tsv", "w") as fh:
        for h in hits[:30]:
            fh.write("\t".join(h) + "\n")
    r = sh(f"blastn -query {out}/{name}.locus.fna -db {out}/gdb -evalue 1e-10 -outfmt '6 sseqid sstart send pident length qlen evalue bitscore' ")
    bl = [l.split("\t") for l in r.stdout.strip().splitlines()]
    bl.sort(key=lambda h: -float(h[7]))
    with open(out / f"{name}.blastn.tsv", "w") as fh:
        for h in bl[:20]:
            fh.write("\t".join(h) + "\n")
    top = hits[0] if hits else None
    print(f"[{name}] tblastn top: {top}")
    print(f"[{name}] blastn top: {bl[0] if bl else None}  n_hsp={len(bl)}")
    # miniprot
    for G in ("5000", "200000"):
        r = sh(f"miniprot -u -G {G} --gff -N 3 {genome} {out}/{name}.faa")
        (out / f"{name}.miniprot_G{G}.gff").write_text(r.stdout)
        mrna = [l for l in r.stdout.splitlines() if "\tmRNA\t" in l]
        for l in mrna[:3]:
            f = l.split("\t")
            attrs = dict(kv.split("=") for kv in f[8].split(";") if "=" in kv)
            print(f"[{name}] miniprot G{G}: {f[0]}:{f[3]}-{f[4]}{f[6]} Rank={attrs.get('Rank')} Identity={attrs.get('Identity')} Positive={attrs.get('Positive')} Target={attrs.get('Target')}")
    # exonerate on best window
    if top:
        c, s, e = top[0], int(top[1]), int(top[2])
        lo, hi = min(s, e) - 3000, max(s, e) + 3000
        lo = max(lo, 1); hi = min(hi, len(gen[c]))
        win = gen[c][lo - 1:hi]
        (out / f"{name}.win.fna").write_text(f">{c}:{lo}-{hi}\n{win}\n")
        r = sh(f"exonerate --model protein2genome --refine region --maxintron 5000 --showtargetgff yes --showalignment no --showvulgar yes --bestn 1 {out}/{name}.faa {out}/{name}.win.fna")
        (out / f"{name}.exonerate.txt").write_text(r.stdout)
        vul = [l for l in r.stdout.splitlines() if l.startswith("vulgar")]
        print(f"[{name}] exonerate window {c}:{lo}-{hi}: {vul[0] if vul else 'NO MODEL'}")
        # CDS from exonerate gff
        cds = []
        strand = "+"
        for l in r.stdout.splitlines():
            f = l.split("\t")
            if len(f) > 8 and f[2] == "exon":
                cds.append((int(f[3]), int(f[4]))); strand = f[6]
        cds.sort()
        if cds:
            seq = "".join(win[a - 1:b] for a, b in cds)
            if strand == "-":
                seq = rc(seq)
            prot = tr(seq)
            introns = [cds[i + 1][0] - cds[i][1] - 1 for i in range(len(cds) - 1)]
            ref = prots[name]
            print(f"[{name}] exonerate model: strand {strand} exons {len(cds)} introns {introns} protein {len(prot)} aa (ref {len(ref)}) starts_M={prot.startswith('M')} internal_stops={prot.rstrip('*').count('*')}")
            (out / f"{name}.model.faa").write_text(f">{name}_model_{genome.stem}\n{prot.rstrip('*')}\n")
            # global-ish identity by simple alignment via parasail-free approach: use blast pairing
            rp = sh(f"printf '>q\\n{ref}\\n' > {out}/_q.faa; printf '>s\\n{prot.rstrip('*')}\\n' > {out}/_s.faa; blastp -query {out}/_q.faa -subject {out}/_s.faa -outfmt '6 pident length qstart qend sstart send' -max_hsps 1")
            print(f"[{name}] model vs record blastp: {rp.stdout.strip()}")

# Pfam on the models
try:
    import pyhmmer
    from pyhmmer import easel, plan7
    with plan7.HMMFile(pfam) as hf:
        hmms = list(hf)
    for name in RECS:
        mp = out / f"{name}.model.faa"
        if not mp.exists():
            continue
        with easel.SequenceFile(str(mp), digital=True, alphabet=easel.Alphabet.amino()) as sf:
            seqs = sf.read_block()
        for q in (seqs,):
            for hit_list in pyhmmer.hmmsearch(hmms, seqs, cpus=1, bit_cutoffs="gathering"):
                for h in hit_list:
                    for d in h.domains.included:
                        print(f"[{name}] Pfam {hit_list.query.name.decode() if isinstance(hit_list.query.name, bytes) else hit_list.query.name} {hit_list.query.accession} score={d.score:.1f} env={d.env_from}-{d.env_to}")
except Exception as ex:
    print("pyhmmer step failed:", ex)
