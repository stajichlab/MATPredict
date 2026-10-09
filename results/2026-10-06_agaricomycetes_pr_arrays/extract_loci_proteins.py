#!/usr/bin/env python3
"""Gene models and proteins for every STE3-like locus of one genome.

Re-derives the loci exactly as array_scan.py / scan_genome.py (miniprot of ste3_all.faa, -I --outn=100,
qcov >= 0.5, <= 8 kb, merged per strand), then runs miniprot --gff with the best query of each locus
(same parameters) and, per locus, takes the model of that query overlapping the locus (most overlap, then
identity), builds the CDS (miniprot CDS features include the stop codon) from the genome and translates it (genetic code 1).

Usage: extract_loci_proteins.py ASMID OUTDIR
Writes OUTDIR/<ASMID>.prot.tsv and OUTDIR/<ASMID>.prot.faa
locus start/end: 0-based half-open (PAF), as in loci_all.tsv.gz; model_start/model_end: 1-based GFF.
"""
import importlib.util, os, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("scan_genome", os.path.join(HERE, "scan_genome.py"))
sg = importlib.util.module_from_spec(spec); spec.loader.exec_module(sg)
sg.STE3 = os.path.join(HERE, "ste3_all.faa")
THREADS = os.environ.get("SLURM_CPUS_PER_TASK", "8")
COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
COLS = ["locus_id", "genome", "contig", "start", "end", "strand", "best_query", "best_qcov", "locus_ident",
        "model_query", "model_start", "model_end", "model_identity", "model_positive", "model_qcov", "n_cds",
        "prot_len", "complete", "partial_reason", "start_met", "has_stop", "n_internal_stop", "frameshift"]


def parse_gff(path):
    models = {}
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        ctg, typ, s, e, strand, phase, attr = f[0], f[2], int(f[3]), int(f[4]), f[6], f[7], f[8]
        a = dict(x.split("=", 1) for x in attr.split(";") if "=" in x)
        if typ == "mRNA":
            t = a.get("Target", "").split()
            models[a["ID"]] = dict(contig=ctg, start=s, end=e, strand=strand, ident=float(a.get("Identity", 0)),
                                   pos=float(a.get("Positive", 0)), q=t[0] if t else "",
                                   qs=int(t[1]) if len(t) > 2 else 0, qe=int(t[2]) if len(t) > 2 else 0,
                                   cds=[], stop=None, attr=a)
        elif typ == "CDS" and a.get("Parent") in models:
            models[a["Parent"]]["cds"].append((s, e, int(phase) if phase.isdigit() else 0))
        elif typ == "stop_codon" and a.get("Parent") in models:
            models[a["Parent"]]["stop"] = (s, e)
    return models


def rc(s):
    return s.translate(COMP)[::-1]


def build(m, seq):
    """CDS nt in transcript orientation (miniprot CDS features include the stop codon), frameshift flag."""
    segs = sorted(m["cds"])
    if m["strand"] == "-":
        tsegs = [(s, e, p) for s, e, p in reversed(segs)]
        cds = "".join(rc(seq[s - 1:e]) for s, e, _ in tsegs)
    else:
        tsegs = segs
        cds = "".join(seq[s - 1:e] for s, e, _ in tsegs)
    cum, fs = 0, False   # segment phase must be (3 - cumulative length % 3) % 3, else a frameshift
    for s, e, p in tsegs:
        fs |= p != (3 - cum % 3) % 3
        cum += e - s + 1
    return cds, fs or cum % 3 != 0 or "Frameshift" in m["attr"]


def main():
    asm, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)
    tmp = tempfile.mkdtemp(prefix="prot_" + asm + "_", dir=os.environ.get("SCRATCH", "/tmp"))
    gz = os.path.join(sg.LIB, asm + ".fa.gz")
    loci = sg.miniprot_loci(gz, tmp, threads=int(THREADS))
    q_all = sg.read_fasta(sg.STE3)
    qlen = {k: len(v) for k, v in q_all.items()}
    tsv, faa = open(os.path.join(outdir, asm + ".prot.tsv"), "w"), open(os.path.join(outdir, asm + ".prot.faa"), "w")
    tsv.write("\t".join(COLS) + "\n")
    models, ncomp = {}, 0
    if loci:
        qf = os.path.join(tmp, "bq.faa")
        with open(qf, "w") as fo:
            for q in sorted({L["best_query"] for L in loci}):
                fo.write(f">{q}\n{q_all[q]}\n")
        gff = os.path.join(tmp, "m.gff")
        with open(gff, "w") as fo:
            subprocess.run([f"{sg.BIN}/miniprot", "-t", THREADS, "-I", "--outn=100", "--gff", gz, qf], stdout=fo,
                           stderr=subprocess.DEVNULL, check=True)
        models = parse_gff(gff)
    by = {}
    for m in models.values():
        by.setdefault((m["contig"], m["strand"]), []).append(m)
    seqs = sg.read_fasta(gz) if loci else {}
    for L in sorted(loci, key=lambda x: (x["contig"], x["start"])):
        ctg, st, en, strand = L["contig"], L["start"], L["end"], L["strand"]
        lid = f"{asm}|{ctg}:{st}-{en}{strand}"
        cand = []
        for m in by.get((ctg, strand), []):
            ov = min(en, m["end"]) - max(st, m["start"] - 1)
            if ov > 0 and m["cds"]:
                cand.append((m["q"] == L["best_query"], ov, m["ident"], m))
        row = dict(locus_id=lid, genome=asm, contig=ctg, start=st, end=en, strand=strand, best_query=L["best_query"],
                   best_qcov=f"{L['best_qcov']:.2f}", locus_ident=f"{L['best_ident']:.3f}")
        if not cand:
            row.update(partial_reason="no_model", complete="False")
            tsv.write("\t".join(str(row.get(c, "")) for c in COLS) + "\n")
            continue
        m = max(cand, key=lambda x: x[:3])[3]
        cds, fs = build(m, seqs[ctg])
        prot = sg.translate(cds)
        has_stop = prot.endswith("*")
        if has_stop:
            prot = prot[:-1]
        nint = prot.count("*")
        met = prot.startswith("M")
        reasons = [r for r, bad in (("no_start", not met), ("no_stop", not has_stop), ("internal_stop", nint > 0),
                                    ("frameshift", fs)) if bad]
        complete = not reasons
        ncomp += complete
        row.update(model_query=m["q"], model_start=m["start"], model_end=m["end"], model_identity=f"{m['ident']:.3f}",
                   model_positive=f"{m['pos']:.3f}",
                   model_qcov=f"{(m['qe'] - m['qs'] + 1) / max(1, qlen.get(m['q'], 1)):.2f}", n_cds=len(m["cds"]),
                   prot_len=len(prot), complete=complete, partial_reason=",".join(reasons), start_met=met,
                   has_stop=has_stop, n_internal_stop=nint, frameshift=fs)
        tsv.write("\t".join(str(row.get(c, "")) for c in COLS) + "\n")
        p = prot.replace("*", "X")
        faa.write(f">{lid} q={m['q']} complete={complete}\n" + "\n".join(p[i:i + 80] for i in range(0, len(p), 80)) + "\n")
    tsv.close(); faa.close()
    subprocess.run(["rm", "-r", tmp])
    print(asm, len(loci), "loci", ncomp, "complete", flush=True)


main()
