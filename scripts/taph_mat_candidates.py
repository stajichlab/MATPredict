#!/usr/bin/env python3
"""Collect every candidate MAT protein of a set of genomes, low confidence included, for a gene tree.

For each genome the reference MAT proteins of the database are aligned to the genome with miniprot (spliced protein-to-genome alignment); every alignment that passes --min-identity and
--min-coverage is kept, overlapping alignments of one locus are merged to the best, the CDS is read from the genome and translated. The reference sets are
  PM        mat1-P/M and matPi/Mc/Mi proteins of the Schizosaccharomycetales, Pneumocystidales and Taphrinales records (the Taphrinomycotina P and M genes)
  alpha_a   alpha-box and HMG-box MAT proteins of the other Ascomycota records, one record per order for each gene name (MAT1-1-1, MAT1-2-1, MATalpha1, MATa1, MATA1 ...)
  plus_minus  sexP / sexM proteins of the Mucoromycota records (HMG, Plus / Minus)
Each candidate says whether it lies in a locus the detection reports of the genome list (detected, or withheld as low confidence). Hits of non-MAT HMG-box and alpha-box genes are expected: the
best reference name and its identity are in the table, so they can be dropped or kept as outgroups; nothing is filtered by confidence.

    taph_mat_candidates.py --db db --genomes GENOMES.tsv --genome-dir DIR --out OUTDIR [--reports-glob 'RUNS/*/GENOME/detection_report.yaml'] [--min-identity 0.25] [--min-coverage 0.4]
GENOMES.tsv: a column `genome` (file DIR/<genome>.fa.gz) and `species`. Needs miniprot, Biopython and PyYAML. Output: candidates.tsv, candidates.faa, references.faa.
"""
import argparse, collections, csv, glob, gzip, os, re, subprocess, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from idiomorph_labels import gene_label, idiomorph_label, idiomorph_of_gene, lineage_of

PM_ORDERS = ("Schizosaccharomycetales", "Pneumocystidales", "Taphrinales")
ALPHA_A = re.compile(r"^(MAT1-1-1|MAT1-2-1|MATALPHA1|MATA1|MATALPHA2|MATA2|MTLALPHA1|MTLA1|MTLALPHA2|MTLA2)$", re.I)
PLUS_MINUS = re.compile(r"^(sexP|sexM)$", re.I)


def read_fasta(path):
    name, seq = None, []
    for ln in open(path):
        if ln.startswith(">"):
            if name: yield name, "".join(seq)
            name, seq = ln[1:].strip(), []
        else: seq.append(ln.strip())
    if name: yield name, "".join(seq)


def label_cols(gene, detected):
    """(gene label, idiomorph label, lineage, basis) of a hit on reference gene `gene`; `detected` = the hit lies in a locus the pipeline reported."""
    basis = "pipeline_detected" if detected else "similarity_only"; m = idiomorph_of_gene(gene)
    return gene_label(gene, basis), idiomorph_label(m, ["P", "M"]) if m else "", lineage_of(["P", "M"]) if m else "", basis


def gene_of(header):
    m = re.search(r"name=([^|\s]+)", header); return m.group(1) if m else header.split("|")[0]


def references(db):
    """[(id, class, gene, protein)]: PM = every protein of the Taphrinomycotina records; alpha_a = one protein per Ascomycota order and gene name; plus_minus = up to 6 per sexP / sexM."""
    refs, per = [], collections.Counter()
    for f in sorted(glob.glob(os.path.join(db, "*", "*", "*", "proteins.faa"))):
        parts = f.split(os.sep); top, order, rec = parts[-4], parts[-3], parts[-2]
        for h, sq in read_fasta(f):
            g = gene_of(h); sq = sq.replace("*", "")
            if len(sq) < 40: continue
            if order in PM_ORDERS: refs.append((f"{rec}|{g}", "PM", g, sq))
            elif top == "Mucoromycota" and PLUS_MINUS.match(g):
                if per[("mu", g.upper())] < 6: per[("mu", g.upper())] += 1; refs.append((f"{rec}|{g}", "plus_minus", g, sq))
            elif top == "Ascomycota" and ALPHA_A.match(g):
                if per[(order, g.upper())] < 1: per[(order, g.upper())] += 1; refs.append((f"{rec}|{g}", "alpha_a", g, sq))
    return refs


def load_genome(path):
    op = gzip.open if path.endswith(".gz") else open; seqs, cur, name = {}, [], None
    with op(path, "rt") as fh:
        for ln in fh:
            if ln.startswith(">"):
                if name: seqs[name] = "".join(cur).upper()
                name, cur = ln[1:].split()[0], []
            else: cur.append(ln.strip())
    if name: seqs[name] = "".join(cur).upper()
    return seqs


def parse_gff(text):
    mr = {}
    for ln in text.splitlines():
        if ln.startswith("#") or not ln.strip(): continue
        p = ln.split("\t")
        if len(p) < 9: continue
        at = dict(kv.split("=", 1) for kv in p[8].split(";") if "=" in kv)
        if p[2] == "mRNA": mr[at["ID"]] = dict(contig=p[0], start=int(p[3]), end=int(p[4]), strand=p[6], score=float(p[5]) if p[5] != "." else 0.0, ident=float(at.get("Identity", 0)), pos=float(at.get("Positive", 0)), target=at["Target"].split()[0], tstart=int(at["Target"].split()[1]), tend=int(at["Target"].split()[2]), cds=[])
        elif p[2] == "CDS" and at.get("Parent") in mr: mr[at["Parent"]]["cds"].append((int(p[3]), int(p[4]), int(p[7]) if p[7].isdigit() else 0))
    return list(mr.values())


def revcomp(s): return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def translate(seqs, hit):
    from Bio.Seq import Seq
    sq = seqs[hit["contig"]]; cds = sorted(hit["cds"])
    dna = "".join(sq[a - 1:b] for a, b, _ in cds)
    if hit["strand"] == "-": dna = revcomp(dna)
    dna = dna[: len(dna) // 3 * 3]
    return str(Seq(dna).translate(table=1))


def report_loci(glob_pat, genome):
    """[(contig, start, end, idiomorph, status)] of the detection report of the genome: detected loci and loci withheld as low confidence."""
    import yaml
    L = getattr(yaml, "CSafeLoader", yaml.SafeLoader); out = []
    for f in glob.glob(glob_pat.replace("GENOME", genome)):
        d = yaml.load(open(f), Loader=L) or {}
        for l in d.get("detected") or []: out.append((l.get("contig"), int(l["start"]), int(l["end"]), l.get("idiomorph") or "", "detected"))
        for l in d.get("suppressed_loci") or []: out.append((l.get("contig"), int(l["start"]), int(l["end"]), l.get("idiomorph") or "", "withheld"))
    return out


def domain_pass(a, rows, refs_path, rcls, out):
    """tblastn of the reference proteins against the report locus windows of each genome; one peptide per non-overlapping hit (the best bit score wins)."""
    tsv = open(os.path.join(out, "domain_hits.tsv"), "w"); faa = open(os.path.join(out, "domain_hits.faa"), "w")
    tsv.write("\t".join(["hit", "genome", "species", "contig", "start", "end", "best_ref", "ref_class", "ref_gene", "ref_gene_label", "idiomorph_label", "lineage", "label_basis", "identity", "aligned_aa", "bitscore", "evalue", "window_status"]) + "\n"); n = 0
    for r in rows:
        g = r["genome"]; path = os.path.join(a.genome_dir, g + ".fa.gz"); rl = report_loci(a.reports_glob, g) if a.reports_glob else []
        if not os.path.exists(path) or not rl: continue
        seqs = load_genome(path); tmp = tempfile.mkdtemp(dir=os.environ.get("SCRATCH", "/tmp")); wins = {}
        with open(os.path.join(tmp, "w.fa"), "w") as fh:
            for i, (c, s0, e0, idio, st) in enumerate(rl):
                if c not in seqs: continue
                wid = f"w{i}"; wins[wid] = (c, s0, e0, st, idio); fh.write(f">{wid}\n{seqs[c][s0 - 1:e0]}\n")
        if not wins: continue
        p = subprocess.run([a.tblastn, "-query", refs_path, "-subject", os.path.join(tmp, "w.fa"), "-evalue", "1e-3", "-seg", "no", "-outfmt", "6 qseqid sseqid pident length sstart send evalue bitscore sseq"], capture_output=True, text=True)
        hits = []
        for ln in p.stdout.splitlines():
            q, sid, pid, ln_, ss, se, ev, bs, sseq = ln.split("\t")
            if float(bs) >= a.min_bits: hits.append((float(bs), q, sid, float(pid), int(ln_), int(ss), int(se), ev, sseq))
        hits.sort(reverse=True); kept = []
        for bs, q, sid, pid, ln_, ss, se, ev, sseq in hits:
            lo, hi = sorted((ss, se))
            if any(k[0] == sid and min(hi, k[2]) - max(lo, k[1]) > 0.5 * min(hi - lo, k[2] - k[1]) for k in kept): continue
            kept.append((sid, lo, hi)); c, s0, e0, st, idio = wins[sid]; cls, gene = rcls[q]; n += 1
            hid = f"{g}|{c}:{s0 + lo - 1}-{s0 + hi - 1}"; pep = sseq.replace("-", "").replace("*", "")
            tsv.write("\t".join(map(str, [hid, g, r.get("species", ""), c, s0 + lo - 1, s0 + hi - 1, q, cls, gene, *label_cols(gene, st == "detected"), pid, len(pep), bs, ev, st])) + "\n")
            faa.write(f">{hid} {r.get('species', '')} best={label_cols(gene, st == 'detected')[0]}({cls}) id={pid:.0f} bits={bs:.0f} window={st}\n{pep}\n")
        subprocess.run(["rm", "-rf", tmp])
    tsv.close(); faa.close(); print(n, "domain-level hits")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0]); ap.add_argument("--db", required=True); ap.add_argument("--genomes", required=True); ap.add_argument("--genome-dir", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--reports-glob"); ap.add_argument("--min-identity", type=float, default=0.25); ap.add_argument("--min-coverage", type=float, default=0.4)
    ap.add_argument("--miniprot", default="miniprot"); ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--domain-pass", action="store_true", help="also tblastn the references against every locus window of the detection reports (detected and withheld) and keep each hit of at least --min-bits as an aligned peptide (domain level: what the pipeline saw before it tried to model a gene)")
    ap.add_argument("--tblastn", default="tblastn"); ap.add_argument("--min-bits", type=float, default=33.0)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    refs = references(a.db); rlen = {f"{r[0]}|{r[1]}": len(r[3]) for r in refs}; rcls = {f"{r[0]}|{r[1]}": (r[1], r[2]) for r in refs}      # the FASTA header (miniprot's target name) is id|class
    with open(os.path.join(a.out, "references.faa"), "w") as fh:
        for rid, c, g, s in refs: fh.write(f">{rid}|{c}\n{s}\n")
    print(len(refs), "reference proteins:", dict(collections.Counter(r[1] for r in refs)), flush=True)
    rows = list(csv.DictReader(open(a.genomes), delimiter="\t")); tsv = open(os.path.join(a.out, "candidates.tsv"), "w"); faa = open(os.path.join(a.out, "candidates.faa"), "w")
    cols = ["candidate", "genome", "species", "contig", "start", "end", "strand", "best_ref", "ref_class", "ref_gene", "ref_gene_label", "idiomorph_label", "lineage", "label_basis", "identity", "positive", "query_coverage", "protein_aa", "has_stop", "report_locus"]
    tsv.write("\t".join(cols) + "\n"); n_all = 0
    for r in rows:
        g = r["genome"]; path = os.path.join(a.genome_dir, g + ".fa.gz")
        if not os.path.exists(path): continue
        seqs = load_genome(path)
        p = subprocess.run([a.miniprot, "-t", str(a.threads), "--gff", "-N", "10", "--outs", "0.3", "--outc", "0.3", path, os.path.join(a.out, "references.faa")], capture_output=True, text=True)
        hits = []
        for h in parse_gff(p.stdout):
            ql = rlen.get(h["target"], 0) or 1; cov = (h["tend"] - h["tstart"] + 1) / ql
            if h["ident"] >= a.min_identity and cov >= a.min_coverage and h["cds"] and h["contig"] in seqs: h["cov"] = cov; hits.append(h)
        hits.sort(key=lambda h: -(h["ident"] * h["cov"]))
        loci = []      # merge hits that overlap an already kept hit by more than half of the shorter one, same contig
        for h in hits:
            for k in loci:
                if k["contig"] == h["contig"] and min(k["end"], h["end"]) - max(k["start"], h["start"]) > 0.5 * min(k["end"] - k["start"], h["end"] - h["start"]): break
            else: loci.append(h)
        rl = report_loci(a.reports_glob, g) if a.reports_glob else []
        for i, h in enumerate(sorted(loci, key=lambda h: (h["contig"], h["start"])), 1):
            prot = translate(seqs, h); stop = prot.endswith("*"); prot = prot.rstrip("*")
            ov = [f"{st}:{idio}" for c, s0, e0, idio, st in rl if c == h["contig"] and min(e0, h["end"]) - max(s0, h["start"]) > 0]
            cid = f"{g}|{h['contig']}:{h['start']}-{h['end']}({h['strand']})"; cls, gene = rcls[h["target"]]
            tsv.write("\t".join(map(str, [cid, g, r.get("species", ""), h["contig"], h["start"], h["end"], h["strand"], h["target"], cls, gene, *label_cols(gene, any(x.startswith("detected") for x in ov)), round(h["ident"], 3), round(h["pos"], 3), round(h["cov"], 3), len(prot), int(stop), ",".join(sorted(set(ov))) or "none"])) + "\n")
            faa.write(f">{cid} {r.get('species', '')} best={label_cols(gene, any(x.startswith('detected') for x in ov))[0]}({cls}) id={h['ident']:.2f} cov={h['cov']:.2f}\n{prot}\n"); n_all += 1
        print(g, len(loci), "candidate loci", flush=True)
    tsv.close(); faa.close(); print(n_all, "candidates")
    if a.domain_pass: domain_pass(a, rows, os.path.join(a.out, "references.faa"), rcls, a.out)


if __name__ == "__main__":
    main()
