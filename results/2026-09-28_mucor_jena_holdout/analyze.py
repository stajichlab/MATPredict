"""Score the Mucor_Jena held-out runs and compare detect's sexM/sexP models to
the funannotate annotation. Keyed by strain folder only; species names from
file names are never used.

Usage: python analyze.py   (reads inputs.tsv, runs_scaffolds/, runs_contigs/,
hmm/*.tbl; writes per_strain.tsv, gene_models.tsv, summary.txt)
"""
import collections, csv, glob, os, re, sys
import yaml
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner, substitution_matrices

O = os.path.dirname(os.path.abspath(__file__))
inputs = {}
for line in open(f"{O}/inputs.tsv"):
    s, sc, ct, pr, gf, ag = (line.rstrip("\n").split("\t") + [""] * 6)[:6]
    inputs[s] = dict(scaffolds=sc, contigs=ct, proteins=pr, gff=gf)
leak = collections.defaultdict(list)
for r in csv.DictReader(open(f"{O}/leakage.tsv"), delimiter="\t"):
    leak[r["strain"]].append(f"{r['source_type']}:{r['match'].split(' | ')[0]}")
# classifier-training leakage found by hand (manifest lists GCA_052058895.1 = CBS 293.63)
TRAIN_LEAK = {"CBS293_63": "BFD copy GCA_052058895.1 sexP protein is in classifier training_extra",
              "CBS210_80": "BFD copy GCA_060309335.1 was one of the 7 cases used to set min_margin 25"}

hmm = {}
for f in glob.glob(f"{O}/hmm/*.tbl"):
    s = os.path.basename(f)[:-4]
    d = collections.defaultdict(dict)
    for l in open(f):
        if l.startswith("#"):
            continue
        x = l.split()
        d[x[0]][x[2]] = float(x[5])
    hmm[s] = d

aligner = PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score, aligner.extend_gap_score = -10, -0.5
aligner.mode = "global"


def ident(a, b):
    if not a or not b:
        return None, None
    aln = aligner.align(a, b)[0]
    m = sum(1 for x, y in zip(*aln[:2]) if x == y and x != "-")
    L = sum(1 for x, y in zip(*aln[:2]) if x != "-" and y != "-")
    return round(100 * m / max(1, L), 1), round(100 * L / max(len(a), len(b)), 1)


def load_gff(path):
    cds = collections.defaultdict(list)
    for l in open(path):
        if l.startswith("#") or "\tCDS\t" not in l:
            continue
        f = l.rstrip("\n").split("\t")
        par = re.search(r"Parent=([^;]+)", f[8]).group(1)
        cds[par].append((f[0], int(f[3]), int(f[4]), f[6]))
    return cds


def model_protein(genome, ev):
    ex = ev.get("exons") or [{"start": ev["start"], "end": ev["end"]}]
    seq = genome.get(ev["contig"])
    if seq is None:
        return ""
    parts = [seq[e["start"] - 1:e["end"]] for e in sorted(ex, key=lambda e: e["start"])]
    nt = Seq("".join(str(p) for p in parts))
    if ev.get("strand") == "-":
        nt = nt.reverse_complement()
    best = ""
    for fr in range(3):
        sub = nt[fr:]
        sub = sub[: len(sub) // 3 * 3]
        aa = str(sub.translate()).split("*")
        cand = max(aa, key=len)
        if len(cand) > len(best):
            best = cand
    return best


def read_report(path):
    if not os.path.exists(path):
        return None
    return yaml.safe_load(open(path))


per, genes = [], []
for s, inp in inputs.items():
    base = dict(strain=s, leakage="; ".join(leak.get(s, [])), training_leak=TRAIN_LEAK.get(s, ""))
    for mode in ("scaffolds", "contigs"):
        r = read_report(f"{O}/runs_{mode}/{s}/detection_report.yaml")
        wall = open(f"{O}/runs_{mode}/{s}/wall_seconds").read().strip() if os.path.exists(f"{O}/runs_{mode}/{s}/wall_seconds") else ""
        if r is None:
            base[f"{mode}_status"] = "no_report" if inp[mode] else "no_input"
            continue
        det = r.get("detected") or []
        base[f"{mode}_status"] = "called" if det else "not_called"
        base[f"{mode}_n_calls"] = len(det)
        base[f"{mode}_routing"] = r.get("routing_mode")
        base[f"{mode}_gate_withheld"] = r.get("suppressed_mat_gene_gate")
        base[f"{mode}_wall_s"] = wall
        calls = []
        for d in det:
            c = d.get("idiomorph_classifier") or {}
            calls.append("{}|{}|{}|{}|{}:{}-{}|clf={}:{}|split={}|verif={}|genes={}".format(
                d["family"], d["idiomorph"], d["confidence"], d["locus_class"], d["contig"], d["start"], d["end"],
                c.get("classifier_input"), c.get("margin"), bool(d.get("split_locus")),
                (d.get("verification") or {}).get("status") if isinstance(d.get("verification"), dict) else d.get("verification"),
                ",".join(d.get("genes_found", []))))
        base[f"{mode}_calls"] = " ;; ".join(calls)
        if mode != "scaffolds" or not det:
            continue
        # gene-model vs annotation, scaffolds run only
        genome = {rec.id: str(rec.seq) for rec in SeqIO.parse(inp["scaffolds"], "fasta")}
        prots = {rec.id: str(rec.seq).rstrip("*") for rec in SeqIO.parse(inp["proteins"], "fasta")}
        cds = load_gff(inp["gff"])
        for d in det:
            for ev in d.get("gene_evidence", []):
                if ev["gene"] not in ("sexM", "sexP"):
                    continue
                mp = model_protein(genome, ev)
                # overlapping annotated mRNA
                best_id, best_ov = "", 0
                for mid, segs in cds.items():
                    if segs[0][0] != ev["contig"]:
                        continue
                    lo, hi = min(x[1] for x in segs), max(x[2] for x in segs)
                    ov = min(hi, ev["end"]) - max(lo, ev["start"])
                    if ov > best_ov:
                        best_id, best_ov = mid, ov
                ap = prots.get(best_id, "")
                idn, cov = ident(mp, ap)
                hs = hmm.get(s, {}).get(best_id, {})
                am, aP = hs.get("sexM", 0.0), hs.get("sexP", 0.0)
                acall = "sexP" if aP - am >= 25 else "sexM" if am - aP >= 25 else ("undetermined" if best_id else "")
                c = d.get("idiomorph_classifier") or {}
                genes.append(dict(
                    strain=s, call_idiomorph=d["idiomorph"], gene=ev["gene"], status=ev.get("status"),
                    gene_matches_call={"sexP": "Plus", "sexM": "Minus"}[ev["gene"]] == d["idiomorph"],
                    contig=ev["contig"], start=ev["start"], end=ev["end"], strand=ev.get("strand"),
                    model_exons=len(ev.get("exons") or []), model_aa=len(mp),
                    annotated=best_id or "NOT_ANNOTATED",
                    annot_exons=len(cds.get(best_id, [])), annot_aa=len(ap),
                    model_vs_annot_identity=idn, model_vs_annot_coverage=cov,
                    detect_classifier_input=c.get("classifier_input"), detect_classifier_call=d["idiomorph"],
                    annot_sexM_bits=am, annot_sexP_bits=aP, annot_classifier_call=acall,
                    agree=("" if not best_id else str(
                        {"sexP": "Plus", "sexM": "Minus"}.get(acall, acall) == d["idiomorph"]))))
    # best annotated classifier hit, independent of detect
    d = hmm.get(s, {})
    if d:
        top = max(d.items(), key=lambda kv: max(kv[1].values()))
        am, aP = top[1].get("sexM", 0.0), top[1].get("sexP", 0.0)
        base["annot_top_protein"] = top[0]
        base["annot_top_sexM_bits"], base["annot_top_sexP_bits"] = am, aP
        base["annot_top_call"] = "sexP" if aP - am >= 25 else "sexM" if am - aP >= 25 else "undetermined"
    per.append(base)

cols = sorted({k for r in per for k in r}, key=lambda k: (k != "strain", k))
with open(f"{O}/per_strain.tsv", "w", newline="") as fo:
    w = csv.DictWriter(fo, fieldnames=["curator_taxonomy", "curator_known_mating_type"] + cols, delimiter="\t")
    w.writeheader()
    for r in per:
        w.writerow({"curator_taxonomy": "", "curator_known_mating_type": "", **r})
if genes:
    with open(f"{O}/gene_models.tsv", "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=list(genes[0]), delimiter="\t")
        w.writeheader()
        w.writerows(genes)

# summary
out = []
for mode in ("scaffolds", "contigs"):
    st = collections.Counter(r.get(f"{mode}_status") for r in per)
    out.append(f"{mode}: {dict(st)}")
    idi = collections.Counter()
    for r in per:
        for c in (r.get(f"{mode}_calls") or "").split(" ;; "):
            if c:
                f = c.split("|")
                idi[(f[1], f[2])] += 1
    out.append(f"  calls by (idiomorph, confidence): {dict(idi)}")
    out.append(f"  genomes with >1 call: {sum(1 for r in per if (r.get(f'{mode}_n_calls') or 0) > 1)}")
ch = [(r["strain"], r.get("scaffolds_status"), r.get("contigs_status")) for r in per
      if r.get("scaffolds_status") != r.get("contigs_status")]
out.append(f"scaffold vs contig status differs: {ch}")
if genes:
    genes_c = [g for g in genes if g["gene_matches_call"]]
    out.append(f"(rows for the called idiomorph's gene only: {len(genes_c)}; annotation shorter than 70% of model: "
               f"{sum(1 for g in genes_c if g['annot_aa'] and g['annot_aa'] < 0.7 * g['model_aa'])})")
    na = sum(1 for g in genes if g["annotated"] == "NOT_ANNOTATED")
    ag = collections.Counter(g["agree"] for g in genes)
    hi = sum(1 for g in genes if g["model_vs_annot_identity"] and g["model_vs_annot_identity"] >= 95)
    out.append(f"sexM/sexP gene models in calls: {len(genes)}; not annotated: {na}; "
               f"model vs annotation identity >=95%: {hi}; classifier on annotated protein agrees with detect: {dict(ag)}")
top = collections.Counter(r.get("annot_top_call") for r in per)
out.append(f"top annotated HMM hit per strain (independent of detect): {dict(top)}")
open(f"{O}/summary.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
