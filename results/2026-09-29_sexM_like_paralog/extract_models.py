"""Extract modelled protein sequences for every detected Mucoromycota call in the LCG
held-out run (results/2026-09-28_lcg_holdout/runs2), from the genome FASTA and the
report's exon coordinates. Writes models.tsv and models.faa.

Usage: python extract_models.py [N_WORKERS]
"""
import csv, glob, os, sys
from concurrent.futures import ProcessPoolExecutor
import yaml
from Bio import SeqIO
from Bio.Seq import Seq

RUNS = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-28_lcg_holdout/runs2"
GEN = "/bigdata/stajichlab/shared/projects/ZyGoLife/LCG/Annotation/genomes"
KEEP = {"sexM", "sexP", "rnhA", "tptA", "glrA", "algA", "btbA"}


def one(org):
    rep = os.path.join(RUNS, org, "detection_report.yaml")
    try:
        r = yaml.safe_load(open(rep))
    except Exception:
        return []
    calls = r.get("detected") or []
    if not calls:
        return []
    code = r.get("genetic_code") or 1
    need = set()
    for c in calls:
        for g in c.get("gene_evidence") or []:
            if g["gene"] in KEEP and g.get("exons") and str(g.get("status", "")).startswith("polished"):
                need.add(g["contig"])
    fa = os.path.join(GEN, f"{org}.sorted.fasta")
    if not os.path.exists(fa) or not need:
        return []
    seqs = {rec.id: rec.seq for rec in SeqIO.parse(fa, "fasta") if rec.id in need}
    out = []
    for ci, c in enumerate(calls):
        clf = c.get("idiomorph_classifier") or {}
        for g in c.get("gene_evidence") or []:
            if g["gene"] not in KEEP or not g.get("exons"):
                continue
            if not str(g.get("status", "")).startswith("polished"):
                continue
            s = seqs.get(g["contig"])
            if s is None:
                continue
            ex = sorted(g["exons"], key=lambda e: e["start"])
            nt = Seq("".join(str(s[e["start"] - 1:e["end"]]) for e in ex))
            if g["strand"] == "-":
                nt = nt.reverse_complement()
            best = None
            for f in range(3):
                sub = nt[f:]
                sub = sub[: len(sub) // 3 * 3]
                p = str(sub.translate(table=code))
                core = p.rstrip("*")
                score = core.count("*")
                if best is None or score < best[0]:
                    best = (score, core)
            out.append(dict(
                org=org, call=ci, family=c["family"], idiomorph=c["idiomorph"],
                confidence=c["confidence"], locus_class=c.get("locus_class"),
                detection_pass=c.get("detection_pass"), call_contig=c["contig"],
                call_start=c["start"], call_end=c["end"],
                margin=clf.get("margin"), score_minus=(clf.get("scores") or {}).get("Minus"),
                score_plus=(clf.get("scores") or {}).get("Plus"),
                clf_input=clf.get("classifier_input"),
                genes_found=",".join(c.get("genes_found") or []),
                gene=g["gene"], role=g["role"], identity=g["identity"], status=g["status"],
                bitscore=g.get("bitscore"), contig=g["contig"], start=g["start"], end=g["end"],
                strand=g["strand"], n_exons=len(ex), internal_stops=best[0], protein=best[1]))
    return out


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    orgs = sorted(os.listdir(RUNS))
    rows = []
    with ProcessPoolExecutor(n) as ex:
        for res in ex.map(one, orgs, chunksize=4):
            rows.extend(res)
    cols = list(rows[0].keys())
    with open("models.tsv", "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    with open("models.faa", "w") as fo:
        for i, r in enumerate(rows):
            fo.write(f">m{i}|{r['org']}|c{r['call']}|{r['idiomorph']}|{r['gene']}\n{r['protein']}\n")
    print(len(rows), "models from", len({r['org'] for r in rows}), "genomes")
