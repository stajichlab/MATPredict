#!/usr/bin/env python
"""Scan v0.6.0 Ascomycota reports for loci where the partial-coverage rule's TRIGGER holds.

    scan_trigger.py --db DB --reports 'GLOB/*/detection_report.yaml' --out DIR [--workers 4]

Report coverage is null for polished genes, so it is ESTIMATED as (sum of exon lengths / 3) / reference protein
length (the same ratio `pipeline.model_reference_coverage` now writes), capped at 100. Reference lengths come from
the curated db passed in (the current db is a superset of v0.6.0's). The trigger is evaluated as in
`partial_core.evaluate` WITHOUT the contig-edge guard and WITHOUT the profile step (neither is in the report).
Writes: coverage_hist.tsv (non-optional core genes), trigger_grid.tsv (counts per threshold pair),
triggered_PARTIAL80_FULL90.tsv (loci at the chosen thresholds).
"""
import argparse, glob, os, sys
from collections import Counter, defaultdict
from multiprocessing import Pool
from pathlib import Path
import yaml

L = getattr(yaml, "CSafeLoader", yaml.SafeLoader)
FAMILY = "Ascomycota:MAT"
GRID_PARTIAL = (0.6, 0.7, 0.8, 0.9)
GRID_FULL = (0.8, 0.9, 0.95)
MODELLED = ("polished_agree", "polished_disagree", "polished_single")
STATE = {}


def init(db):
    from MATPredict.detect.family_registry import load_all_families
    from MATPredict.detect.pipeline import _curated_protein_lengths, load_record_families
    from MATPredict.detect.classifier import classifier_genes
    db = Path(db)
    fams = [f for f in load_all_families(db) if f.key.phylum == "Ascomycota" and f.key.locus_name == "MAT"]
    rf = load_record_families(db)
    by_record = {}
    _curated_protein_lengths(db, fams, rf, by_record=by_record)
    fam = fams[0]
    optional = {g["name"] for g in fam.genes if g.get("optional")}
    STATE["core"] = {g: i for g, i in classifier_genes(fam).items() if g not in optional}
    STATE["len"] = by_record


def coverage(g):
    ref = STATE["len"].get((g.get("reference_record"), g["gene"]))
    ex = g.get("exons")
    if not ref or not ex:
        return None
    aa = sum(int(e["end"]) - int(e["start"]) + 1 for e in ex) // 3
    return min(100.0, 100.0 * aa / ref)


def scan(path):
    genome = os.path.basename(os.path.dirname(path))
    try:
        doc = yaml.load(open(path), Loader=L)
    except Exception as exc:  # noqa: BLE001
        return genome, [], [("ERR", str(exc))]
    rows, genes = [], []
    for c in doc.get("detected") or []:
        if c.get("family") != FAMILY:
            continue
        stats = []
        for g in c.get("gene_evidence") or []:
            if g["gene"] not in STATE["core"] or g.get("lost_to") or g["status"] not in MODELLED:
                continue
            cov = coverage(g)
            if cov is None:
                continue
            stats.append((g["gene"], STATE["core"][g["gene"]], cov, g.get("bitscore") or 0.0, g.get("identity")))
            genes.append((g["gene"], cov))
        vote = c.get("idiomorph_candidates") or []
        if len(vote) < 2 or not vote[0]["score"] > vote[1]["score"]:
            continue
        w, o = vote[0]["idiomorph"], vote[1]["idiomorph"]
        ws = [s for s in stats if s[1] == w]
        os_ = [s for s in stats if s[1] == o]
        if not ws or not os_:
            continue
        wbest = max(ws, key=lambda s: s[3])
        row = dict(genome=genome, contig=c["contig"], start=c["start"], label=c["idiomorph"], confidence=c["confidence"],
                   locus_class=c.get("locus_class"), winner=w, w_gene=wbest[0], w_cov=round(wbest[2], 1),
                   w_bits=wbest[3], w_ident=wbest[4], w_maxcov=round(max(s[2] for s in ws), 1),
                   other=o, o_gene=max(os_, key=lambda s: s[3])[0], o_maxcov=round(max(s[2] for s in os_), 1),
                   edge_dist=min((sg.get("contig_edge_distance") for sg in c.get("segments") or []
                                  if sg.get("contig_edge_distance") is not None), default=None))
        rows.append(row)
    return genome, rows, genes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True); ap.add_argument("--reports", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    paths = sorted(glob.glob(a.reports))
    print(len(paths), "reports", file=sys.stderr)
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    allrows, hist, errs = [], Counter(), []
    with Pool(a.workers, initializer=init, initargs=(a.db,)) as pool:
        for n, (genome, rows, genes) in enumerate(pool.imap_unordered(scan, paths, chunksize=20)):
            allrows += rows
            errs += [(genome, e[1]) for e in genes if isinstance(e, tuple) and e and e[0] == "ERR"]
            for gene, cov in genes if not (genes and isinstance(genes[0], tuple) and genes[0][0] == "ERR") else []:
                hist[(gene, min(int(cov // 10) * 10, 100))] += 1
            if n % 2000 == 0:
                print(n, file=sys.stderr)
    with open(out / "coverage_hist.tsv", "w") as fh:
        fh.write("gene\tcoverage_bin_start\tn\n")
        for (g, b), n in sorted(hist.items()):
            fh.write(f"{g}\t{b}\t{n}\n")
    with open(out / "trigger_grid.tsv", "w") as fh:
        fh.write("partial_below\tfull_at_least\tloci\tgenomes\tlabel_MAT1-1\tlabel_MAT1-2\n")
        for p in GRID_PARTIAL:
            for f in GRID_FULL:
                sel = [r for r in allrows if r["w_cov"] < p * 100 and r["w_maxcov"] < f * 100 and r["o_maxcov"] >= f * 100]
                lab = Counter(r["label"] for r in sel)
                fh.write(f"{p}\t{f}\t{len(sel)}\t{len({r['genome'] for r in sel})}\t{lab['MAT1-1']}\t{lab['MAT1-2']}\n")
    sel = [r for r in allrows if r["w_cov"] < 80 and r["w_maxcov"] < 90 and r["o_maxcov"] >= 90]
    with open(out / "triggered_PARTIAL80_FULL90.tsv", "w") as fh:
        keys = list(sel[0].keys()) if sel else ["genome"]
        fh.write("\t".join(keys) + "\n")
        for r in sorted(sel, key=lambda r: (r["genome"], r["contig"], r["start"])):
            fh.write("\t".join("" if r[k] is None else str(r[k]) for k in keys) + "\n")
    with open(out / "loci_with_contest.tsv", "w") as fh:
        keys = list(allrows[0].keys())
        fh.write("\t".join(keys) + "\n")
        for r in allrows:
            fh.write("\t".join("" if r[k] is None else str(r[k]) for k in keys) + "\n")
    print(f"{len(allrows)} contested loci; {len(sel)} triggered at 80/90; {len(errs)} read errors", file=sys.stderr)


if __name__ == "__main__":
    main()
