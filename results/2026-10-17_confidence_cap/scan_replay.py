#!/usr/bin/env python
"""Per-locus features from v0.6.0 Ascomycota reports for the confidence-cap replay.

    scan_replay.py --db DB --reports 'GLOB/*/detection_report.yaml' --out loci_features.tsv [--workers 4]

One row per detected locus (every Ascomycota family). Report coverage is null for polished genes, so it is ESTIMATED as in
results/2026-10-17_partial_coverage_vote/scan_trigger.py: min(100, 100 * (sum exon lengths / 3) / reference protein length) with the reference
length from the CURRENT db (a superset of v0.6.0's). Raw (unpolished) genes keep their own search coverage.
The called allele's core genes are the `core_MAT` genes whose `present_in_idiomorphs` is empty or names the call (all idiomorphs for
`undetermined`); `lost_to` genes are skipped. `best_*` = the one with the highest bitscore ("best core evidence").
edge_dist = smallest `contig_edge_distance` over the locus segments (a LOCUS measure; gene-level edge distance is not in the report).
"""
import argparse, glob, os, sys
from multiprocessing import Pool
from pathlib import Path
import yaml

L = getattr(yaml, "CSafeLoader", yaml.SafeLoader)
MODELLED = ("polished_agree", "polished_disagree", "polished_single")
S = {}
COLS = ["genome", "family", "contig", "start", "idiomorph", "confidence", "locus_class", "pass_", "polished", "n_found", "genes_found",
        "n_flank", "n_modelled_core", "best_gene", "best_status", "best_ident", "best_cov", "best_bits", "max_cov", "max_ident",
        "edge_dist", "verification", "classifier", "best_aa", "cands", "n_cand", "margin", "note", "best_gene_edge"]


def init(db):
    from MATPredict.detect.family_registry import load_all_families
    from MATPredict.detect.pipeline import _curated_protein_lengths, load_record_families
    db = Path(db)
    fams = [f for f in load_all_families(db) if f.key.phylum == "Ascomycota"]
    by_record = {}
    _curated_protein_lengths(db, fams, load_record_families(db), by_record=by_record)
    S["len"] = by_record
    S["fam"] = {f"{f.key.phylum}:{f.key.locus_name}": {g["name"]: g for g in f.genes} for f in fams}


def aa_of(g):
    ex = g.get("exons")
    return (sum(int(e["end"]) - int(e["start"]) + 1 for e in ex) // 3) if g["status"] in MODELLED and ex else None


def cov_of(g):
    ex = g.get("exons")
    ref = S["len"].get((g.get("reference_record"), g["gene"]))
    if g["status"] in MODELLED and ex and ref:
        return min(100.0, 100.0 * (sum(int(e["end"]) - int(e["start"]) + 1 for e in ex) // 3) / ref)
    return g.get("coverage")


def scan(path):
    genome = os.path.basename(os.path.dirname(path))
    try:
        doc = yaml.load(open(path), Loader=L)
    except Exception:  # noqa: BLE001
        return []
    rows = []
    for c in doc.get("detected") or []:
        fam = S["fam"].get(c.get("family"))
        if fam is None:
            continue
        idio = c["idiomorph"]
        core, flank = [], 0
        for g in c.get("gene_evidence") or []:
            spec = fam.get(g["gene"], {})
            role = g.get("role") or spec.get("role")
            if role in ("flanking_conserved", "flanking_variable"):
                flank += 1
                continue
            if role != "core_MAT" or g.get("lost_to"):
                continue
            pin = spec.get("present_in_idiomorphs") or []
            if pin and idio != "undetermined" and idio not in pin:
                continue
            core.append(dict(gene=g["gene"], aa=aa_of(g), start=g["start"], end=g["end"], contig=g["contig"], status=g["status"], ident=g.get("identity"), cov=cov_of(g), bits=g.get("bitscore") or 0.0))
        mod = [x for x in core if x["status"] in MODELLED]
        best = max(core, key=lambda x: x["bits"]) if core else None
        covs = [x["cov"] for x in core if x["cov"] is not None]
        ids = [x["ident"] for x in core if x["ident"] is not None]
        edge = min((s.get("contig_edge_distance") for s in c.get("segments") or [] if s.get("contig_edge_distance") is not None), default=None)
        clf = c.get("idiomorph_classifier")
        gedge = None  # per-gene edge distance of the best gene, from its segment's contig_edge_distance and which end is nearer
        if best:
            for sg in c.get("segments") or []:
                ed = sg.get("contig_edge_distance")
                if ed is not None and sg["contig"] == best["contig"] and sg["start"] <= best["start"] and best["end"] <= sg["end"]:
                    gedge = ed + (best["start"] - sg["start"]) if sg["start"] - 1 == ed else ed + (sg["end"] - best["end"])
                    break
        rows.append([genome, c["family"], c["contig"], c["start"], idio, c["confidence"], c.get("locus_class"), c.get("detection_pass"),
                     c.get("polished_genes"), len(c.get("genes_found") or []), "|".join(c.get("genes_found") or []), flank, len(mod),
                     best["gene"] if best else "", best["status"] if best else "", best["ident"] if best else "",
                     round(best["cov"], 1) if best and best["cov"] is not None else "", best["bits"] if best else "",
                     round(max(covs), 1) if covs else "", max(ids) if ids else "", edge if edge is not None else "",
                     c.get("verification") or "", "y" if clf else "", best["aa"] if best and best["aa"] is not None else "",
                     ";".join(f"{x['idiomorph']}={x['score']}" for x in c.get("idiomorph_candidates") or []),
                     len(c.get("idiomorph_candidates") or []), c.get("idiomorph_margin"),
                     (c.get("idiomorph_vote_note") or {}).get("decision", "") if isinstance(c.get("idiomorph_vote_note"), dict) else "", gedge if gedge is not None else ""])
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True); ap.add_argument("--reports", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    paths = sorted(glob.glob(a.reports))
    print(len(paths), "reports", file=sys.stderr)
    n = 0
    with Pool(a.workers, initializer=init, initargs=(a.db,)) as pool, open(a.out, "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for rows in pool.imap_unordered(scan, paths, chunksize=20):
            for r in rows:
                fh.write("\t".join("" if v is None else str(v) for v in r) + "\n"); n += 1
    print(n, "loci", file=sys.stderr)


if __name__ == "__main__":
    main()
