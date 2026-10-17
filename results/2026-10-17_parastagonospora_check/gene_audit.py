#!/usr/bin/env python
"""Score every MAT-named gene entry of the Parastagonospora Ascomycota:MAT loci on both core profiles.

    gene_audit.py --profiles DIR --genome-list FILE --genomes DIR --run v060 DIR --run records DIR --out audit.tsv

One row per gene_evidence entry whose gene name starts with MAT (lost_to entries included and flagged).
Needs MATPredict (PYTHONPATH=src) and score_core_profiles.py (branch ascomycota-core-hmms) on the path.
"""
import argparse, csv, glob, os, sys
from collections import defaultdict
from pathlib import Path
import yaml
sys.path.insert(0, str(Path(__file__).resolve().parent))
import score_core_profiles as sp

L = getattr(yaml, "CSafeLoader", yaml.SafeLoader)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profiles", type=Path, required=True); ap.add_argument("--genome-list", required=True)
    ap.add_argument("--genomes", required=True); ap.add_argument("--run", nargs=2, action="append", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    wanted = {l.strip() for l in open(a.genome_list) if l.strip()}
    hmms = sp.load_hmms(a.profiles)
    rows = []
    for name, root in a.run:
        for g in sorted(wanted):
            fs = glob.glob(f"{root}/**/{g}/detection_report.yaml", recursive=True)
            if not fs: continue
            doc = yaml.load(open(fs[0]), Loader=L)
            loci = [c for c in doc.get("detected") or [] if c.get("family") == "Ascomycota:MAT"]
            if not loci: continue
            seqs = sp.genome_contigs(a.genomes, g, {c["contig"] for c in loci})
            for li, c in enumerate(loci):
                for e in c.get("gene_evidence") or []:
                    if not e["gene"].startswith("MAT") or e["contig"] not in seqs: continue
                    prots = sp.proteins(e, seqs[e["contig"]])
                    sc = sp.best_scores(hmms, prots)
                    mod = e["status"] in ("polished_agree", "polished_single", "polished_disagree")
                    aa = len(prots[0]) if mod else (e["end"] - e["start"] + 1) // 3
                    rows.append(dict(run=name, genome=g, locus=li, contig=c["contig"], lstart=c["start"], label=c["idiomorph"],
                        conf=c.get("confidence"), lclass=c.get("locus_class"), gene=e["gene"], ref=e.get("reference_record"),
                        method=e.get("method"), status=e["status"], modelled=mod, lost_to=e.get("lost_to"),
                        identity=e.get("identity"), bitscore=e.get("bitscore"), evalue=e.get("evalue"), aa=aa,
                        alpha=round(sc["MAT1-1-1"], 1), hmg=round(sc["MAT1-2-1"], 1), strand=e["strand"],
                        gstart=e["start"], gend=e["end"]))
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
    print(len(rows), "entries")
main()
