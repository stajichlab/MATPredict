"""One line per run: calls on scaffold_145 / its contigs, and the not_detected reason."""
import glob, os, sys, yaml
LOCUS = {"scaffold_145", "contig_551", "contig_552", "contig_553"}
for p in sorted(glob.glob(os.path.join(sys.argv[1] if len(sys.argv) > 1 else "runs", "*", "*", "detection_report.yaml"))):
    sha, arm = p.split(os.sep)[-3:-1]
    r = yaml.safe_load(open(p)) or {}
    det = [d for d in r.get("detected") or [] if d["contig"] in LOCUS]
    sup = [d for d in r.get("suppressed_loci") or [] if d["contig"] in LOCUS]
    out = [f"{d['contig']}:{d['start']}-{d['end']} {d['idiomorph']}/{d['confidence']}/{d['locus_class']} {d.get('genes_found')}" for d in det]
    out += [f"WITHHELD {d['contig']}:{d['start']}-{d['end']} {d.get('idiomorph')} polished={d.get('polished_genes')} {d.get('genes_found')}" for d in sup]
    nd = [f"ND {n.get('reason','')[:90]} found={n.get('genes_found')}" for n in r.get("not_detected") or []]
    print(f"{sha} {arm:9s} " + (" | ".join(out) if out else " | ".join(nd) or "nothing"))
