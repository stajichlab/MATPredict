"""List every call in a run whose idiomorph was decided by the classifier on
HSP fragments (classifier_input: hsp_fragment), with species and margin, and
the old label for the same call in a previous run.
usage: fragment_labels.py OLD_RUNS NEW_RUNS"""
import csv, glob, os, sys, yaml
old_dir, new_dir = sys.argv[1:3]
samp = {r["ASMID"]: r for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"))}
def load(d):
    out = {}
    for p in glob.glob(os.path.join(d, "*", "detection_report.yaml")):
        g = os.path.basename(os.path.dirname(p))
        for x in (yaml.safe_load(open(p)) or {}).get("detected") or []:
            out.setdefault((g, x["family"], x["contig"]), x)
    return out
old, new = load(old_dir), load(new_dir)
n = 0
for k, x in sorted(new.items()):
    c = x.get("idiomorph_classifier") or {}
    if c.get("classifier_input") not in ("hsp_fragment", "mixed"):
        continue
    n += 1
    s = samp.get(k[0], {})
    o = old.get(k, {})
    print(f"{k[0]}\t{s.get('SPECIES','')}\t{s.get('STRAIN','')}\t{k[2]}\t{c.get('classifier_input')}\t"
          f"old={o.get('idiomorph','(no call)')}\tnew={x['idiomorph']}\tmargin={c.get('margin')}\t"
          f"scores={c.get('scores')}\tconf={x['confidence']}")
print(f"fragment-decided calls: {n}")
