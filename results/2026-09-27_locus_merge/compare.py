"""Locus-merge + CAAX-unverified measurement vs the CAAX 'on' run.

Old: results/2026-09-27_caax_precursor/on_{agaricales_panel,polyrussu_controls}
     (code c8412b7 + basidio-anchors db, run-f017ecb)
New: results/2026-09-27_locus_merge/{agaricales_panel,polyrussu_controls}
     (locus-merge 4bf4672 + basidio-anchors db, run-merge-measure)
"""
import collections, csv, glob, os, yaml

BASE = "/bigdata/stajichlab/jstajich/projects/MATPredict/results"
samp = {r["ASMID"]: r for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"))}


def load(d):
    out = {}
    for f in glob.glob(f"{d}/runs/*/detection_report.yaml"):
        g = f.split("/")[-2]
        try:
            out[g] = yaml.safe_load(open(f)) or {}
        except Exception:
            pass
    return out


for panel in ("agaricales_panel", "polyrussu_controls"):
    old = load(f"{BASE}/2026-09-27_caax_precursor/on_{panel}")
    new = load(f"{BASE}/2026-09-27_locus_merge/{panel}")
    both = sorted(set(old) & set(new))
    print(f"== {panel}: genomes old {len(old)} new {len(new)} both {len(both)}")
    fam_old, fam_new = collections.Counter(), collections.Counter()
    merges = collections.Counter(); unver = collections.Counter(); unver_calls = []
    lost_locus = []
    for g in both:
        for x in old[g].get("detected") or []:
            fam_old[x["family"]] += 1
        for x in new[g].get("detected") or []:
            fam_new[x["family"]] += 1
            if x.get("merged_from"):
                merges["+".join(sorted(m["family"].split(":")[1] for m in x["merged_from"]))] += 1
            v = x.get("verification") or {}
            if v.get("status") == "unverified" and v.get("taxon"):
                unver[v["taxon"]] += 1
                unver_calls.append((g, samp.get(g, {}).get("FAMILY"), x["family"], x["confidence"]))
        # every old locus position must still be covered by some new call (merge only folds)
        new_spans = [(x["contig"], x["start"], x["end"]) for x in new[g].get("detected") or []]
        for x in old[g].get("detected") or []:
            if not any(c == x["contig"] and s <= x["end"] and e >= x["start"] for c, s, e in new_spans):
                lost_locus.append((g, x["family"], x["contig"], x["start"]))
    print(" calls old", sum(fam_old.values()), "new", sum(fam_new.values()))
    for k in sorted(set(fam_old) | set(fam_new)):
        print(f"   {k}: {fam_old[k]} -> {fam_new[k]}")
    print(" merges", dict(merges))
    print(" CAAX-unverified by taxon", dict(unver), "total", sum(unver.values()))
    print(" old loci with no covering new call:", len(lost_locus))
    for l in lost_locus[:10]:
        print("   ", l)
    for g in both:
        sp = samp.get(g, {}).get("SPECIES", "")
        if any(k in sp for k in ("Schizophyllum", "Coprinopsis", "Russula nobilis")):
            print(f" {sp} {g}")
            for x in new[g].get("detected") or []:
                print("    ", x["family"], x["contig"], x["start"], x["end"], x["confidence"],
                      [m["family"].split(":")[1] for m in x.get("merged_from") or []],
                      (x.get("verification") or {}).get("status"))
    rs = f"{BASE}/2026-09-27_locus_merge/{panel}/rollout_summary.yaml"
    if os.path.exists(rs):
        s = yaml.safe_load(open(rs))
        uv = s.get("unverified_calls") or {}
        print(" rollout unverified_calls: genomes", len(uv), "calls", sum(uv.values()))
