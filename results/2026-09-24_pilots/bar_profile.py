"""Is the modelled-gene bar withholding real loci or background?

For each pilot genome, profile three populations:
  CALLED      the reported loci in genomes that got a call (real-locus profile)
  BG          the withheld loci in those SAME called genomes (background profile:
              the genome's real locus is elsewhere, so these are not it)
  WITHHELD    the best withheld locus in genomes with no call (the unknown)
Features come from the report (genes, modelled count, span) joined to the
evidence-diagnostics rows (admitted to polishing, best identity, roles).
"""
import collections, csv, glob, json, os, statistics as st, sys
import yaml
from yaml import CSafeLoader as L

P = os.path.dirname(os.path.abspath(__file__))
meta = {r["ASMID"]: r for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"))}

def diag_index(path):
    idx = []
    if os.path.exists(path):
        for line in open(path):
            d = json.loads(line)
            if d.get("kind") == "evidence":
                idx.append(d)
    return idx

def features(locus, diags, family):
    rows = [d for d in diags if d["contig"] == locus["contig"] and d["family"] == family
            and not (d["cluster_end"] < locus["start"] or d["cluster_start"] > locus["end"])]
    roles = set().union(*[set(d["roles"]) for d in rows]) if rows else set()
    return {
        "genes": len(locus["genes_found"]),
        "modelled": locus.get("polished_genes", locus.get("polished_genes")),
        "admitted": any(d["admitted"] for d in rows),
        "identity": max([d["best_identity"] for d in rows], default=None),
        "flank": any(r.startswith("flanking") for r in roles),
        "core": "core_MAT" in roles,
        "span": locus["end"] - locus["start"] + 1,
    }

def tier(f):
    if f["genes"] >= 2 and f["core"] and f["flank"]:
        return "T2 core+flank"
    if f["genes"] >= 2:
        return "T1 multi-gene"
    return "T0 single-gene"

def summarize(label, fs):
    if not fs:
        return f"    {label:9s} n=0"
    t = collections.Counter(tier(f) for f in fs)
    ids = [f["identity"] for f in fs if f["identity"] is not None]
    adm = sum(f["admitted"] for f in fs)
    return (f"    {label:9s} n={len(fs):4d}  " + "  ".join(f"{k}:{100*v/len(fs):.0f}%" for k, v in sorted(t.items()))
            + f"  admitted:{100*adm/len(fs):.0f}%  identity median {st.median(ids) if ids else float('nan'):.1f}")

allrows = []
for d in sorted(glob.glob(f"{P}/*/runs")):
    name = d.split("/")[-2]
    called, bg, withheld = [], [], []
    for g in sorted(os.listdir(d)):
        rp = f"{d}/{g}/detection_report.yaml"
        if not os.path.exists(rp) or os.path.getsize(rp) == 0:
            continue
        doc = yaml.load(open(rp), Loader=L) or {}
        diags = diag_index(f"{d}/{g}/evidence_diagnostics.jsonl")
        det = doc.get("detected") or []
        sup = doc.get("suppressed_loci") or []
        if det:
            for x in det:
                called.append(features({**x, "genes_found": x.get("genes_found") or []}, diags, x["family"]))
            for x in sup:
                bg.append(features(x, diags, x["family"]))
        elif sup:
            fs = [(features(x, diags, x["family"]), x) for x in sup]
            best_f, best_x = max(fs, key=lambda p: (p[0]["genes"], p[0]["modelled"] or 0, p[0]["identity"] or 0))
            withheld.append(best_f)
            allrows.append({"panel": name, "genome": g, "order": meta.get(g, {}).get("ORDER"),
                            "genus": meta.get(g, {}).get("GENUS"), "tier": tier(best_f),
                            "family": best_x["family"], "contig": best_x["contig"], "start": best_x["start"],
                            "end": best_x["end"], "idiomorph": best_x["idiomorph"],
                            "genes": ",".join(best_x["genes_found"]), **best_f})
    print(f"## {name}")
    print(summarize("CALLED", called)); print(summarize("BG", bg)); print(summarize("WITHHELD", withheld))

with open(f"{P}/withheld_best_loci.tsv", "w") as fh:
    cols = list(allrows[0].keys())
    fh.write("\t".join(cols) + "\n")
    for r in allrows:
        fh.write("\t".join(str(r[c]) for c in cols) + "\n")
print(f"\nwrote {P}/withheld_best_loci.tsv ({len(allrows)} genomes)")
