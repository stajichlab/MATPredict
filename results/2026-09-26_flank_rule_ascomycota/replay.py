"""Replay the flank-carried rule (src/MATPredict/detect/flank_carried.py at 3684c60)
on existing detection reports and write one row per changed call.
usage: replay.py OUT.tsv RUNS_DIR [RUNS_DIR ...]"""
import csv, glob, json, os, sys, yaml
from yaml import CSafeLoader as L

PAD = 3000
MODELLED = {"polished_agree", "polished_disagree", "polished_single"}
meta = {r["ASMID"]: r for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"))}
supp = {l.split(",")[0].strip() for l in open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/data/curation/suppress.txt") if l.strip() and not l.startswith(("#", "ASMID"))}

def modelled(e):
    return e.get("status") in MODELLED or e.get("method") == "diamond_proteome"

def fmt(e):
    return f'{e["gene"]}:{e["start"]}-{e["end"]}:{e.get("identity")}:{e.get("coverage")}:{e.get("status")}'

out = csv.writer(open(sys.argv[1], "w"), delimiter="\t")
out.writerow(["panel", "genome", "phylum", "class", "order", "family_tax", "species", "mat_family", "contig", "start", "end",
              "confidence", "locus_class", "idiomorph", "outcome", "why", "flank_lo", "flank_hi", "n_flank_genes",
              "flank_genes", "core_genes", "max_core_dist_bp", "polish_capped_cluster", "core_detail", "flank_detail"])
tot = {}
for runs in sys.argv[2:]:
    panel = runs.rstrip("/").split("/")[-2]
    for p in sorted(glob.glob(f"{runs}/*/detection_report.yaml")):
        g = p.split("/")[-2]
        if g in supp: continue
        det = (yaml.load(open(p), Loader=L) or {}).get("detected") or []
        caps = {}
        dp = os.path.join(os.path.dirname(p), "evidence_diagnostics.jsonl")
        if os.path.exists(dp):
            for l in open(dp):
                r = json.loads(l)
                if r.get("kind") == "evidence" and r.get("admitted"):
                    caps.setdefault((r["family"], r["contig"]), []).append((r["cluster_start"], r["cluster_end"], r.get("polish_capped")))
        for x in det:
            tot[panel] = tot.get(panel, 0) + 1
            ev = x.get("gene_evidence") or []
            core = [e for e in ev if e["role"] == "core_MAT"]
            fl = [e for e in ev if e["role"].startswith("flanking")]
            if any(modelled(e) for e in core) or not any(modelled(e) for e in fl): continue
            c = fl[0]["contig"]
            lo = min(e["start"] for e in fl if e["contig"] == c) - PAD
            hi = max(e["end"] for e in fl if e["contig"] == c) + PAD
            inside = bool(core) and all(e["contig"] == c and lo <= e["start"] and e["end"] <= hi for e in core)
            if not core: why = "no_core_hit"
            elif any(e["contig"] != c for e in core): why = "core_other_contig"
            elif inside: why = "inside"
            else: why = "core_outside_span"
            dist = max([max(lo - e["end"], e["start"] - hi, 0) for e in core if e["contig"] == c] or [0])
            capped = [cp for (s, e2, cp) in caps.get((x["family"], x["contig"]), []) if s < x["end"] and e2 > x["start"]]
            m = meta.get(g, {})
            out.writerow([panel, g, m.get("PHYLUM"), m.get("CLASS"), m.get("ORDER"), m.get("FAMILY"), m.get("SPECIES"),
                          x["family"], x["contig"], x["start"], x["end"], x["confidence"], x["locus_class"], x["idiomorph"],
                          "low" if inside else "withheld", why, lo + PAD, hi - PAD, len({e["gene"] for e in fl}),
                          ",".join(sorted({e["gene"] for e in fl})), ",".join(sorted({e["gene"] for e in core})), dist,
                          ",".join(str(v) for v in capped), ";".join(fmt(e) for e in core), ";".join(fmt(e) for e in fl)])
print("calls per panel:", tot, file=sys.stderr)
