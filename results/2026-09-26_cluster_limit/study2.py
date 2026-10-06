"""Projected wall time under a top-N limit, ranking per genome vs per family."""
import glob, json, os, statistics as st
import yaml
from yaml import CSafeLoader as L
R = "/bigdata/stajichlab/jstajich/projects/MATPredict/results"
PANELS = ["2026-09-25_serinales_all_8dec9e1", "2026-09-23_saccharomyces_bar", "2026-09-24_pilots/dothideomycetes",
          "2026-09-24_pilots/pezizomycotina_uncurated", "2026-09-24_pilots/orbiliomycetes",
          "2026-09-24_pilots/dipodascomycetes_etc", "2026-09-24_pilots/saccharomycetes_other"]
def ov(a, b): return a["contig"] == b["contig"] and not (a["end"] < b["start"] or a["start"] > b["end"])
key = lambda c: (-c["genes"], -c["ident"], -c["hits"])
for panel in PANELS:
    rows = []
    for rp in glob.glob(f"{R}/{panel}/runs/*/detection_report.yaml"):
        g = os.path.dirname(rp); dp = f"{g}/evidence_diagnostics.jsonl"
        if not os.path.exists(dp): continue
        det = (yaml.load(open(rp), Loader=L) or {}).get("detected") or []
        adm = [dict(fam=e["family"], contig=e["contig"], start=e["cluster_start"], end=e["cluster_end"], genes=e["gene_count"],
                    ident=e["best_identity"], hits=e["hit_count"]) for e in map(json.loads, open(dp))
               if e.get("kind") == "evidence" and e.get("admitted")]
        w = int(open(f"{g}/wall_seconds").read()) if os.path.exists(f"{g}/wall_seconds") else None
        rows.append((det, adm, w))
    ww = [(sum(c["genes"] for c in a), w) for d, a, w in rows if w is not None]
    if ww:
        xs, ys = zip(*ww); mx, my = st.mean(xs), st.mean(ys); b = sum((x-mx)*(y-my) for x, y in ww)/sum((x-mx)**2 for x in xs); a0 = my - b*mx
        print(f"\n## {panel}  (wall ~ {a0:.0f} + {b:.2f}/gene; current median wall {st.median(ys):.0f}s)")
    else:
        a0, b = float("nan"), float("nan"); print(f"\n## {panel}  (no wall_seconds recorded; call loss only)")
    for mode in ("per-genome", "per-family"):
        for N in (2, 3, 4, 6):
            lost = 0; kept_work = []; walls = []
            for det, adm, w in rows:
                if mode == "per-genome":
                    kept = sorted(adm, key=key)[:N]
                else:
                    kept = []
                    for f in {c["fam"] for c in adm}:
                        kept += sorted([c for c in adm if c["fam"] == f], key=key)[:N]
                lost += sum(1 for x in det if any(ov(c, x) for c in adm) and not any(ov(c, x) for c in kept))
                kw = sum(c["genes"] for c in kept); walls.append(a0 + b*kw)
            print(f"   {mode:10s} N={N}: calls lost {lost:4d}/{sum(len(d) for d,_,_ in rows):<5d}  projected median wall {st.median(walls):5.0f}s")
