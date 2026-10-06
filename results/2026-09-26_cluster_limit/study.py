"""Would polishing only the top-N admitted clusters per genome save substantial time,
and what would it lose? Read-only study over existing runs.

Per genome: admitted clusters (evidence rows with admitted=true) are the ones
sent to exonerate/miniprot. A reported locus 'came from' an admitted cluster when
they overlap on the same contig. Clusters are ranked by what is known BEFORE
polishing: distinct genes (gene_count), then best_identity, then hit_count.
"""
import collections, glob, json, os, statistics as st, sys
import yaml
from yaml import CSafeLoader as L
R = "/bigdata/stajichlab/jstajich/projects/MATPredict/results"
PANELS = sys.argv[1:] or ["2026-09-25_serinales_all_8dec9e1", "2026-09-25_serinales_all_2c6c3a7",
    "2026-09-23_saccharomyces_bar", "2026-09-25_pezizomycetes_d689647", "2026-09-24_pilots/dothideomycetes",
    "2026-09-24_pilots/pezizomycotina_uncurated", "2026-09-24_pilots/orbiliomycetes",
    "2026-09-24_pilots/dipodascomycetes_etc", "2026-09-24_pilots/saccharomycetes_other"]
NS = (1, 2, 3, 4, 6)
def ov(a, b): return a["contig"] == b["contig"] and not (a["end"] < b["start"] or a["start"] > b["end"])
for panel in PANELS:
    D = f"{R}/{panel}/runs"
    rows = []
    for rp in glob.glob(f"{D}/*/detection_report.yaml"):
        g = os.path.dirname(rp); dp = f"{g}/evidence_diagnostics.jsonl"
        if not os.path.exists(dp): continue
        det = (yaml.load(open(rp), Loader=L) or {}).get("detected") or []
        adm = []
        for l in open(dp):
            e = json.loads(l)
            if e.get("kind") == "evidence" and e.get("admitted"):
                adm.append(dict(contig=e["contig"], start=e["cluster_start"], end=e["cluster_end"],
                                genes=e["gene_count"], ident=e["best_identity"], hits=e["hit_count"]))
        adm.sort(key=lambda c: (-c["genes"], -c["ident"], -c["hits"]))
        # rank (1-based) of the admitted cluster behind each reported locus; None = no admitted cluster overlaps
        ranks = [next((i + 1 for i, c in enumerate(adm) if ov(c, x)), None) for x in det]
        w = None
        if os.path.exists(f"{g}/wall_seconds"):
            try: w = int(open(f"{g}/wall_seconds").read())
            except ValueError: pass
        rows.append(dict(adm=len(adm), work=sum(c["genes"] for c in adm), calls=len(det), ranks=ranks, wall=w,
                         work_by_rank=[c["genes"] for c in adm]))
    if not rows: continue
    n = len(rows); A = [r["adm"] for r in rows]; C = sum(r["calls"] for r in rows)
    used = sum(len({rk for rk in r["ranks"] if rk}) for r in rows)
    print(f"\n## {panel}: {n} genomes, {C} calls, admitted clusters/genome median {st.median(A)} mean {st.mean(A):.1f} max {max(A)}; "
          f"admitted clusters that became a call: {used}/{sum(A)} ({100*used/max(1,sum(A)):.0f}%)")
    print(f"   {'top-N':>6} {'calls lost':>11} {'genomes losing a call':>22} {'polish work saved':>18}")
    for N in NS:
        lost = sum(1 for r in rows for rk in r["ranks"] if rk and rk > N)
        glost = sum(1 for r in rows if any(rk and rk > N for rk in r["ranks"]))
        saved = sum(sum(r["work_by_rank"][N:]) for r in rows); tot = sum(r["work"] for r in rows)
        print(f"   {N:>6} {lost:>5}/{C:<5} {glost:>12}/{n:<8} {100*saved/max(1,tot):>15.0f}%")
    ww = [(r["work"], r["wall"]) for r in rows if r["wall"] is not None]
    if len(ww) > 20:
        xs, ys = zip(*ww); mx, my = st.mean(xs), st.mean(ys)
        sxx = sum((x - mx) ** 2 for x in xs)
        if sxx:
            b = sum((x - mx) * (y - my) for x, y in ww) / sxx; a = my - b * mx
            r2 = 1 - sum((y - (a + b * x)) ** 2 for x, y in ww) / max(1e-9, sum((y - my) ** 2 for y in ys))
            print(f"   wall_s ~ {a:.1f} + {b:.2f} x (admitted-cluster genes); R^2 {r2:.2f}; n={len(ww)} (median wall {st.median(ys):.0f}s)")
