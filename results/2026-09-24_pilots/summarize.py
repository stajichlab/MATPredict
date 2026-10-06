"""Summarise the 2026-09-24 pilot panels: call rate, loci/genome, bar losses, wall time."""
import collections, csv, glob, os, statistics as st, sys
import yaml
try:
    from yaml import CSafeLoader as L
except ImportError:
    from yaml import SafeLoader as L

P = os.path.dirname(os.path.abspath(__file__))
meta = {r["ASMID"]: r for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"))}

def q(v, p):
    v = sorted(v); return v[min(len(v) - 1, int(p * len(v)))]

for d in sorted(glob.glob(f"{P}/*/runs")):
    name = d.split("/")[-2]
    by_order = collections.defaultdict(lambda: collections.Counter())
    tot = collections.Counter(); conf = collections.Counter(); idio = collections.Counter()
    fam = collections.Counter(); route = collections.Counter(); wall = []; loci = []
    bar_only = 0; reasons = collections.Counter()
    for g in sorted(os.listdir(d)):
        rp = f"{d}/{g}/detection_report.yaml"
        if not os.path.exists(rp) or os.path.getsize(rp) == 0:
            tot["no_report"] += 1; continue
        doc = yaml.load(open(rp), Loader=L) or {}
        det = doc.get("detected") or []; sup = doc.get("suppressed_loci") or []
        o = meta.get(g, {}).get("ORDER") or "?"
        tot["genomes"] += 1; route[doc.get("routing_mode")] += 1; loci.append(len(det))
        by_order[o]["genomes"] += 1
        if det:
            tot["called"] += 1; by_order[o]["called"] += 1
        elif sup:
            bar_only += 1; by_order[o]["bar_only"] += 1
        else:
            for n in doc.get("not_detected") or []:
                reasons[n["reason"].split(";")[0][:70]] += 1
        for x in det:
            conf[x["confidence"]] += 1; idio[str(x["idiomorph"])] += 1; fam[x["family"]] += 1
        w = f"{d}/{g}/wall_seconds"
        if os.path.exists(w):
            try: wall.append(int(open(w).read()))
            except ValueError: pass
    n = tot["genomes"] or 1
    print(f"\n## {name}: genomes={tot['genomes']} no_report={tot['no_report']} routing={dict(route)}")
    print(f"  called={tot['called']} ({100*tot['called']/n:.0f}%)  bar-withheld only={bar_only}  "
          f"loci/genome={st.mean(loci) if loci else 0:.2f} dist={dict(sorted(collections.Counter(loci).items()))}")
    print(f"  confidence={dict(conf)} families={dict(fam)}")
    print(f"  idiomorph={dict(idio.most_common(8))}")
    if wall:
        print(f"  wall_s n={len(wall)} median={st.median(wall):.0f} p90={q(wall,.9)} max={max(wall)} "
              f"sum_slot_h={sum(wall)/3600:.1f}")
    print("  by order: " + "; ".join(f"{o} {c['called']}/{c['genomes']}" + (f" (+{c['bar_only']} bar)" if c['bar_only'] else "")
                                    for o, c in sorted(by_order.items(), key=lambda x: -x[1]['genomes'])))
    print(f"  top not-detected reasons (uncalled, not bar): {dict(reasons.most_common(4))}")
