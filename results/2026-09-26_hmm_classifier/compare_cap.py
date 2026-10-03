"""Cap-rank check: calls on the same 34 genomes, live-gene rank (541a658) vs
all-distinct-gene rank (7277e10), both with the signed-off Dothideomycetes
records (db/Ascomycota from b7bd2a6). usage: compare_cap.py BASE_RUNS NEW_RUNS"""
import collections, json, os, sys
import yaml
from yaml import CSafeLoader as L

B, N = sys.argv[1], sys.argv[2]


def calls(runs, g):
    p = f"{runs}/{g}/detection_report.yaml"
    if not os.path.exists(p):
        return None
    return [(c["family"], c["contig"], c["start"], c["idiomorph"], c["confidence"])
            for c in (yaml.load(open(p), Loader=L) or {}).get("detected") or []]


def wall(runs, g):
    try:
        return float(open(f"{runs}/{g}/wall_seconds").read())
    except Exception:
        return None


def ov(a, b):
    return a[0] == b[0] and a[1] == b[1] and abs(a[2] - b[2]) < 20000


tot = collections.Counter()
wb, wn = [], []
for g in sorted(set(os.listdir(B)) | set(os.listdir(N))):
    b, n = calls(B, g), calls(N, g)
    if b is None or n is None:
        print(f"  MISSING {g}: base={b is not None} new={n is not None}")
        continue
    if wall(B, g) and wall(N, g):
        wb.append(wall(B, g)); wn.append(wall(N, g))
    lost = [x for x in b if not any(ov(x, y) for y in n)]
    gained = [y for y in n if not any(ov(x, y) for x in b)]
    same = [(x, next(y for y in n if ov(x, y))) for x in b if any(ov(x, y) for y in n)]
    changed = [(x, y) for x, y in same if x[3:] != y[3:]]
    tot["calls_base"] += len(b); tot["calls_new"] += len(n)
    tot["lost"] += len(lost); tot["gained"] += len(gained); tot["label_or_conf_changed"] += len(changed)
    if lost or gained or changed:
        print(f"  {g}: lost {lost} gained {gained} changed {changed}")
print(dict(tot))
if wb:
    wb.sort(); wn.sort()
    print(f"wall s median base {wb[len(wb)//2]:.0f} new {wn[len(wn)//2]:.0f}; total h base {sum(wb)/3600:.1f} new {sum(wn)/3600:.1f}")
