"""Gate C: the Basidiomycota pilot on run-ad1f865 vs the earlier anchors-searched arm."""
import os, statistics, yaml, collections
R = os.path.dirname(os.path.abspath(__file__))
NEW = os.path.join(R, "runs_out", "runs")
OLD = os.path.join(R, "..", "2026-09-26_basidio_anchors", "detect_anchor", "runs")
meta = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in open(os.path.join(R, "..", "2026-09-26_basidio_anchors", "pilot.tsv"))}

def load(d, g):
    p = os.path.join(d, g, "detection_report.yaml")
    return yaml.safe_load(open(p)) or {} if os.path.exists(p) else None

def wall(d, g):
    p = os.path.join(d, g, "wall_seconds")
    return int(open(p).read()) if os.path.exists(p) else None

rows = []; homo = 0; fc_low = 0; fc_withheld = 0; walls = []; walls_by = collections.defaultdict(list)
called = collections.Counter(); total = collections.Counter()
for g in sorted(os.listdir(NEW)):
    r = load(NEW, g); o = load(OLD, g)
    if r is None: rows.append((g, "NO REPORT")); continue
    cls = meta.get(g, [g, "", "Tremellomycetes*"])[2]
    det = r.get("detected") or []
    homo += sum(1 for d in det if d.get("homothallic_candidate"))
    fc_low += sum(1 for d in det if d.get("idiomorph_unmodelled"))
    fc_withheld += r.get("suppressed_flank_carried") or 0
    w = wall(NEW, g); walls.append(w); walls_by[r.get("routing_mode")].append(w)
    total[cls] += 1; called[cls] += bool(det)
    oc = None if o is None else sorted((d["family"], d["idiomorph"], d["confidence"]) for d in (o.get("detected") or []))
    nc = sorted((d["family"], d["idiomorph"], d["confidence"]) for d in det)
    rows.append((g, cls, r.get("routing_mode"), w, oc, nc, r.get("suppressed_flank_carried") or 0))
print("called by class (new):")
for k in total: print(f"  {k:20s} {called[k]}/{total[k]}")
ws = [w for w in walls if w is not None]
print(f"median wall {statistics.median(ws)} s (n={len(ws)}), max {max(ws)}")
for k, v in walls_by.items(): print(f"  routing {k}: n={len(v)} median {statistics.median(v)} s max {max(v)}")
print(f"homothallic_candidate calls: {homo}")
print(f"flank-carried: kept at low {fc_low}, withheld {fc_withheld}")
print("\nper genome (old anchor arm -> new):")
for r in rows:
    if len(r) == 2: print(" ", *r); continue
    g, cls, rm, w, oc, nc, fcw = r
    flag = "" if oc == nc else "  <-- CHANGED"
    print(f"  {g[:40]:40s} {cls[:16]:16s} {rm:16s} {w}s old={oc} new={nc} fc_withheld={fcw}{flag}")
