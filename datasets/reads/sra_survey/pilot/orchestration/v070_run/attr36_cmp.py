import glob, os, sys, yaml
sys.path.insert(0, "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-26_v070_submit")
from full_compare import call_sig, model_sig, L, R
A = f"{R}/2026-10-26_v070_submit/attr"
def rep(path):
    return yaml.load(open(path), Loader=L) or {}
def find(pats, g):
    for p in pats:
        m = glob.glob(p.format(g=g))
        if m: return m[0]
PN = [f"{R}/2026-10-26_v070_ascomycota/wave_*/small_*/runs/{{g}}/detection_report.yaml", f"{R}/2026-10-26_v070_basidiomycota/runs/*/runs/{{g}}/detection_report.yaml"]
PO = [f"{R}/2026-10-23_ascomycota_cleanrun/wave_*/small_*/runs/{{g}}/detection_report.yaml", f"{R}/2026-10-23_basidiomycota_cleanrun/runs/*/runs/{{g}}/detection_report.yaml"]
def sig(d): return (call_sig(d), model_sig(d))
def label(a, b):
    return "same" if sig(a) == sig(b) else ("model" if call_sig(a) == call_sig(b) else "call")
for ctl, lst in (("cap5k_all36", "list_all36.tsv"), ("old10b_all36", "list_all36.tsv"), ("new_r1_asco5", "list_asco5.tsv"), ("new_r2_asco5", "list_asco5.tsv")):
    print("==", ctl)
    rows = []
    for g in [l.split("\t")[0] for l in open(f"{A}/{lst}")]:
        cp = f"{A}/{ctl}/runs/{g}/detection_report.yaml"
        if not os.path.exists(cp) or os.path.getsize(cp) == 0:
            rows.append((g, "NO_REPORT", "", "")); continue
        c, n, o = rep(cp), rep(find(PN, g)), rep(find(PO, g))
        rows.append((g, "control_vs_new=" + label(c, n), "control_vs_old10b=" + label(c, o), "new_vs_old=" + label(n, o)))
    import collections
    print(" ", dict(collections.Counter((r[1], r[2], r[3]) for r in rows)))
    if "asco5" in ctl or "old10b" in ctl:
        for r in rows[:5]: print("  ", r)
print("\n== split: Basidiomycota 31 vs Ascomycota 5; r1 vs r2")
import collections
a5 = [l.split("\t")[0] for l in open(f"{A}/list_asco5.tsv")]
b31 = [l.split("\t")[0] for l in open(f"{A}/list_basidio31.tsv")]
def ctlrep(ctl, g):
    p = f"{A}/{ctl}/runs/{g}/detection_report.yaml"
    return rep(p) if os.path.exists(p) and os.path.getsize(p) else None
for ctl in ("cap5k_all36", "old10b_all36"):
    for nm, gl in (("basidio31", b31), ("asco5", a5)):
        c = collections.Counter()
        for g in gl:
            cr = ctlrep(ctl, g); o = rep(find(PO, g)); n = rep(find(PN, g))
            c["ctl==old10b:" + label(cr, o)] += 1
            c["ctl==new:" + label(cr, n)] += 1
        print(ctl, nm, dict(c))
for g in a5:
    r1, r2 = ctlrep("new_r1_asco5", g), ctlrep("new_r2_asco5", g)
    print("r1 vs r2", g, label(r1, r2), "| r1 vs new", label(r1, rep(find(PN, g))), "| r2 vs old", label(r2, rep(find(PO, g))))
print("\n== basidio call diffs and r1/r2 detail")
def calls(d):
    return sorted((l.get("family"), l.get("contig"), l.get("start"), l.get("end"), l.get("idiomorph"), l.get("confidence"), l.get("polished_genes")) for l in d.get("detected") or [])
for g in [l.strip() for l in open(f"{R}/2026-10-26_v070_submit/compare/diff_basidio_call.txt")]:
    n, o, c5, ob = rep(find(PN, g)), rep(find(PO, g)), ctlrep("cap5k_all36", g), ctlrep("old10b_all36", g)
    print(g)
    for nm, d in (("old10b original", o), ("new", n), ("cap5k control", c5), ("old10b rerun", ob)):
        print("   ", nm, calls(d))
for g in a5[:1]:
    for ctl in ("new_r1_asco5", "new_r2_asco5"):
        d = ctlrep(ctl, g)
        print(g, ctl, [(e["gene"], e.get("start"), e.get("end"), e.get("strand")) for l in d["detected"] for e in l["gene_evidence"]][:4])
