import yaml, glob
L = yaml.CSafeLoader
R = "/bigdata/stajichlab/jstajich/projects/MATPredict/results"
A = f"{R}/2026-10-26_v070_submit/attr"
def ev(p):
    d = yaml.load(open(p), Loader=L)
    return {(l.get("contig"), l.get("start"), e["gene"]): (e.get("status"), e.get("start"), e.get("end"), e.get("strand"), e.get("identity")) for l in d.get("detected") or [] for e in l.get("gene_evidence") or []}
for g in ("GCA_030564055.1_ASM3056405v1", "GCA_057550925.1_JCM35544_assembly01"):
    a = ev(f"{A}/new_r1_asco5/runs/{g}/detection_report.yaml"); b = ev(f"{A}/new_r2_asco5/runs/{g}/detection_report.yaml")
    print(g)
    for k in sorted(set(a) | set(b), key=str):
        if a.get(k) != b.get(k): print("  ", k, "r1", a.get(k), "r2", b.get(k))
g = "GCF_022376455.1_Epityp1"
o = glob.glob(f"{R}/2026-10-23_basidiomycota_cleanrun/runs/*/runs/{g}/detection_report.yaml")[0]
d = yaml.load(open(o), Loader=L)
print("old Epityp1 keys", [k for k in d]); print("suppressed:", [(s.get("family"), s.get("contig"), s.get("start"), s.get("reason") or s.get("note")) for s in d.get("suppressed_loci") or []][:6])
print("run.wall", d["run"].get("wall_seconds"))
