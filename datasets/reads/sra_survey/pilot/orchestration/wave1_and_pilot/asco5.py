import yaml, glob, sys
L = yaml.CSafeLoader
R = "/bigdata/stajichlab/jstajich/projects/MATPredict/results"
D = f"{R}/2026-10-26_v070_submit/compare"
for g in open(f"{D}/diff_asco_model.txt").read().split():
    n = glob.glob(f"{R}/2026-10-26_v070_ascomycota/wave_*/small_*/runs/{g}/detection_report.yaml")[0]
    o = glob.glob(f"{R}/2026-10-23_ascomycota_cleanrun/wave_*/small_*/runs/{g}/detection_report.yaml")[0]
    dn, do = yaml.load(open(n), Loader=L), yaml.load(open(o), Loader=L)
    def ev(d):
        return {(l.get("contig"), l.get("start"), e["gene"]): (e.get("status"), e.get("start"), e.get("end"), e.get("strand"), tuple((x["start"], x["end"]) for x in e.get("exons") or []), e.get("identity"), e.get("frameshifts"), e.get("polish_tool"))
                for l in d.get("detected") or [] for e in l.get("gene_evidence") or []}
    a, b = ev(dn), ev(do)
    print("==", g, n.split("/")[-5:-3])
    for k in sorted(set(a) | set(b), key=str):
        if a.get(k) != b.get(k):
            print(" ", k, "\n    new", a.get(k), "\n    old", b.get(k))
