"""Score the Zygo 23 run against zygo_truth.tsv (scaffold coordinates).

The genomes searched are the funannotate .contigs.fsa files (contig_N names);
each organism's .agp maps contigs onto the scaffolds the truth table uses.
A call scores when, translated to scaffold coordinates, it lies on the truth
scaffold and overlaps the truth span, with the truth idiomorph.
"""
import glob, os, sys, yaml
R = os.path.dirname(os.path.abspath(__file__))
RUNS = sys.argv[1] if len(sys.argv) > 1 else os.path.join(R, "zygo", "runs")
T = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-23_zygo_bar/zygo_truth.tsv"

def agp_map(fsa):
    m = {}
    for agp in glob.glob(os.path.join(os.path.dirname(fsa), "*.agp")):
        for line in open(agp):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 9 and f[4] == "W":
                m[f[5]] = (f[0], int(f[1]), int(f[6]), f[8])
    return m

loc = idio = n = 0
for line in open(T):
    org, scaf, s, e, truth, fsa = line.rstrip("\n").split("\t")[:6]
    s, e = int(s), int(e); n += 1
    amap = agp_map(fsa)
    p = os.path.join(RUNS, org, "detection_report.yaml")
    det = ((yaml.safe_load(open(p)) or {}).get("detected") or []) if os.path.exists(p) else []
    hits = []
    for d in det:
        sc, obj_start, comp_start, strand = amap.get(d["contig"], (d["contig"], 1, 1, "+"))
        if strand != "+":
            continue  # none of the 23 truth contigs is reverse-placed; flag if one is
        a = obj_start + d["start"] - comp_start; b = obj_start + d["end"] - comp_start
        if sc == scaf and a <= e and b >= s:
            hits.append(d)
    ok_loc = bool(hits); ok_id = any(d["idiomorph"] == truth for d in hits)
    loc += ok_loc; idio += ok_id
    if not (ok_loc and ok_id):
        print("MISS", org, truth, [(d["contig"], d["start"], d["idiomorph"]) for d in det])
print(f"locus on truth scaffold {loc}/{n}; idiomorph correct {idio}/{n}")
