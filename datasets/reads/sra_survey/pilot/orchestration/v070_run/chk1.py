import sys, glob, yaml
sys.path.insert(0, "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-26_v070_submit")
from full_compare import call_sig, model_sig, L, R
A = f"{R}/2026-10-26_v070_submit/attr"
g = [x for x in open(f"{A}/list_basidio31.tsv").read().split("\n") if x.startswith("GCA_011057475")][0].split("\t")[0]
def rep(p): return yaml.load(open(p), Loader=L)
n = rep(glob.glob(f"{R}/2026-10-26_v070_basidiomycota/runs/*/runs/{g}/detection_report.yaml")[0])
o = rep(glob.glob(f"{R}/2026-10-23_basidiomycota_cleanrun/runs/*/runs/{g}/detection_report.yaml")[0])
c = rep(f"{A}/cap5k_all36/runs/{g}/detection_report.yaml"); ob = rep(f"{A}/old10b_all36/runs/{g}/detection_report.yaml")
S = lambda d: model_sig(d)
print(g, "cap5k==old", S(c)==S(o), "cap5k==new", S(c)==S(n), "old10brerun==old", S(ob)==S(o), "old10brerun==new", S(ob)==S(n))
