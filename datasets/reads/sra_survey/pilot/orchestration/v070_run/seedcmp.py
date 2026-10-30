import sys, yaml, os
sys.path.insert(0, "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-26_v070_submit")
from full_compare import model_sig, call_sig, L, R
A = f"{R}/2026-10-26_v070_submit/attr"
g5 = [l.split("\t")[0] for l in open(f"{A}/list_asco5.tsv")]
def sg(ctl, g):
    return model_sig(yaml.load(open(f"{A}/{ctl}/runs/{g}/detection_report.yaml"), Loader=L))
for s in ("0", "1"):
    print("seed", s, "a==b per genome:", [sg(f"seed{s}_a", g) == sg(f"seed{s}_b", g) for g in g5])
print("seed0 vs seed1 per genome same:", [sg("seed0_a", g) == sg("seed1_a", g) for g in g5])
