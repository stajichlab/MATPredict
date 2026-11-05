"""shared_flank_check on every pair of called loci with complementary idiomorphs (MATa-type vs MATalpha-type) in each of the 80 genomes, whatever their
family (run on the HPCC with the detector code of PR #264: PYTHONPATH=<worktree>/src python3 flank_pairs.py). Writes flank_pairs.tsv."""
import glob, yaml, sys, itertools
from types import SimpleNamespace
from MATPredict.detect.shared_flank import shared_flank_check
R = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-11-04_189_secondary_loci/run/out/runs/"
LIB = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes/"
side = lambda v: "a" if (v or "").lower() in ("a", "mata") else "alpha" if (v or "").lower() in ("alpha", "matalpha") else None
rows = ["genome\tlocus_a\tlocus_b\tsame_family\tstatus\taligned_bp\tbest_identity"]
for f in sorted(glob.glob(R + "*/detection_report.yaml")):
    g = f.split("/")[-2]; d = yaml.safe_load(open(f))
    L = [l for l in d["detected"] if side(l["idiomorph"])]
    for a, b in itertools.combinations(L, 2):
        if side(a["idiomorph"]) == side(b["idiomorph"]) or (a["contig"] == b["contig"] and a["start"] <= b["end"] and b["start"] <= a["end"]): continue
        c = shared_flank_check(LIB + g + ".fa.gz", SimpleNamespace(contig=a["contig"], start=a["start"], end=a["end"]), SimpleNamespace(contig=b["contig"], start=b["start"], end=b["end"]))
        rows.append("\t".join(str(x) for x in (g, f"{a['family'].split(':')[1]}:{a['idiomorph']}:{a['contig']}:{a['start']}-{a['end']}", f"{b['family'].split(':')[1]}:{b['idiomorph']}:{b['contig']}:{b['start']}-{b['end']}",
                                                a["family"] == b["family"], c["status"], c.get("aligned_bp", ""), c.get("best_identity", ""))))
open("flank_pairs.tsv", "w").write("\n".join(rows) + "\n"); print(len(rows) - 1, "pairs")
