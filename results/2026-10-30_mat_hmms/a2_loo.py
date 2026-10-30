"""a2 profiles with a clade held out (the control and the tier-3 genes are then truly external to training). HPCC (mafft)."""
import sys
sys.argv = ["x"]
exec(open("loco.py").read().split("CL = [")[0])   # reuse seqs, build, score
for hold in ("Pichiaceae", "Saccharomycetaceae", "Serinales"):
    items = [(i, v[3]) for i, v in seqs.items() if v[0] == "a2" and v[2] != hold]
    h = build("a2_no_" + hold, items)
    with open(f"hmm/a2_no_{hold}.hmm", "wb") as fh: h.write(fh)
    print(hold, len(items), h.M)
