"""How many copies are needed: random k sequences per class (k = 2..20, 8 draws) build the a1 and alpha2 profiles; accuracy on all other a1/alpha2 sequences. HPCC (mafft)."""
import random, sys
sys.argv = ["x"]
exec(open("loco.py").read().split("CL = [")[0])
random.seed(7)
a1 = [i for i, v in seqs.items() if v[0] == "a1"]; al2 = [i for i, v in seqs.items() if v[0] == "alpha2"]
out = open("data/ablation.tsv", "w"); out.write("k\tdraw\tn_test\tcorrect\tacc\tmin_margin\n")
for k in (2, 3, 4, 6, 8, 12, 18):
    for d in range(8):
        tr1 = random.sample(a1, k); tr2 = random.sample(al2, k)
        h1 = build("a1", [(i, seqs[i][3]) for i in tr1]); h2 = build("alpha2", [(i, seqs[i][3]) for i in tr2])
        ok = n = 0; mm = 999
        for i, v in seqs.items():
            if v[0] in ("a1", "alpha2") and i not in tr1 and i not in tr2:
                s1, s2 = score(h1, v[3]), score(h2, v[3]); n += 1
                good = (s1 > s2) == (v[0] == "a1"); ok += good
                mm = min(mm, (s1 - s2) if v[0] == "a1" else (s2 - s1))
        out.write(f"{k}\t{d}\t{n}\t{ok}\t{ok / n:.3f}\t{mm:.1f}\n"); out.flush()
