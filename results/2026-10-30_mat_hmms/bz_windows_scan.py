"""All five Barnettozyma windows (four B. californica assemblies, B. discipulorum): best stop-to-stop segment per profile (a2, a1, alpha2, alpha1)."""
import glob, sys
sys.argv = ["x"]
src = open("bc_alpha2_scan.py").read().split('D = "../2026')[0]
exec(src)
for n in ("a2", "sacch_a2"):
    p = "../../db/Ascomycota/classifiers/MATa2_saccharomycotina/MATa2_saccharomycotina.hmm" if n == "sacch_a2" else "hmm/a2.hmm"
    with pyhmmer.plan7.HMMFile(p) as hf: H[n] = hf.read()
D = "../2026-10-30_barnettozyma/data/"
out = open("data/bz_windows_scan.tsv", "w"); out.write("window\tprofile\tbest_segment\tlen\tbits\n")
for fa in sorted(glob.glob(D + "GCA_*_*.fa")):
    seq = next(iter(fasta(fa).values())); sg = [(k, v) for k, v in segs(seq, 1, len(seq), 1) if len(v) >= 90]
    rr = score(sg)
    for pr in ("a2", "sacch_a2", "a1", "alpha2", "alpha1"):
        k, v = max(rr.items(), key=lambda kv: kv[1].get(pr, -9))
        line = f"{fa.split('/')[-1][:28]}\t{pr}\t{k}\t{len(dict(sg)[k])}\t{v.get(pr, 0):.1f}"
        print(line); out.write(line + "\n")
