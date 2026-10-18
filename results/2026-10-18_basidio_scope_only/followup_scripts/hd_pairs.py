import json, glob, collections
W = "/bigdata/stajichlab/jstajich/projects/mat-fam18/rc/hd/"
out = {}
for f in sorted(glob.glob(W + "*.json")):
    d = json.load(open(f))
    hs = sorted(d["hits"], key=lambda h: (h[0], h[2]))
    # merge Homeodomain/KN hits of one gene: same contig+strand within 300 bp
    loci = []
    for h in hs:
        if loci and loci[-1][0] == h[0] and loci[-1][1] == h[1] and h[2] - loci[-1][3] <= 300:
            loci[-1][3] = max(loci[-1][3], h[3]); loci[-1][4].add(h[4])
        else: loci.append([h[0], h[1], h[2], h[3], {h[4]}])
    pairs = []
    for i, a in enumerate(loci):
        for b in loci[i + 1:]:
            if b[0] != a[0] or b[2] - a[3] > 4000: break
            if a[1] == "-" and b[1] == "+" and b[2] > a[3]:
                pairs.append((a[0], a[2], a[3], "".join(sorted(a[4])), b[2], b[3], "".join(sorted(b[4])), b[2] - a[3]))
    print(f.split("/")[-1][:-5], "orfs", d["orfs"], "hd_hits", d["hd_hits"], "hd_loci", len(loci), "divergent_pairs<=4kb", len(pairs))
    for p in pairs: print("    ", p)
