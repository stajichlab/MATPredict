"""Synteny check of every phylum_fallback call where an anchor applies.

Anchor = the genome's own ortholog, taken as the top tblastn hit (bitscore)
of the anchor protein, found independently of the pipeline:
  Ascomycota MAT / MATsc / MATyl / MATtub / PM: SLA2 or APN2
      (hits from results/2026-09-24_bar_synteny/hits, A. fumigatus + Diaporthales queries)
  Basidiomycota HD: MIP1 or beta-fg (Agaricus H97 queries, hd_hits/)
A call is SUPPORTED when an anchor top hit lies on the same contig within
WINDOW bp of the call span; NOT SUPPORTED otherwise; NOT TESTABLE when the
family has no anchor (bLocus, aLocus, Tremellales-type MAT, Mucoromycota,
MTL -- its PAP1/OBP1/PIK1 flanks are roster genes, so the test is circular)
or the genome has no anchor hit file.
Writes adjacency.tsv and adjacency_summary.txt.
"""
import csv, collections, os, sys

WINDOW = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
R = "/bigdata/stajichlab/jstajich/projects/MATPredict/results"
ASCO_HITS = f"{R}/2026-09-24_bar_synteny/hits"
HD_HITS = f"{R}/2026-09-27_fallback_confidence_floor/hd_hits"
ANCHORS = {
    "Ascomycota:MAT": (ASCO_HITS, ("SLA2", "APN2")),
    "Ascomycota:MATsc": (ASCO_HITS, ("SLA2", "APN2")),
    "Ascomycota:MATyl": (ASCO_HITS, ("SLA2", "APN2")),
    "Ascomycota:MATtub": (ASCO_HITS, ("SLA2", "APN2")),
    "Ascomycota:PM": (ASCO_HITS, ("SLA2", "APN2")),
    "Basidiomycota:HD": (HD_HITS, ("MIP1", "beta_fg")),
}


def top_hits(path, genes):
    best = {}
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if len(f) < 8 or f[0] == "none":
            continue
        gene = f[0].split("|")[0]
        if gene not in genes:
            continue
        bs = float(f[7])
        if gene not in best or bs > best[gene][0]:
            s, e = sorted((int(f[4]), int(f[5])))
            best[gene] = (bs, f[1], s, e)
    return best


rows = [r for r in csv.DictReader(open("calls.tsv"), delimiter="\t") if r["routing"] == "phylum_fallback"]
out = []
for r in rows:
    fam = r["family"]
    status, dist, anchor = "not_testable", "", ""
    if fam in ANCHORS:
        d, genes = ANCHORS[fam]
        p = os.path.join(d, r["genome"] + ".tsv")
        if not os.path.exists(p):
            status = "no_anchor_file"
        else:
            best = top_hits(p, genes)
            # distance from the best-identity CORE gene, not the call span:
            # the span includes SLA2/APN2 when they are roster flanks, which
            # would make the test circular.
            if r["core_start"] == "":
                out.append({**r, "synteny": "no_core_coords", "anchor": "", "anchor_gap_bp": ""}); continue
            s0, e0 = sorted((int(r["core_start"]), int(r["core_end"])))
            ds = []
            for g, (bs, contig, s, e) in best.items():
                if contig != (r["core_contig"] or r["contig"]):
                    continue
                gap = max(0, s - e0, s0 - e)
                ds.append((gap, g))
            if ds:
                gap, g = min(ds)
                dist, anchor = gap, g
                status = "supported" if gap <= WINDOW else "not_supported"
            else:
                status = "not_supported" if best else "no_anchor_hit"
    out.append({**r, "synteny": status, "anchor": anchor, "anchor_gap_bp": dist})

with open(f"adjacency_{WINDOW}.tsv", "w", newline="") as fo:
    w = csv.DictWriter(fo, fieldnames=list(out[0]), delimiter="\t"); w.writeheader(); w.writerows(out)

# summary: synteny by identity band, per run group and family
def band(v):
    if v == "":
        return "none"
    v = float(v)
    return "<30" if v < 30 else "30-35" if v < 35 else "35-40" if v < 40 else "40-50" if v < 50 else ">=50"


BANDS = ["<30", "30-35", "35-40", "40-50", ">=50"]
with open(f"adjacency_summary_{WINDOW}.txt", "w") as fo:
    def p(*a):
        print(*a); print(*a, file=fo)
    p(f"window {WINDOW} bp; best_core_any identity bands; cells = supported/tested (not_testable excluded)")
    groups = collections.defaultdict(lambda: collections.Counter())
    for r in out:
        if r["synteny"] not in ("supported", "not_supported"):
            continue
        grp = "Ascomycota" if r["family"].startswith("Ascomycota") else r["family"]
        for key in ((grp, r["confidence"]), (grp, "all")):
            b = band(r["best_core_any"])
            groups[key][(b, "n")] += 1
            if r["synteny"] == "supported":
                groups[key][(b, "s")] += 1
    for key in sorted(groups):
        c = groups[key]
        cells = [f"{b}:{c[(b,'s')]}/{c[(b,'n')]}" for b in BANDS]
        p(f"  {key[0]:18s} {key[1]:6s} " + "  ".join(cells))
    p("\nstatus counts by family:")
    sc = collections.Counter((r["family"], r["synteny"]) for r in out)
    for k, v in sorted(sc.items()):
        p(f"  {k[0]:24s} {k[1]:15s} {v}")
