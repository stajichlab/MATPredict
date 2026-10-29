#!/usr/bin/env python3
"""Flag the genomes whose loci come from two Saccharomycotina MAT systems (for example MATsc and MTL), and say what looks wrong, from the committed clean-run tables.

    idiomorph_crosscall_flags.py --genomes results/<run>/genomes.tsv --loci results/<run>/loci.tsv --out flags.tsv [--genera Ogataea Pichia]
One row per genome that has loci of two or more of the families MATsc, MTL, MATog, MATkp. `same_locus_calls` counts the loci of the secondary family whose contig, start and end are identical to a locus of the primary
family (one locus called by two panels, not a second locus); `routing` is how the genome was routed (`phylum_fallback`: no lineage record matched, so every family of the phylum was tried). Standard library only.
Suggested action: `duplicate_family_call` (dedupe: keep the better-scoring family, flag the other), `reroute` (a yeast genus with its own family but a phylum-fallback genome), `curate` (secondary loci at other positions).
"""
import argparse, collections, csv

FAMS = ("Ascomycota:MATsc", "Ascomycota:MTL", "Ascomycota:MATog", "Ascomycota:MATkp")
OWN_FAMILY = {"Ogataea": "Ascomycota:MATog", "Komagataella": "Ascomycota:MATkp"}


def rd(p):
    return list(csv.DictReader(open(p, newline=""), delimiter="\t"))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0]); ap.add_argument("--genomes", required=True); ap.add_argument("--loci", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--genera", nargs="*"); a = ap.parse_args()
    G = {r["genome"]: r for r in rd(a.genomes)}; by = collections.defaultdict(list)
    for r in rd(a.loci):
        if r["family_called"] in FAMS and r["genome"] in G: by[r["genome"]].append(r)
    rows = []
    for gn, v in sorted(by.items()):
        fams = collections.Counter(x["family_called"] for x in v)
        if len(fams) < 2: continue
        g = G[gn]; genus = (g["species"].split() or [""])[0]
        if a.genera and genus not in a.genera: continue
        prim = fams.most_common(1)[0][0]; sec = [f for f in fams if f != prim]
        pcoords = {(x["contig"], x["start"], x["end"]) for x in v if x["family_called"] == prim}
        sl = [x for x in v if x["family_called"] in sec]
        same = sum(1 for x in sl if (x["contig"], x["start"], x["end"]) in pcoords)
        own = OWN_FAMILY.get(genus)
        act = "duplicate_family_call" if same == len(sl) else "curate"
        if g["routing"] == "phylum_fallback" and own and own != prim: act += "+reroute"
        rows.append(dict(genome=gn, species=g["species"], order=g["order"], routing=g["routing"], primary_family=prim.replace("Ascomycota:", ""), secondary_families=",".join(f.replace("Ascomycota:", "") for f in sec),
                         secondary_loci=len(sl), same_locus_calls=same, secondary_classes=",".join(sorted({x["locus_class"] for x in sl})), secondary_genes="|".join(sorted({x["genes_found"] for x in sl}))[:80], suggested_action=act))
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
    c = collections.Counter(r["suggested_action"] for r in rows)
    print(len(rows), "genomes;", dict(c), "; routing", dict(collections.Counter(r["routing"] for r in rows)))


if __name__ == "__main__":
    main()
