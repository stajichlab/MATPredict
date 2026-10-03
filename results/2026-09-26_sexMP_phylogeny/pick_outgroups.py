"""Pick same-genome non-MAT HMG-box copies as outgroups.

One genome per genus (the one with the most sexM/sexP loci), and in it the three
highest-scoring HMG copies that overlap no reported locus. Chytrids are excluded.
Also writes per-genome counts of HMG copies (E <= 1e-3) inside/outside loci.
"""
import csv, collections
tax = {r["ASMID"]: r for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"))}
nloci = collections.Counter(r["genome"] for r in csv.DictReader(open("loci.tsv"), delimiter="\t"))
cl = list(csv.DictReader(open("clusters.tsv"), delimiter="\t"))
counts = collections.defaultdict(lambda: [0, 0])
for r in cl:
    counts[r["genome"]][0 if r["locus_ids"] else 1] += 1
with open("hmg_copies_per_genome.tsv", "w") as fo:
    fo.write("genome\tphylum\tfamily\tgenus\tcopies_in_loci\tcopies_not_included\n")
    for g, (a, b) in sorted(counts.items()):
        t = tax.get(g, {})
        fo.write(f"{g}\t{t.get('PHYLUM','')}\t{t.get('FAMILY','')}\t{t.get('GENUS','')}\t{a}\t{b}\n")
by_genus = collections.defaultdict(list)
for g in nloci:
    by_genus[tax.get(g, {}).get("GENUS", "?")].append(g)
chosen = {max(gs, key=lambda g: (nloci[g], g)) for gs in by_genus.values()}
out = [r for r in cl if r["genome"] in chosen and not r["locus_ids"]]
pick = []
for g in sorted(chosen):
    rs = sorted((r for r in out if r["genome"] == g), key=lambda r: -float(r["best_bits"]))[:3]
    pick += rs
with open("outgroup_clusters.tsv", "w", newline="") as fo:
    w = csv.DictWriter(fo, fieldnames=list(cl[0]), delimiter="\t"); w.writeheader(); w.writerows(pick)
print("genera", len(by_genus), "outgroup copies", len(pick))
