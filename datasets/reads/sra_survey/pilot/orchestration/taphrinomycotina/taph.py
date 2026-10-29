import csv,collections,statistics
W="/bigdata/stajichlab/jstajich/projects/MATPredict/results/sra_reads_typing/wave_inputs"
G={r["genome"]:r for r in csv.DictReader(open(W+"/asco_genomes.tsv"),delimiter="\t")}
for sp in ("Schizosaccharomyces pombe","Schizosaccharomyces japonicus","Pneumocystis jirovecii","Debaryomyces hansenii","Zygosaccharomyces sapae"):
    rows=[r for r in csv.DictReader(open(W+"/asco_loci.tsv"),delimiter="\t") if r["genome"] in G and G[r["genome"]]["species"].startswith(sp)]
    n=sum(1 for g in G.values() if g["species"].startswith(sp))
    sp_=[int(r["end"])-int(r["start"]) for r in rows]
    print(sp,"genomes",n,"loci",len(rows),"classes",dict(collections.Counter(r["locus_class"] for r in rows)),"idiomorphs",dict(collections.Counter(r["idiomorph"] for r in rows)),"span min/median/max",(min(sp_),int(statistics.median(sp_)),max(sp_)) if sp_ else None, "families",dict(collections.Counter(r["family_called"] for r in rows)))
