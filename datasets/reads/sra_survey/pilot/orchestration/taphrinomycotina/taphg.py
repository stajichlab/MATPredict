import csv,collections
W="/bigdata/stajichlab/jstajich/projects/MATPredict/results/sra_reads_typing/wave_inputs"
G=list(csv.DictReader(open(W+"/asco_genomes.tsv"),delimiter="\t"))
cls=collections.Counter((r["class_"],r["order"]) for r in G if r["class_"] in ("Taphrinomycetes","Schizosaccharomycetes","Pneumocystidomycetes","Neolectomycetes","Archaeorhizomycetes","Taphrinomycotina","NA") or r["order"] in ("Taphrinales","Schizosaccharomycetales","Pneumocystidales","Neolectales","Archaeorhizomycetales","Protomycetales","Saitoellales"))
print(cls.most_common(20))
tp=[r for r in G if r["order"] in ("Taphrinales","Schizosaccharomycetales","Pneumocystidales","Neolectales","Archaeorhizomycetales","Protomycetales","Saitoellales") or r["class_"] in ("Taphrinomycetes","Schizosaccharomycetes","Pneumocystidomycetes","Neolectomycetes","Archaeorhizomycetes") or r["genus"] in ("Saitoella","Taphrina","Protomyces","Neolecta","Schizosaccharomyces","Pneumocystis","Archaeorhizomyces","Saitozyma","Pneumocystis")]
print(len(tp),collections.Counter(r["species"] for r in tp).most_common(40))
print(collections.Counter(r["status"] for r in tp))
