"""Map BFD Basidiomycota genomes to their funannotate proteome files.
Dir name = '<SPECIES_IN>_<STRAIN>' with spaces -> '_' (as funannotate names it);
falls back to '<SPECIES_IN>' when strain is empty. Suppressed ASMIDs skipped."""
import csv, os, re, collections
B = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD"
sup = set()
for l in open(f"{B}/data/curation/suppress.txt"):
    a = l.split(",")[0].strip()
    if a and a != "ASMID": sup.add(a)
rows = [r for r in csv.DictReader(open(f"{B}/samples.csv")) if r["PHYLUM"] == "Basidiomycota" and r["ASMID"] not in sup]
def name(r):
    s = r["SPECIES_IN"].strip()
    st = r["STRAIN"].strip()
    n = f"{s} {st}" if st else s
    return re.sub(r"[ /]+", "_", n)
out = open("proteomes.tsv", "w"); out.write("ASMID\tLOCUSTAG\tSPECIES\tSTRAIN\tCLASS\tORDER\tproteome\n")
cnt = collections.Counter(); tot = collections.Counter()
for r in rows:
    n = name(r)
    p = f"{B}/genome_annotation/{n}/predict_results/{n}.proteins.fa"
    ok = os.path.exists(p)
    tot[(r["CLASS"], r["ORDER"])] += 1
    if ok: cnt[(r["CLASS"], r["ORDER"])] += 1
    out.write("\t".join([r["ASMID"], r["LOCUSTAG"], r["SPECIES_IN"], r["STRAIN"], r["CLASS"], r["ORDER"], p if ok else ""]) + "\n")
print("Basidiomycota genomes", sum(tot.values()), "with proteome", sum(cnt.values()))
for k, v in sorted(tot.items(), key=lambda kv: -kv[1])[:40]:
    print(f"{k[0]:25s} {k[1]:25s} {cnt[k]:5d}/{v}")
