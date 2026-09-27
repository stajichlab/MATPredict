"""Pick ~50 uncalled, capped, phylum-fallback Basidiomycota genomes (<750 Mb),
stratified by order, most-capped first; write ASMID<TAB>TAXID lists sized for
2 h short jobs (16-way), assuming cap-off wall = 3x the recorded capped wall."""
import csv, collections, os
B = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-26_basidiomycota_full"
O = os.path.dirname(os.path.abspath(__file__))
QUOTA = {"Sporidiobolales": 6, "Trichosporonales": 6, "Wallemiales": 5, "Microbotryales": 5,
         "Filobasidiales": 5, "Cystofilobasidiales": 5, "Cystobasidiales": 4, "Auriculariales": 4,
         "Cantharellales": 4, "Holtermanniales": 2, "Kriegeriales": 2, "Dacrymycetales": 2}
tax = {r["ASMID"]: r["NCBI_TAXONID"] for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"))}
rows = [r for r in csv.DictReader(open(B + "/genomes.tsv"), delimiter="\t")
        if r["routing"] == "phylum_fallback" and r["status"] == "uncalled"
        and int(r["capped_clusters"] or 0) > 0 and int(r["size_bp"] or 0) < 750e6]
by = collections.defaultdict(list)
for r in rows: by[r["order"]].append(r)
pick = []
for o, q in QUOTA.items():
    pick += sorted(by[o], key=lambda r: (-int(r["capped_clusters"]), r["genome"]))[:q]
with open(O + "/selected.tsv", "w") as fo:
    fo.write("genome\torder\tspecies\tsize_bp\twall_s_capped\tcapped_clusters\tadmitted_clusters\tchunk_capped\n")
    for r in pick:
        fo.write("\t".join([r["genome"], r["order"], r["species"], r["size_bp"], r["wall_s"], r["capped_clusters"], r["admitted_clusters"], r["chunk"]]) + "\n")
# (chunking below superseded: lists were re-split into 3 balanced chunks by greedy load)
pick.sort(key=lambda r: -float(r["wall_s"]))
chunks, load = [], []
for r in pick:
    est = 3 * float(r["wall_s"])
    for i in range(len(chunks)):
        if load[i] + est <= 1.2 * 3600 * 16 and max(3 * float(x["wall_s"]) for x in chunks[i] + [r]) < 1.6 * 3600:
            chunks[i].append(r); load[i] += est; break
    else:
        chunks.append([r]); load.append(est)
for i, ch in enumerate(chunks):
    with open(f"{O}/lists/capoff_{i:02d}.tsv", "w") as fo:
        for r in ch: fo.write(f"{r['genome']}\t{tax[r['genome']]}\n")
    print(f"capoff_{i:02d}: {len(ch)} genomes, est {load[i]/16/3600:.2f} h at 16-way, max single {max(3*float(x['wall_s']) for x in ch)/3600:.2f} h")
print("selected", len(pick))
