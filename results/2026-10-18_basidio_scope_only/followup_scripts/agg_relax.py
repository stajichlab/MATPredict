import sys, glob, yaml, csv, collections
sys.path.insert(0, "/bigdata/stajichlab/jstajich/projects/mat-fam18/scripts")
import family_performance_table as f
W = "/bigdata/stajichlab/jstajich/projects/mat-fam18/"
grp = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in open(W + "rc/groups.tsv")}
n50 = {r["ASM_FOLDER"]: int(r["N50"] or 0) for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/ncbi_accessions.csv"))}


def loci(p):
    d = yaml.load(open(p), Loader=yaml.CSafeLoader)
    return [x for x in d.get("detected") or [] if x["family"].endswith("rustHD")], d


rows = []; newl = []
for p in sorted(glob.glob(W + "runs_relax/bar1/*/detection_report.yaml")):
    g = p.split("/")[-2]; acc = g.split("_", 1)[1]
    base_p = glob.glob(W + f"runs/pucc/*_{acc}/detection_report.yaml")[0]
    b, _ = loci(base_p); r, d = loci(p)
    def ov(x, ys): return any(x["contig"] == y["contig"] and x["start"] <= y["end"] and x["end"] >= y["start"] for y in ys)
    new = [x for x in r if not ov(x, b)]
    for x in new:
        kind, gap = f.pair_arrangement(x, ("bE", "bW"))
        newl.append(dict(acc=acc, genus=grp[acc][1], group=grp[acc][3], n50=n50.get(acc, 0), contig=x["contig"], start=x["start"], end=x["end"], span=x["end"] - x["start"] + 1,
                         polished=x["polished_genes"], conf=x["confidence"], cls=x["locus_class"], genes="+".join(x["genes_found"]), arrangement=kind, gap=gap))
    rows.append(dict(acc=acc, genus=grp[acc][1], group=grp[acc][3], n50=n50.get(acc, 0), base_called=bool(b), base_loci=len(b), relaxed_loci=len(r), new_loci=len(new)))
with open(W + "rc/relax_genomes.tsv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
with open(W + "rc/relax_new_loci.tsv", "w") as fh:
    if newl:
        w = csv.DictWriter(fh, fieldnames=list(newl[0]), delimiter="\t"); w.writeheader(); w.writerows(newl)
print("genomes", len(rows))
unc = [r for r in rows if not r["base_called"]]
print("baseline-uncalled genomes in this set:", len(unc), "gained a call with bar 1:", sum(r["relaxed_loci"] > 0 for r in unc))
print("called genomes with extra loci:", sum(1 for r in rows if r["base_called"] and r["new_loci"]), "extra loci", sum(r["new_loci"] for r in rows if r["base_called"]))
print("new loci total", len(newl), "by polished", collections.Counter(x["polished"] for x in newl), "arrangement", collections.Counter(x["arrangement"] for x in newl),
      "classes", collections.Counter(x["cls"] for x in newl), "conf", collections.Counter(x["conf"] for x in newl))
print("new loci with BOTH genes on one contig and divergent:", sum(1 for x in newl if x["arrangement"] == "divergent"), "gaps", sorted(x["gap"] for x in newl if x["gap"] is not None)[:20])
for x in newl: print(x)
