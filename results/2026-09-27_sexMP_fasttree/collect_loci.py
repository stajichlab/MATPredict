"""Step 1: sexM/sexP loci for the FastTree redraw.

Mucoromycota loci come from the final classifier scan (Mucoromycota_1a00b0a).
Mortierellomycota and Kickxellomycota loci are carried over unchanged from
results/2026-09-26_sexMP_phylogeny/loci.tsv. Each new Mucoromycota locus is
matched to an old one by (genome, status, contig, start, end); a match keeps
the old locus_id so its extracted proteins can be reused. Writes loci.tsv and
changed_genomes.txt (Mucoromycota genomes whose locus set changed at all).
"""
import collections, csv, glob, os, yaml

NEW = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-26_hmm_classifier/Mucoromycota_1a00b0a/runs"
OLD = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-26_sexMP_phylogeny/loci.tsv"
SAMPLES = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"
SUPP = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/data/curation/suppress.txt"
supp = {l.split(",")[0].strip() for l in open(SUPP) if l.strip() and not l.startswith(("#", "ASMID"))}
tax = {r["ASMID"]: r for r in csv.DictReader(open(SAMPLES))}

old = list(csv.DictReader(open(OLD), delimiter="\t"))
old_key = {(r["genome"], r["status"], r["contig"], r["start"], r["end"]): r["locus_id"] for r in old}
old_by_g = collections.defaultdict(set)
for r in old:
    if r["group"] == "Mucoromycota":
        old_by_g[r["genome"]].add((r["status"], r["contig"], r["start"], r["end"]))

rows = [r for r in old if r["group"] != "Mucoromycota"]
new_by_g = collections.defaultdict(set)
for rep in sorted(glob.glob(f"{NEW}/*/detection_report.yaml")):
    g = os.path.basename(os.path.dirname(rep))
    if g in supp:
        continue
    r = yaml.safe_load(open(rep)) or {}
    t = tax.get(g, {})
    n = collections.Counter()
    for status, key in (("called", "detected"), ("withheld", "suppressed_loci")):
        for x in r.get(key) or []:
            genes = x.get("genes_found") or []
            if not ({"sexM", "sexP"} & set(genes)):
                continue
            k = (g, status, x["contig"], str(x["start"]), str(x["end"]))
            lid = old_key.get(k)
            if lid is None:
                lid = f"{g}|n{status[0]}{n[status]}"
            n[status] += 1
            new_by_g[g].add(k[1:])
            rows.append(dict(
                locus_id=lid, genome=g, group="Mucoromycota",
                phylum=t.get("PHYLUM", ""), order=t.get("ORDER", ""),
                family=t.get("FAMILY", ""), species=t.get("SPECIES", ""),
                genetic_code=r.get("genetic_code", 1), status=status,
                contig=x["contig"], start=x["start"], end=x["end"],
                idiomorph=x.get("idiomorph", ""),
                confidence=x.get("confidence", ""), locus_class=x.get("locus_class", ""),
                genes_found=",".join(genes)))

changed = sorted(g for g in set(new_by_g) | set(old_by_g) if new_by_g[g] != old_by_g[g])
with open("loci.tsv", "w", newline="") as fo:
    w = csv.DictWriter(fo, fieldnames=list(old[0]), delimiter="\t")
    w.writeheader(); w.writerows(rows)
open("changed_genomes.txt", "w").write("\n".join(changed) + "\n")
c = collections.Counter((r["group"], r["status"]) for r in rows)
print(len(rows), "loci"); [print(" ", k, v) for k, v in sorted(c.items())]
print("Mucoromycota genomes with a changed locus set:", len(changed))
print("new locus ids:", sum(1 for r in rows if "|n" in r["locus_id"]))
