"""Step 1: list every detected or withheld locus whose genes_found includes sexM or sexP.

Scans the early-diverging detection reports. Adds taxonomy from BFD samples.csv.
Skips ASMIDs on the BFD suppress list. Writes loci.tsv.
"""
import csv, glob, os, yaml
R = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-26_early_diverging"
SAMPLES = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"
SUPP = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/data/curation/suppress.txt"
supp = {l.split(",")[0].strip() for l in open(SUPP) if l.strip() and not l.startswith(("#", "ASMID"))}
tax = {r["ASMID"]: r for r in csv.DictReader(open(SAMPLES))}
rows = []
for group in ("Mucoromycota", "Mortierellomycota", "Kickxellomycota"):
    for rep in sorted(glob.glob(f"{R}/{group}/runs/*/detection_report.yaml")):
        g = os.path.basename(os.path.dirname(rep))
        if g in supp:
            continue
        r = yaml.safe_load(open(rep)) or {}
        t = tax.get(g, {})
        for status, key in (("called", "detected"), ("withheld", "suppressed_loci")):
            for i, x in enumerate(r.get(key) or []):
                genes = x.get("genes_found") or []
                if not ({"sexM", "sexP"} & set(genes)):
                    continue
                rows.append(dict(
                    locus_id=f"{g}|{status[0]}{i}", genome=g, group=group,
                    phylum=t.get("PHYLUM", ""), order=t.get("ORDER", ""),
                    family=t.get("FAMILY", ""), species=t.get("SPECIES", ""),
                    genetic_code=r.get("genetic_code", 1), status=status,
                    contig=x["contig"], start=x["start"], end=x["end"],
                    idiomorph=x.get("idiomorph", ""),
                    confidence=x.get("confidence", ""), locus_class=x.get("locus_class", ""),
                    genes_found=",".join(genes)))
with open("loci.tsv", "w", newline="") as fo:
    w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader(); w.writerows(rows)
import collections
c = collections.Counter((r["group"], r["status"]) for r in rows)
print(len(rows), "loci"); [print(k, v) for k, v in sorted(c.items())]
print("genetic codes:", collections.Counter(r["genetic_code"] for r in rows))
print("genomes:", len({r["genome"] for r in rows}))
