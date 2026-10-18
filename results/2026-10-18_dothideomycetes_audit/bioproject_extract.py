"""BioProject (and strain) per Dothideomycetes genome from the BFD samples.csv, for the species-ratio check.
Run on the HPCC (reads /bigdata/.../Fungi_BFD/samples.csv): python bioproject_extract.py genome_list.tsv out.tsv"""
import csv, sys
want = {r["ASM_ACCESSION"] for r in csv.DictReader(open(sys.argv[1]), delimiter="\t")}
with open(sys.argv[2], "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(["genome", "bioproject", "strain"])
    for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv")):
        if r["ASMID"] in want:
            w.writerow([r["ASMID"], r["BIOPROJECT"], r["STRAIN"]])
