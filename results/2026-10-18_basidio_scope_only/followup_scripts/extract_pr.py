"""Per-locus redPR call table: label, identity/bitscore of the best STE3a1 and STE3a2 hit, nearest reference record."""
import csv, sys, yaml
from pathlib import Path
rundir = Path(sys.argv[1]); out = sys.argv[2]
B = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/ncbi_accessions_taxonomy.csv"
sp = {r["ASM_ACCESSION"]: (r["GENUS"], r["SPECIES"], r["NCBI_TAXID"]) for r in csv.DictReader(open(B))}
cols = ["genome", "taxid", "genus", "species", "n_loci", "contig", "start", "end", "label", "conf", "id_a1", "cov_a1", "bits_a1", "ref_a1",
        "id_a2", "cov_a2", "bits_a2", "ref_a2", "bits_margin", "id_margin", "ribl6", "resolution_basis"]
with open(out, "w") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(cols)
    for p in sorted(rundir.glob("*/detection_report.yaml")):
        g = p.parent.name; acc = g.split("_", 1)[1]
        genus, species, taxid = sp[acc]
        d = yaml.load(p.read_text(), Loader=yaml.CSafeLoader)
        loci = [x for x in d.get("detected") or [] if x["family"].endswith("redPR")]
        if not loci:
            w.writerow([g, taxid, genus, species, 0] + [""] * (len(cols) - 5)); continue
        for x in loci:
            best = {}
            for e in x["gene_evidence"]:
                if e["gene"] in ("STE3a1", "STE3a2") and (e.get("bitscore") or 0) >= (best.get(e["gene"], {}).get("bitscore") or 0):
                    best[e["gene"]] = e
            a1, a2 = best.get("STE3a1", {}), best.get("STE3a2", {})
            b1, b2 = a1.get("bitscore") or 0, a2.get("bitscore") or 0
            i1, i2 = a1.get("identity") or 0, a2.get("identity") or 0
            res = x.get("idiomorph_resolutions") or []
            w.writerow([g, taxid, genus, species, len(loci), x["contig"], x["start"], x["end"], x["idiomorph"], x["confidence"],
                        a1.get("identity"), a1.get("coverage"), a1.get("bitscore"), a1.get("reference_record"),
                        a2.get("identity"), a2.get("coverage"), a2.get("bitscore"), a2.get("reference_record"),
                        round(abs(b1 - b2), 1), round(abs(i1 - i2), 1), any(e["gene"] == "RibL6" for e in x["gene_evidence"]),
                        ",".join(sorted({r.get("basis") or "" for r in res}))])
