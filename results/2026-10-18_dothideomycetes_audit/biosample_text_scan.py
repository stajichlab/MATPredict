"""Scan every attribute (name and value) plus title/description of the BioSample of each Dothideomycetes genome for
mating-type statements, beyond the `mating_type` attribute the 2026-10-17 harvest used.

    python biosample_text_scan.py --genomes genome_list.tsv --out-dir DIR --email you@example.org
Writes biosample_dothideo_attrs.tsv (one row per biosample with a hit: biosample, genome, organism, attribute, value)
and biosample_dothideo_coverage.tsv (genome, biosample, fetched, n_attributes, hit). Submitter text is not ground truth.
"""
import argparse, csv, re, sys, time, xml.etree.ElementTree as ET
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from harvest_biosample_mating_types import EFETCH, assembly_to_biosample, fetch

PAT = re.compile(r"(mating|mat[-_ ]?type|mat[-_ ]?locus|\bmat\s?[-_]?[12]\b|mat1|homothall|heterothall|idiomorph|"
                 r"self[- ]fertile|self[- ]sterile|sexual)", re.I)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genomes", required=True); ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--email", required=True)
    a = ap.parse_args()
    a.out_dir.mkdir(parents=True, exist_ok=True)
    genomes = [r["ASM_ACCESSION"] for r in csv.DictReader(open(a.genomes), delimiter="\t")]
    acc = {g: "_".join(g.split("_")[:2]) for g in genomes}
    amap = assembly_to_biosample(set(acc.values()), a.out_dir)
    bss = sorted({m["biosample"] for m in amap.values() if m["biosample"].startswith("SAM")})
    docs = {}
    for i in range(0, len(bss), 100):
        xml = fetch(EFETCH, {"db": "biosample", "id": ",".join(bss[i:i+100]), "retmode": "xml",
                             "tool": "matpredict-truth", "email": a.email})
        for bs in ET.fromstring(xml).iter("BioSample"):
            items = [(x.get("attribute_name") or "", (x.text or "").strip()) for x in bs.iter("Attribute")]
            t = bs.find(".//Title"); items.append(("title", (t.text or "") if t is not None else ""))
            for p in bs.iter("Paragraph"):
                items.append(("description", (p.text or "")))
            docs[bs.get("accession")] = items
        time.sleep(0.4)
    with open(a.out_dir / "biosample_dothideo_attrs.tsv", "w", newline="") as fh, \
         open(a.out_dir / "biosample_dothideo_coverage.tsv", "w", newline="") as fc:
        w = csv.writer(fh, delimiter="\t"); c = csv.writer(fc, delimiter="\t")
        w.writerow(["genome", "biosample", "organism", "strain", "attribute", "value"])
        c.writerow(["genome", "biosample", "fetched", "n_attributes", "hit"])
        for g in genomes:
            m = amap.get(acc[g])
            items = docs.get(m["biosample"]) if m else None
            if items is None:
                c.writerow([g, m["biosample"] if m else "", False, 0, False]); continue
            hits = [(n, v) for n, v in items if v and v.lower() not in ("missing", "not applicable", "na", "not collected")
                    and (PAT.search(n) or PAT.search(v))]
            c.writerow([g, m["biosample"], True, len(items), bool(hits)])
            for n, v in hits:
                w.writerow([g, m["biosample"], m["organism"], m["strain"], n, v])
main()
