#!/usr/bin/env python3
"""Harvest submitter-stated mating types from NCBI BioSample for a set of genome assemblies.

    harvest_biosample_mating_types.py --genomes genomes.tsv --out-dir DIR --email you@example.org

`--genomes` is a TSV with a `genome` column holding BFD ids that start with an assembly accession
(GCA_000150145.1_ASM15014v1). Steps: map each accession to its BioSample through the NCBI fungal assembly
summaries (GenBank and RefSeq), fetch each BioSample once with E-utilities, and keep every attribute whose name
or harmonised name says mating type (never `material`). Writes

  biosample_map.tsv         genome, accession, biosample, organism, strain, isolate
  biosample_mating.tsv      genome, biosample, organism, attribute, value   (one row per matching attribute)
  biosample_coverage.tsv    genome, biosample, n_attributes, has_mating_attribute

The values are the submitter's, free text, and unverified: they are a truth source independent of this pipeline, not
ground truth. Normalising them to MAT1-1 / MAT1-2 is a separate step (`normalise_mating_values`).
Standard library only. Respects the NCBI limit (3 requests per second without a key).
"""
import argparse
import csv
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

SUMMARIES = {
    "GCA": "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/fungi/assembly_summary.txt",
    "GCF": "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/assembly_summary.txt",
}
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
#: attribute names that mean mating type; `material` and similar are not matched
NAME_RE = re.compile(r"(?<![a-z])(mating|mat[-_ ]?type|mat[-_ ]?locus|mat1|matingtype)", re.I)
BATCH = 100


def fetch(url: str, params: dict | None = None, retries: int = 4) -> bytes:
    if params:
        url = url + "?" + urllib.parse.urlencode(params)
    for k in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return r.read()
        except Exception as exc:  # noqa: BLE001 - retry any transient network error
            if k == retries - 1:
                raise
            print(f"retry {k + 1} after {type(exc).__name__}", file=sys.stderr)
            time.sleep(3 * (k + 1))


def assembly_to_biosample(accessions: set[str], out: Path) -> dict[str, dict]:
    """accession (GCA_/GCF_ with version) -> {biosample, organism, strain, isolate}."""
    found = {}
    for prefix, url in SUMMARIES.items():
        wanted = {a for a in accessions if a.startswith(prefix)}
        if not wanted:
            continue
        path = out / f"assembly_summary_{prefix}.txt"
        if not path.exists():
            path.write_bytes(fetch(url))
        for line in path.read_text().splitlines():
            if line.startswith("#") or not line:
                continue
            f = line.split("\t")
            if f[0] in wanted:
                found[f[0]] = dict(biosample=f[2], organism=f[7], strain=f[8].replace("strain=", ""),
                                   isolate=f[9])
    return found


def parse_biosamples(xml: bytes) -> list[dict]:
    rows = []
    root = ET.fromstring(xml)
    for bs in root.iter("BioSample"):
        acc = bs.get("accession")
        attrs = [(a.get("attribute_name") or "", a.get("harmonized_name") or "", (a.text or "").strip())
                 for a in bs.iter("Attribute")]
        org = bs.find(".//Organism")
        rows.append(dict(biosample=acc, organism=(org.get("taxonomy_name") if org is not None else ""),
                         n_attributes=len(attrs), attrs=attrs))
    return rows


def mating_attributes(row: dict) -> list[tuple[str, str]]:
    return [(name or harm, value) for name, harm, value in row["attrs"]
            if (NAME_RE.search(name) or NAME_RE.search(harm)) and value]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--email", required=True)
    ap.add_argument("--limit", type=int, default=0, help="only the first N genomes (testing)")
    a = ap.parse_args(argv)
    a.out_dir.mkdir(parents=True, exist_ok=True)
    genomes = [r["genome"] for r in csv.DictReader(open(a.genomes), delimiter="\t")]
    if a.limit:
        genomes = genomes[:a.limit]
    acc_of = {g: "_".join(g.split("_")[:2]) for g in genomes}
    amap = assembly_to_biosample(set(acc_of.values()), a.out_dir)
    print(f"{len(genomes)} genomes, {len(amap)} with a BioSample in the summaries", file=sys.stderr)
    with open(a.out_dir / "biosample_map.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["genome", "accession", "biosample", "organism", "strain", "isolate"])
        for g in genomes:
            m = amap.get(acc_of[g])
            if m:
                w.writerow([g, acc_of[g], m["biosample"], m["organism"], m["strain"], m["isolate"]])
    biosamples = sorted({m["biosample"] for m in amap.values() if m["biosample"].startswith("SAM")})
    fetched = {}
    for i in range(0, len(biosamples), BATCH):
        batch = biosamples[i:i + BATCH]
        xml = fetch(EFETCH, {"db": "biosample", "id": ",".join(batch), "retmode": "xml",
                             "tool": "matpredict-truth", "email": a.email})
        for row in parse_biosamples(xml):
            fetched[row["biosample"]] = row
        time.sleep(0.4)
        if (i // BATCH) % 20 == 0:
            print(f"  fetched {min(i + BATCH, len(biosamples))}/{len(biosamples)}", file=sys.stderr)
    by_bs = defaultdict_list = {}
    with open(a.out_dir / "biosample_mating.tsv", "w", newline="") as fm, \
            open(a.out_dir / "biosample_coverage.tsv", "w", newline="") as fc:
        wm, wc = csv.writer(fm, delimiter="\t"), csv.writer(fc, delimiter="\t")
        wm.writerow(["genome", "biosample", "organism", "attribute", "value"])
        wc.writerow(["genome", "biosample", "n_attributes", "has_mating_attribute"])
        for g in genomes:
            m = amap.get(acc_of[g])
            row = fetched.get(m["biosample"]) if m else None
            if not row:
                continue
            ma = mating_attributes(row)
            wc.writerow([g, row["biosample"], row["n_attributes"], bool(ma)])
            for name, value in ma:
                wm.writerow([g, row["biosample"], row["organism"], name, value])
    print(f"wrote {a.out_dir}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
