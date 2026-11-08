#!/usr/bin/env python3
"""Manifest of one species' samples for the sketch store: one row per public BioSample (one run each) and one per local
Population_Genomics sample, with the subsample size that gives COVERAGE x of the genome.

    build_manifest.py --runs fungi_runs_ena_<date>.tsv.gz --taxmap taxonomy_species_map.csv.gz --species-taxid 746128 \
        --genome-mb 29.4 --local-samples .../Afumigatus_Global/samples.csv --local-input .../Afumigatus_Global/input \
        --out manifest.tsv

Public: PopSeq's kingdom run table (ENA read_run, Illumina WGS), every taxon that rolls up to the species. Per BioSample the run kept
is the paired one with the most bases (single-end only if no paired run). Studies flagged by survey_ena.classify_study (crosses,
mutants, evolution lines, not WGS) are left out unless --include-flagged; the flag is kept in the row either way.
Local: samples.csv of a Population_Genomics subproject (strain, ';'-separated FASTQ patterns relative to input/, platform);
rows whose files are not FASTQ (Pre_assembled) are left out and listed on stderr.
sample_id: the BioSample accession for public rows, 'local__<strain>' for local rows. A local strain whose name matches a public
strain/isolate name is kept (it is a different sequencing of a possibly identical strain) and named in also_public."""
import argparse, csv, gzip, os, re, sys, collections

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "sra_survey"))
from survey_ena import classify_study  # noqa: E402

csv.field_size_limit(10 ** 9)
COLS = ["sample_id", "source", "run_accession", "biosample", "strain", "study_accession", "study_flag", "layout", "read_count",
        "base_count", "spot_bases", "target_spots", "full_coverage", "first_public", "country", "center_name", "local_pattern",
        "also_public", "ena_fastq"]


def norm(s):
    return re.sub(r"[^a-z0-9]", "", (s or "").lower())


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--runs", required=True); ap.add_argument("--taxmap", required=True); ap.add_argument("--species-taxid", required=True)
    ap.add_argument("--genome-mb", type=float, required=True); ap.add_argument("--coverage", type=float, default=30)
    ap.add_argument("--local-samples"); ap.add_argument("--local-input"); ap.add_argument("--include-flagged", action="store_true")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    target_bases = a.genome_mb * 1e6 * a.coverage
    tax = {r["taxid"]: r["species_taxid"] for r in csv.DictReader(gzip.open(a.taxmap, "rt"))}
    runs = [r for r in csv.DictReader(gzip.open(a.runs, "rt"), delimiter="\t") if tax.get(r["tax_id"], r["tax_id"]) == a.species_taxid]
    bs = lambda r: r["sample_accession"].split(";")[0] or r["run_accession"]
    by_study = collections.defaultdict(list)
    for r in runs:
        by_study[r["study_accession"]].append(r)
    flag = {}
    for st, rs in by_study.items():
        texts = {((r.get("strain") or r.get("isolate") or bs(r)).lower()): " ".join(r.get(f, "") for f in ("sample_title", "strain", "isolate")) for r in rs}
        flag[st] = classify_study(rs[0]["study_title"], texts)
    best = {}
    for r in runs:
        if not int(r["read_count"] or 0):
            continue
        key = (r["library_layout"] == "PAIRED", int(r["base_count"] or 0))
        b = bs(r)
        if b not in best or key > best[b][0]:
            best[b] = (key, r)
    rows, public_names = [], {}
    for b, (_, r) in sorted(best.items()):
        f = flag[r["study_accession"]]
        if f and not a.include_flagged:
            continue
        reads, bases = int(r["read_count"]), int(r["base_count"] or 0)
        spot = bases / reads if reads else 0
        tgt = int(min(reads, -(-target_bases // spot))) if spot else 0
        name = r.get("strain") or r.get("isolate") or ""
        if name:
            public_names.setdefault(norm(name), b)
        rows.append(dict(sample_id=b, source="public", run_accession=r["run_accession"], biosample=b, strain=name,
                         study_accession=r["study_accession"], study_flag=f, layout=r["library_layout"], read_count=reads, base_count=bases,
                         spot_bases=f"{spot:.1f}", target_spots=tgt, full_coverage=f"{bases / (a.genome_mb * 1e6):.1f}",
                         first_public=r["first_public"], country=r["country"], center_name=r["center_name"], local_pattern="", also_public="",
                         ena_fastq=r.get("fastq_ftp", "")))
    n_pub = len(rows)
    skipped = []
    if a.local_samples:
        for rec in csv.reader(open(a.local_samples)):
            if not rec or not rec[0].strip():
                continue
            strain, pat = rec[0].strip(), rec[1].strip()
            if pat.startswith("Pre_assembled") or not re.search(r"\.f(ast)?q(\.gz)?$", pat.split(";")[0]):
                skipped.append(strain); continue
            rows.append(dict(sample_id="local__" + re.sub(r"[^A-Za-z0-9._-]", "_", strain), source="local", run_accession="", biosample="",
                             strain=strain, study_accession=pat.split("/")[0], study_flag="", layout="PAIRED" if "[12]" in pat or "R?" in pat else "",
                             read_count="", base_count="", spot_bases="", target_spots="", full_coverage="", first_public="", country="",
                             center_name="Population_Genomics", local_pattern=pat, also_public=public_names.get(norm(strain), ""), ena_fastq=""))
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, COLS, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
    fl = collections.Counter(flag[r["study_accession"]] for _, r in best.values())
    print(f"public: {len(best)} BioSamples with reads ({dict(fl)}); kept {n_pub}; local: {len(rows) - n_pub} "
          f"({sum(1 for r in rows if r['also_public'])} name-match a public strain); local skipped (not FASTQ): {len(skipped)} {skipped[:12]}", file=sys.stderr)


if __name__ == "__main__":
    main()
