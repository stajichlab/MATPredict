#!/usr/bin/env python3
"""Survey of public fungal population read sets for `matpredict reads-type` (issue: see datasets/reads/sra_survey/README.md).

Counts, per species, how many distinct strains have paired-end Illumina whole-genome reads in ENA (which mirrors SRA), pooled
across studies, after leaving out studies that are mainly not natural isolates (mutant collections, experimental evolution
and mutation accumulation lines, crosses and their progeny). Reads are never downloaded here.

Inputs (made with the ENA portal API, see README.md): --runs (one row per run: run_accession, study_accession, sample_accession,
secondary_sample_accession, scientific_name, tax_id, library_layout, base_count, read_count, center_name, collection_date) and --meta
(run_accession, study_title, sample_title, sample_alias, strain, isolate, experiment_title, first_public, country, host).
Writes species_summary.tsv and studies.tsv next to --out. Deterministic: same inputs, same bytes.
"""
import argparse, csv, collections, os, re

MUTANT = re.compile(r"mutant|knock-?out|deletion|disrupt|insertion library|transposon|t-dna|crispr|overexpress|rnai|gene replacement|"
                    r"deleted|Δ|delta[- ]?\w+|\bko\b|library of (mutants|strains)", re.I)
EVOLVE = re.compile(r"experimental evolution|evolution experiment|evolved|serial passage|mutation[- ]accumulation|passaged|adaptive evolution|"
                    r"laboratory evolution|time[- ]?course|selection experiment", re.I)
CROSS = re.compile(r"progeny|segregant|backcross|\bf[12]\b|tetrad|meiotic|recombinant (inbred|strains)|genetic cross|offspring", re.I)
NOT_SPECIES = re.compile(r"metagenom|uncultured|environmental|unidentified|mixed|fungal sp|fungus sp", re.I)
MIN_BASES = 300_000_000          # about 10x of a 30 Mb genome; a per-species genome size would be better (open item)


def species(name):
    t = name.replace("[", "").replace("]", "").split()
    if len(t) >= 2 and re.fullmatch(r"[A-Z][a-z-]+", t[0]) and re.fullmatch(r"[a-z-]+\.?", t[1]) and t[1] not in ("sp.", "cf.", "aff."):
        return f"{t[0]} {t[1]}"
    return None


def db_tiers(db):
    """({species: set(idiomorphs)}, {genus: set(idiomorphs)}) from the accepted records' metadata.yaml (candidates are left out)."""
    import glob
    import yaml
    sp, gen = collections.defaultdict(set), collections.defaultdict(set)
    for f in sorted(glob.glob(os.path.join(db, "*", "*", "*", "metadata.yaml"))):
        d = yaml.safe_load(open(f)) or {}
        name = (d.get("organism") or {}).get("species", "")
        ids = (d.get("mating_type") or {}).get("idiomorphs") or []
        if name:
            sp[name].update(ids)
            gen[name.split()[0]].update(ids)
    return sp, gen


def tier(name, sp, gen):
    """A: the species has records for two or more idiomorphs; A1: one idiomorph; B: only another species of the genus has records; C: none.
    Name matching is by the scientific name as deposited, so a synonym (Pyricularia oryzae / Magnaporthe oryzae) reads as B or C."""
    if name in sp:
        return "A" if len(sp[name]) >= 2 else "A1"
    return "B" if name.split()[0] in gen else "C"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--runs", required=True)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--out", required=True, help="output folder")
    ap.add_argument("--min-samples", type=int, default=4)
    ap.add_argument("--db", help="MATPredict db folder: adds the reference tier of each species")
    a = ap.parse_args()
    meta = {}
    with open(a.meta, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            meta[r["run_accession"]] = r
    strains = collections.defaultdict(lambda: collections.defaultdict(dict))     # species -> study -> strain key -> (run, bases)
    title = {}
    n_runs = n_used = 0
    with open(a.runs, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            n_runs += 1
            sp = species(r["scientific_name"])
            if not sp or NOT_SPECIES.search(r["scientific_name"]):
                continue
            try:
                bases = int(r["base_count"] or 0)
            except ValueError:
                bases = 0
            if r["library_layout"] != "PAIRED" or bases < MIN_BASES:
                continue
            m = meta.get(r["run_accession"], {})
            key = (m.get("strain") or m.get("isolate") or m.get("sample_alias") or r["sample_accession"] or r["run_accession"]).strip().lower()
            st = r["study_accession"]
            title.setdefault(st, m.get("study_title", ""))
            prev = strains[sp][st].get(key)
            if prev is None or bases > prev[1]:
                strains[sp][st][key] = (r["run_accession"], bases, " ".join(m.get(k, "") for k in ("sample_title", "strain", "isolate", "experiment_title")))
            n_used += 1
    os.makedirs(a.out, exist_ok=True)
    spdb, gendb = db_tiers(a.db) if a.db else ({}, {})
    study_rows, species_rows = [], []
    for sp in sorted(strains):
        pooled, excluded = {}, collections.Counter()
        for st in sorted(strains[sp]):
            ks = strains[sp][st]
            txt = {k: v[2] for k, v in ks.items()}
            flag = ""
            tt = title.get(st, "")
            frac = lambda rx: sum(1 for k in ks if rx.search(txt[k]) or rx.search(k)) / len(ks)
            if MUTANT.search(tt) or frac(MUTANT) >= .5:
                flag = "mutant"
            elif EVOLVE.search(tt) or frac(EVOLVE) >= .5:
                flag = "evolution"
            elif CROSS.search(tt) or frac(CROSS) >= .5:
                flag = "cross"
            study_rows.append([sp, st, len(ks), flag, tt[:160]])
            if flag:
                excluded[flag] += len(ks)
                continue
            for k, v in ks.items():
                pooled.setdefault(k, (st, v[0]))
        species_rows.append([sp, len(pooled), len({v[0] for v in pooled.values()}), sum(len(v) for v in strains[sp].values()),
                             excluded["mutant"], excluded["evolution"], excluded["cross"]])
    with open(os.path.join(a.out, "studies.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["species", "study", "strains", "excluded_as", "study_title"])
        w.writerows(r for r in study_rows if r[2] >= 2)
    species_rows.sort(key=lambda r: (-r[1], r[0]))
    with open(os.path.join(a.out, "species_summary.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["species", "strains_pooled", "studies", "strains_all_studies", "excluded_mutant", "excluded_evolution", "excluded_cross", "reference_tier"])
        w.writerows(r + [tier(r[0], spdb, gendb) if a.db else ""] for r in species_rows if r[1] >= a.min_samples or r[3] >= a.min_samples)
    keep = [r for r in species_rows if r[1] >= a.min_samples]
    print(f"{n_runs} runs; {n_used} paired runs at >= {MIN_BASES // 10**6} Mb with a species name; species with >= {a.min_samples} pooled strains: "
          f"{len(keep)} ({sum(r[1] for r in keep)} strains)")
    if a.db:
        for t in ("A", "A1", "B", "C"):
            sel = [r for r in keep if tier(r[0], spdb, gendb) == t]
            print(f"tier {t}: {len(sel)} species, {sum(r[1] for r in sel)} strains")
    print("species with >= 4 strains only before the exclusions:", sum(1 for r in species_rows if r[3] >= a.min_samples), "; exclusions (strains):",
          {k: sum(r[i] for r in species_rows) for k, i in (("mutant", 4), ("evolution", 5), ("cross", 6))})


if __name__ == "__main__":
    main()
