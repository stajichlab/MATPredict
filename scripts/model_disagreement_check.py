"""Splice-motif and model-quality check of genes whose two polish tools disagree on the model (#290, #349).

A gene with `status: polished_disagree` carries two models in the report: the primary (`method`, `exons`) and the `alternate_model`
(the other tool). For each such gene this script reads the two models' exons, takes the genome sequence, and reports for each model:

* the introns, with the donor and acceptor dinucleotides read on the gene's strand and a class (`GT-AG`, `GC-AG`, `AT-AC` are
  canonical; anything else is `noncanonical`; an intron under 30 bp is flagged `short`);
* the coding sequence length, whether it is a multiple of 3, the translation (genetic code from the report: 1, or 12 for the
  CTG-clade yeasts), internal stop codons, whether it starts with Met and ends in a stop, its length in aa and the ratio to the
  reference protein of the record the gene was matched to;
* which introns the two models share and which only one has.

A recommendation per gene follows from fixed rules (see `recommend`): a model is CLEAN when it has no internal stop, a length that is a
multiple of 3 and only canonical introns of at least 30 bp. One clean model -> that model; both clean -> both plausible (the one closer to
the reference length is named, but the choice needs RNA-seq or a curator); neither -> neither clean. When the model the motifs favour has an
identity more than 5 points lower AND is farther from the reference length than the other, the answer is `conflict` (see `reference_conflict`).
It is a recommendation, not a ruling:
a canonical splice site is necessary, not sufficient (a wrong junction can sit on another GT...AG), and an unusual intron can be real.

Usage: model_disagreement_check.py GENES.tsv RUNS_DIR GENOME_DIR GENOME_TEMPLATE DB_ROOT OUT
  GENES.tsv       columns `sample` and `gene` (and optionally `proposed_action`; rows whose action does not say "tools disagree" are
                  skipped when that column is present)
  RUNS_DIR        folder of `<sample>/detection_report.yaml`
  GENOME_TEMPLATE e.g. `FungiDB-71_{sample}_Genome.fasta.gz`, resolved in GENOME_DIR
"""
from __future__ import annotations

import csv
import gzip
import sys
from pathlib import Path

import yaml

CANONICAL = {"GT-AG", "GC-AG", "AT-AC"}
MIN_INTRON = 30
COMP = str.maketrans("ACGTacgtN", "TGCAtgcaN")
BASES = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: AA[16 * i + 4 * j + k] for i, a in enumerate(BASES) for j, b in enumerate(BASES) for k, c in enumerate(BASES)}


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def translate(cds: str, table: int = 1) -> str:
    codon = dict(CODON)
    if table == 12:          # alternative yeast nuclear: CTG is Ser
        codon["CTG"] = "S"
    elif table != 1:
        raise ValueError(f"genetic code {table} is not supported by this check (1 and 12 are)")
    return "".join(codon.get(cds[i:i + 3].upper(), "X") for i in range(0, len(cds) - len(cds) % 3, 3))


def read_contig(fa: Path, contig: str) -> str | None:
    op = gzip.open if str(fa).endswith(".gz") else open
    buf, on = [], False
    with op(fa, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if on:
                    break
                on = line[1:].split()[0] == contig
            elif on:
                buf.append(line.strip())
    return "".join(buf) or None


def intron_class(donor: str, acceptor: str) -> str:
    c = f"{donor.upper()}-{acceptor.upper()}"
    return c if c in CANONICAL else "noncanonical"


def analyse_model(seq: str, exons: list[tuple[int, int]], strand: str, table: int = 1, ref_aa: int | None = None) -> dict:
    """Introns, CDS and translation of one model (1-based inclusive exon coordinates on the contig, in contig order)."""
    exons = sorted(exons)
    introns = []
    for (s1, e1), (s2, e2) in zip(exons, exons[1:]):
        i0, i1 = e1 + 1, s2 - 1
        raw = seq[i0 - 1:i1]
        site = raw if strand == "+" else revcomp(raw)
        donor, acceptor = site[:2], site[-2:]
        introns.append({"start": i0, "end": i1, "length": i1 - i0 + 1, "donor": donor.upper(), "acceptor": acceptor.upper(),
                        "class": intron_class(donor, acceptor), "short": (i1 - i0 + 1) < MIN_INTRON})
    cds = "".join(seq[s - 1:e] for s, e in exons)
    if strand == "-":
        cds = revcomp(cds)
    prot = translate(cds, table)
    ends_stop = prot.endswith("*")
    body = prot[:-1] if ends_stop else prot
    out = {"exons": exons, "introns": introns, "cds_len": len(cds), "in_frame": len(cds) % 3 == 0, "protein": body,
           "aa": len(body), "internal_stops": body.count("*"), "starts_met": body.startswith("M"), "ends_stop": ends_stop,
           "ref_ratio": round(len(body) / ref_aa, 2) if ref_aa else None}
    return out


def issues(m: dict) -> list[str]:
    """Why a model is not CLEAN (empty list: clean)."""
    out = []
    if m["internal_stops"]:
        out.append(f'{m["internal_stops"]} internal stop(s)')
    if not m["in_frame"]:
        out.append("length not a multiple of 3")
    bad = [i for i in m["introns"] if i["class"] == "noncanonical"]
    if bad:
        out.append("noncanonical intron(s): " + ",".join(f'{i["start"]}-{i["end"]} {i["donor"]}..{i["acceptor"]}' for i in bad))
    short = [i for i in m["introns"] if i["short"]]
    if short:
        out.append("intron(s) under 30 bp: " + ",".join(f'{i["start"]}-{i["end"]}' for i in short))
    return out


def shared_introns(a: dict, b: dict) -> tuple[list, list, list]:
    ka = {(i["start"], i["end"]) for i in a["introns"]}
    kb = {(i["start"], i["end"]) for i in b["introns"]}
    return sorted(ka & kb), sorted(ka - kb), sorted(kb - ka)


IDENTITY_MARGIN = 5.0


def reference_conflict(chosen: dict, other: dict, chosen_identity, other_identity, ref_aa) -> bool:
    """True when the model the splice motifs favour is clearly worse against the reference: lower identity by more than
    IDENTITY_MARGIN points AND farther from the reference length. Identity comes from the tool that made the model (exonerate aligns
    against the reference, so it favours its own model); the length test is the less biased of the two, so both must point the same way."""
    if None in (chosen_identity, other_identity) or not ref_aa:
        return False
    return chosen_identity < other_identity - IDENTITY_MARGIN and abs(chosen["aa"] - ref_aa) > abs(other["aa"] - ref_aa)


def recommend(primary: dict, alt: dict, ref_aa: int | None = None, primary_identity=None, alt_identity=None) -> tuple[str, str]:
    """(recommendation, reason). A recommendation, not a ruling."""
    ip, ia = issues(primary), issues(alt)
    if not ip and ia:
        if reference_conflict(primary, alt, primary_identity, alt_identity, ref_aa):
            return "conflict", "motifs favour the primary (the alternate has: " + "; ".join(ia) + ") but the primary is clearly worse against the reference"
        return "primary model", "the alternate has: " + "; ".join(ia)
    if ip and not ia:
        if reference_conflict(alt, primary, alt_identity, primary_identity, ref_aa):
            return "conflict", ("motifs favour the alternate (the primary has: " + "; ".join(ip) + f") but identity ({primary_identity} against {alt_identity}) and length "
                                f"({primary['aa']} aa against {alt['aa']}, reference {ref_aa}) favour the primary")
        return "alternate model", "the primary has: " + "; ".join(ip)
    if ip and ia:
        return "neither clean", "primary: " + "; ".join(ip) + " | alternate: " + "; ".join(ia)
    if ref_aa:
        closer = "primary" if abs(primary["aa"] - ref_aa) <= abs(alt["aa"] - ref_aa) else "alternate"
        if abs(primary["aa"] - alt["aa"]) == 0:
            return "both plausible", f"same length ({primary['aa']} aa); they differ in intron position only"
        return "both plausible", f"{primary['aa']} aa vs {alt['aa']} aa (reference {ref_aa}); the {closer} is closer in length"
    return "both plausible", f"{primary['aa']} aa vs {alt['aa']} aa"


def model_exons(entry: dict) -> list[tuple[int, int]]:
    exs = entry.get("exons") or [{"start": entry["start"], "end": entry["end"]}]
    return [(e["start"], e["end"]) for e in exs]


def reference_aa(db_root: Path, record: str | None, gene: str) -> int | None:
    """Length of `gene`'s protein in the reference record's proteins.faa (None when not found)."""
    if not record:
        return None
    for faa in db_root.glob(f"*/*/{record}/proteins.faa"):
        name, n = None, 0
        for line in list(faa.read_text().splitlines()) + [">end"]:
            if line.startswith(">"):
                if name == gene:
                    return n
                name = next((f.split("=", 1)[1] for f in line[1:].split("|") if f.startswith("name=")), None)
                n = 0
            else:
                n += len(line.strip())
    return None


def check_gene(entry: dict, seq: str, table: int, ref_aa: int | None) -> dict:
    alt = entry["alternate_model"]
    p = analyse_model(seq, model_exons(entry), entry["strand"], table, ref_aa)
    a = analyse_model(seq, model_exons(alt), alt.get("strand", entry["strand"]), table, ref_aa)
    shared, only_p, only_a = shared_introns(p, a)
    rec, why = recommend(p, a, ref_aa, entry.get("identity"), alt.get("identity"))
    return {"primary": p, "alternate": a, "shared": shared, "only_primary": only_p, "only_alternate": only_a, "recommendation": rec, "reason": why}


def run(genes_tsv: Path, runs: Path, genome_dir: Path, template: str, db_root: Path, out: Path) -> list[dict]:
    wanted = []
    with genes_tsv.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if "proposed_action" in r and "tools disagree" not in r["proposed_action"]:
                continue
            wanted.append((r["sample"], r["gene"]))
    rows, intron_rows, faa = [], [], []
    for sample, gene in wanted:
        doc = yaml.safe_load((runs / sample / "detection_report.yaml").read_text())
        table = int(doc.get("genetic_code") or 1)
        genome = genome_dir / template.format(sample=sample)
        k = 0
        for locus in doc.get("detected") or []:
            for e in locus.get("gene_evidence") or []:
                if e["gene"] != gene or e.get("status") != "polished_disagree" or not e.get("alternate_model"):
                    continue
                k += 1
                seq = read_contig(genome, e["contig"])
                if seq is None:
                    rows.append({"sample": sample, "gene": gene, "contig": e["contig"], "recommendation": "not_checked", "reason": "contig not found in " + genome.name})
                    continue
                ref = reference_aa(db_root, e.get("reference_record"), gene)
                r = check_gene(e, seq, table, ref)
                p, a = r["primary"], r["alternate"]
                rows.append({"sample": sample, "gene": gene, "family": locus["family"], "tier": locus["confidence"], "contig": e["contig"],
                             "strand": e["strand"], "code": table, "reference": e.get("reference_record"), "ref_aa": ref,
                             "primary_tool": e["method"], "primary_aa": p["aa"], "primary_exons": len(p["exons"]), "primary_introns": ",".join(f'{i["class"]}' for i in p["introns"]),
                             "primary_stops": p["internal_stops"], "primary_met": p["starts_met"], "primary_ref_ratio": p["ref_ratio"], "primary_identity": e.get("identity"),
                             "alternate_tool": e["alternate_model"].get("method"), "alternate_aa": a["aa"], "alternate_exons": len(a["exons"]),
                             "alternate_introns": ",".join(f'{i["class"]}' for i in a["introns"]), "alternate_stops": a["internal_stops"], "alternate_met": a["starts_met"],
                             "alternate_ref_ratio": a["ref_ratio"], "alternate_identity": e["alternate_model"].get("identity"),
                             "shared_introns": len(r["shared"]), "only_primary_introns": len(r["only_primary"]), "only_alternate_introns": len(r["only_alternate"]),
                             "recommendation": r["recommendation"], "reason": r["reason"]})
                for which, m in (("primary", p), ("alternate", a)):
                    for i in m["introns"]:
                        intron_rows.append([sample, gene, k, which, e["contig"], i["start"], i["end"], i["length"], i["donor"], i["acceptor"], i["class"]])
                    faa.append(f">{sample}|{gene}|{which}|{e['contig']}\n{m['protein']}\n")
    out.mkdir(parents=True, exist_ok=True)
    cols = ["sample", "gene", "family", "tier", "contig", "strand", "code", "reference", "ref_aa", "primary_tool", "primary_aa", "primary_exons", "primary_introns",
            "primary_stops", "primary_met", "primary_ref_ratio", "primary_identity", "alternate_tool", "alternate_aa", "alternate_exons", "alternate_introns",
            "alternate_stops", "alternate_met", "alternate_ref_ratio", "alternate_identity", "shared_introns", "only_primary_introns", "only_alternate_introns",
            "recommendation", "reason"]
    with open(out / "model_check.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", restval="")
        w.writeheader()
        w.writerows(rows)
    with open(out / "introns.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "gene", "n", "model", "contig", "start", "end", "length", "donor", "acceptor", "class"])
        w.writerows(intron_rows)
    (out / "models.faa").write_text("".join(faa))
    return rows


def main(argv: list[str]) -> int:
    genes, runs, gdir, tmpl, db, out = argv
    rows = run(Path(genes), Path(runs), Path(gdir), tmpl, Path(db), Path(out))
    from collections import Counter
    print(dict(Counter(r["recommendation"] for r in rows)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
