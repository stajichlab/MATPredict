#!/usr/bin/env python
"""Score the alpha-box and HMG-box profiles on Ascomycota loci from a detection campaign.

    score_core_profiles.py --profiles db/Ascomycota/classifiers/MAT_core --reports DIR --genomes DIR \
        --lists ambiguous_loci.tsv label_disagrees.tsv --out scored.tsv [--controls 200]

For each target locus it translates the MAT1-1-1 and MAT1-2-1 gene-evidence entries (exons when modelled; the
hit span in all three frames when not), scores each protein on both profiles, and calls the locus by the profile:
MAT1-1 if the alpha-box entry's alpha score beats the HMG-box entry's HMG score by at least MIN_MARGIN bits,
MAT1-2 for the reverse, otherwise `undetermined`. A locus with one core gene is called by which profile that
protein scores higher on. Controls are random loci with a clear call, to show what the profile does when the
campaign label is not in doubt. The campaign label is never used to score.
"""
import argparse
import csv
import glob
import gzip
import os
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

import yaml
from Bio import SeqIO
from Bio.Seq import Seq

from MATPredict.detect.classifier_build import score

MIN_MARGIN = 25.0       # provisional: the Mucoromycota classifier's min_margin
GENES = ("MAT1-1-1", "MAT1-2-1")
L = getattr(yaml, "CSafeLoader", yaml.SafeLoader)


def load_hmms(d: Path):
    import pyhmmer
    out = {}
    for gene in GENES:
        with pyhmmer.plan7.HMMFile(str(d / f"{gene}.hmm")) as fh:
            out[gene] = fh.read()
    return out


_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def coding_sequence(exons, strand: str, contig_seq: str) -> str:
    """The spliced coding sequence of a model, reverse-complemented on the minus strand."""
    dna = "".join(contig_seq[a - 1:b] for a, b in sorted(exons))
    return dna.translate(_COMP)[::-1] if strand == "-" else dna


def modelled(g): return g["status"] in ("polished_agree", "polished_single", "polished_disagree")


def read_reports(root: str):
    """(genome, contig, start) -> {label, confidence, locus_class, entries{gene: entry}} for Ascomycota:MAT loci."""
    out = {}
    for path in glob.glob(root + "/**/detection_report.yaml", recursive=True):
        genome = os.path.basename(os.path.dirname(path))
        doc = yaml.load(open(path), Loader=L)
        for c in doc.get("detected") or []:
            if c.get("family") != "Ascomycota:MAT":
                continue
            ev = {}
            for g in c.get("gene_evidence") or []:
                if g["gene"] in GENES and not g.get("lost_to"):
                    if g["gene"] not in ev or (g.get("bitscore") or 0) > (ev[g["gene"]].get("bitscore") or 0):
                        ev[g["gene"]] = g
            out[(genome, c["contig"], int(c["start"]))] = dict(
                label=c["idiomorph"], confidence=c["confidence"], locus_class=c["locus_class"], entries=ev)
    return out


def proteins(entry, seq: str) -> list[str]:
    """Modelled: the spliced CDS translation. Otherwise the span translated in the three frames of its strand."""
    strand = entry["strand"]
    if entry.get("exons"):
        ex = tuple(sorted((int(x["start"]), int(x["end"])) for x in entry["exons"]))
        dna = coding_sequence(ex, strand, seq)
        return [str(Seq(dna[: len(dna) // 3 * 3]).translate()).replace("*", "")]
    dna = coding_sequence(((int(entry["start"]), int(entry["end"])),), strand, seq)
    return [str(Seq(dna[f:][: (len(dna) - f) // 3 * 3]).translate()).replace("*", "X") for f in range(3)]


def best_scores(hmms, prots: list[str]) -> dict[str, float]:
    res = {}
    for gene, h in hmms.items():
        sc = score(h, {f"p{i}": p for i, p in enumerate(prots) if p})
        res[gene] = max(sc.values(), default=0.0)
    return res


def genome_contigs(genomes: str, genome: str, wanted: set[str]) -> dict[str, str]:
    for suffix in (".fa.gz", ".fa", ".fna.gz", ".fna"):
        p = Path(genomes) / (genome + suffix)
        if p.exists():
            opener = gzip.open if p.suffix == ".gz" else open
            with opener(p, "rt") as fh:
                return {r.id: str(r.seq).upper() for r in SeqIO.parse(fh, "fasta") if r.id in wanted}
    return {}


def call(a_own, h_own, margin=MIN_MARGIN):
    if a_own - h_own >= margin: return "MAT1-1"
    if h_own - a_own >= margin: return "MAT1-2"
    return "undetermined"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--profiles", type=Path, required=True)
    ap.add_argument("--reports", required=True)
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--lists", nargs="+", required=True, help="TSVs with genome, contig, start, group columns")
    ap.add_argument("--controls", type=int, default=200, help="random clear loci per control group")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    hmms = load_hmms(a.profiles)
    loci = read_reports(a.reports)
    targets: dict[tuple, set] = defaultdict(set)
    for lst in a.lists:
        for r in csv.DictReader(open(lst), delimiter="\t"):
            targets[(r["genome"], r["contig"], int(r["start"]))].add(
                r.get("groups") or ("label_disagrees" if r.get("label_disagrees") == "True" else "listed"))
            if r.get("label_disagrees") == "True":
                targets[(r["genome"], r["contig"], int(r["start"]))].add("label_disagrees")
    rng = random.Random(1)
    pools = defaultdict(list)
    for k, v in loci.items():
        if k in targets:
            continue
        e = v["entries"]
        a_, h_ = e.get("MAT1-1-1"), e.get("MAT1-2-1")
        if a_ and not h_ and modelled(a_): pools["control_alpha_only"].append(k)
        elif h_ and not a_ and modelled(h_): pools["control_hmg_only"].append(k)
        elif a_ and h_ and modelled(a_) and modelled(h_):
            gap = (a_.get("bitscore") or 0) - (h_.get("bitscore") or 0)
            if abs(gap) >= 100 and v["locus_class"] == "mat_locus": pools["control_both_clear"].append(k)
    for name, keys in pools.items():
        for k in rng.sample(sorted(keys), min(a.controls, len(keys))):
            targets[k].add(name)
    by_genome = defaultdict(list)
    for k in targets:
        if k in loci:
            by_genome[k[0]].append(k)
    print(f"{len(targets)} target loci in {len(by_genome)} genomes", file=sys.stderr)
    fields = ["genome", "contig", "start", "groups", "label", "confidence", "locus_class", "bitscore_winner",
              "alpha_on_alpha", "alpha_on_hmg", "hmg_on_hmg", "hmg_on_alpha", "profile_call", "profile_margin",
              "alpha_modelled", "hmg_modelled"]
    rows = []
    for gi, (genome, keys) in enumerate(sorted(by_genome.items())):
        wanted = {k[1] for k in keys}
        seqs = genome_contigs(a.genomes, genome, wanted)
        for k in keys:
            v = loci[k]; e = v["entries"]
            if k[1] not in seqs:
                continue
            sc = {g: best_scores(hmms, proteins(e[g], seqs[k[1]])) if g in e else None for g in GENES}
            a_on_a = sc["MAT1-1-1"]["MAT1-1-1"] if sc["MAT1-1-1"] else None
            a_on_h = sc["MAT1-1-1"]["MAT1-2-1"] if sc["MAT1-1-1"] else None
            h_on_h = sc["MAT1-2-1"]["MAT1-2-1"] if sc["MAT1-2-1"] else None
            h_on_a = sc["MAT1-2-1"]["MAT1-1-1"] if sc["MAT1-2-1"] else None
            if a_on_a is not None and h_on_h is not None:
                pc, margin = call(a_on_a, h_on_h), round(a_on_a - h_on_h, 1)
            else:  # one core gene: which profile its protein scores higher on
                own = sc["MAT1-1-1"] or sc["MAT1-2-1"]
                margin = round(own["MAT1-1-1"] - own["MAT1-2-1"], 1)
                pc = call(own["MAT1-1-1"], own["MAT1-2-1"])
            bw = None
            if "MAT1-1-1" in e and "MAT1-2-1" in e:
                bw = "MAT1-1" if (e["MAT1-1-1"].get("bitscore") or 0) > (e["MAT1-2-1"].get("bitscore") or 0) else "MAT1-2"
            rows.append(dict(genome=genome, contig=k[1], start=k[2], groups=",".join(sorted(targets[k])),
                             label=v["label"], confidence=v["confidence"], locus_class=v["locus_class"],
                             bitscore_winner=bw, alpha_on_alpha=a_on_a, alpha_on_hmg=a_on_h, hmg_on_hmg=h_on_h,
                             hmg_on_alpha=h_on_a, profile_call=pc, profile_margin=margin,
                             alpha_modelled=("MAT1-1-1" in e and modelled(e["MAT1-1-1"])),
                             hmg_modelled=("MAT1-2-1" in e and modelled(e["MAT1-2-1"]))))
        if gi % 50 == 0:
            print(f"  {gi + 1}/{len(by_genome)} genomes", file=sys.stderr)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: (round(v, 1) if isinstance(v, float) else v) for k, v in r.items()})
    print(f"wrote {len(rows)} loci -> {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
