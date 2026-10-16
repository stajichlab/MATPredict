#!/usr/bin/env python
"""Re-check the Parastagonospora MAT1-1 / MAT1-2 calls of two runs with the alpha-box and HMG-box profiles.

    check_parastagonospora_profiles.py --profiles DIR --genome-list FILE --genomes DIR \
        --run v060 REPORTS_DIR --run records REPORTS_DIR --out checked.tsv

`--genome-list` has one genome id per line. For every Ascomycota:MAT locus of those genomes in each run it
scores the MAT1-1-1 and MAT1-2-1 entries on both profiles (score_core_profiles) and calls the locus by the
profile. The run's own label is reported next to it; it is never used to score.
"""
import argparse
import csv
import statistics as st
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import score_core_profiles as sp  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--profiles", type=Path, required=True)
    ap.add_argument("--genome-list", required=True)
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--run", nargs=2, action="append", metavar=("NAME", "REPORTS"), required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    wanted = {l.strip() for l in open(a.genome_list) if l.strip()}
    hmms = sp.load_hmms(a.profiles)
    rows = []
    for name, root in a.run:
        loci = {k: v for k, v in sp.read_reports(root).items() if k[0] in wanted}
        by_genome = defaultdict(list)
        for k in loci:
            by_genome[k[0]].append(k)
        print(f"{name}: {len(loci)} loci in {len(by_genome)} of {len(wanted)} genomes", file=sys.stderr)
        for genome, keys in sorted(by_genome.items()):
            seqs = sp.genome_contigs(a.genomes, genome, {k[1] for k in keys})
            for k in keys:
                e = loci[k]["entries"]
                if k[1] not in seqs or not e:
                    continue
                sc = {g: sp.best_scores(hmms, sp.proteins(e[g], seqs[k[1]])) if g in e else None for g in sp.GENES}
                aa = sc["MAT1-1-1"]["MAT1-1-1"] if sc["MAT1-1-1"] else None
                hh = sc["MAT1-2-1"]["MAT1-2-1"] if sc["MAT1-2-1"] else None
                if aa is not None and hh is not None:
                    pc, margin = sp.call(aa, hh), aa - hh
                else:
                    own = sc["MAT1-1-1"] or sc["MAT1-2-1"]
                    pc, margin = sp.call(own["MAT1-1-1"], own["MAT1-2-1"]), own["MAT1-1-1"] - own["MAT1-2-1"]
                rows.append(dict(run=name, genome=genome, contig=k[1], start=k[2], label=loci[k]["label"],
                                 genes="|".join(g for g in sp.GENES if g in e),
                                 alpha_modelled=("MAT1-1-1" in e and sp.modelled(e["MAT1-1-1"])),
                                 hmg_modelled=("MAT1-2-1" in e and sp.modelled(e["MAT1-2-1"])),
                                 alpha_on_alpha=aa, hmg_on_hmg=hh, profile_call=pc, profile_margin=round(margin, 1)))
    fields = list(rows[0].keys())
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: (round(v, 1) if isinstance(v, float) else v) for k, v in r.items()})
    print(f"\n{'run':9s} {'label':14s} {'n':>4s}  profile call (MAT1-1 / MAT1-2 / undetermined)   median alpha-on-alpha  median HMG-on-HMG")
    for name, _ in a.run:
        for label in sorted({r["label"] for r in rows if r["run"] == name}):
            sel = [r for r in rows if r["run"] == name and r["label"] == label]
            c = Counter(r["profile_call"] for r in sel)
            med = lambda f: (round(st.median([r[f] for r in sel if r[f] is not None]), 1)
                             if any(r[f] is not None for r in sel) else None)
            print(f"{name:9s} {label:14s} {len(sel):4d}  {c['MAT1-1']:4d} / {c['MAT1-2']:4d} / {c['undetermined']:4d}"
                  f"{'':20s}{med('alpha_on_alpha')!s:>8s}{'':14s}{med('hmg_on_hmg')!s:>8s}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
