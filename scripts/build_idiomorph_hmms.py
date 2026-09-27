#!/usr/bin/env python
"""Build a family's profile-HMM idiomorph classifier from the curated database.

The only supported way to make the files in db/<Phylum>/classifiers/<family>/
(curator's ruling 2026-09-26). Never edit a .hmm by hand; rerun this.

    python scripts/build_idiomorph_hmms.py --db-root db --family Mucoromycota:MAT \
        --out db/Mucoromycota/classifiers/MAT [--extra-fasta FILE]

--extra-fasta copies a documented extra training file into the classifier
directory as training_extra.faa (headers `>{id}|{gene}|genus={genus}`); later
rebuilds reuse it. Writes <gene>.hmm per idiomorph-specific core gene and
manifest.yaml (training IDs, versions, checksums, training diversity,
leave-one-genus-out accuracy and margins, recommended min_margin).
"""
import argparse
import shutil
import sys
from pathlib import Path

from MATPredict.detect.classifier_build import build


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db-root", type=Path, default=Path("db"))
    ap.add_argument("--family", required=True, help="Phylum:Locus, e.g. Mucoromycota:MAT")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--extra-fasta", type=Path)
    ap.add_argument("--notes", default="", help="free text recorded in manifest.yaml")
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    if args.extra_fasta:
        shutil.copyfile(args.extra_fasta, args.out / "training_extra.faa")
    m = build(args.db_root, args.family, args.out, notes=args.notes)
    loo = m["leave_one_genus_out"]
    print(f"{args.family}: LOO {loo['n_correct']}/{loo['n_tested']} correct "
          f"({loo['n_untestable']} untestable), worst correct margin {loo['worst_correct_margin']} bits; "
          f"recommended min_margin {m['recommended_min_margin']}")
    for g, d in m["genes"].items():
        print(f"  {g}: n={d['n_sequences']} unique={d['n_unique']} genera={d['n_genera']} "
              f"identity(min/median/max)={d['pairwise_identity_min_median_max']}")
    for w in m["warnings"]:
        print(f"  WARNING {w}", file=sys.stderr)


if __name__ == "__main__":
    main()
