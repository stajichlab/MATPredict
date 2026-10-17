#!/usr/bin/env python
"""Build profile HMMs for the two Ascomycota core mating-type genes: alpha box (MAT1-1-1) and HMG box (MAT1-2-1).

    python scripts/build_asco_core_hmms.py --db-root db --out db/Ascomycota/classifiers/MAT_core

These are gene-identity profiles, not allele profiles: the two genes are not homologous, so a protein scores
high on one profile and low on the other. They are built from the curated Ascomycota records only (plus
`training_extra.faa` in the output directory if present), deterministically, with the same options as the
Mucoromycota classifier (MATPredict.detect.classifier_build). Leave-one-genus-out is run and recorded.
The family's `order.yml` does not reference these profiles, so building them changes no call.
"""
import argparse
import sys
import tempfile
from pathlib import Path

import yaml

from MATPredict.detect import classifier_build as cb
from MATPredict.detect.family_registry import FamilyKey, load_all_families

GENES = {"MAT1-1-1": "MAT1-1", "MAT1-2-1": "MAT1-2"}
FAMILY = "Ascomycota:MAT"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db-root", type=Path, default=Path("db"))
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--notes", default="")
    a = ap.parse_args(argv)
    import pyhmmer

    family = next(f for f in load_all_families(a.db_root) if f.key == FamilyKey("Ascomycota", "MAT"))
    a.out.mkdir(parents=True, exist_ok=True)
    rows = [r for r in cb.training_set(a.db_root, family, a.out) if r["gene"] in GENES]
    per_gene, warnings = {}, []
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        for gene, idiomorph in GENES.items():
            seqs = {r["id"]: r["sequence"] for r in rows if r["gene"] == gene}
            if not seqs:
                raise SystemExit(f"no training sequences for {gene}")
            hmm, aln = cb.build_hmm(seqs, gene, work)
            with open(a.out / f"{gene}.hmm", "wb") as fh:
                hmm.write(fh)
            uniq, spread = cb.identity_spread(aln)
            per_gene[gene] = dict(
                idiomorph=idiomorph, n_sequences=len(seqs), n_unique=uniq,
                pairwise_identity_min_median_max=spread,
                genera=sorted({r["genus"] for r in rows if r["gene"] == gene}),
                sha256=cb._sha256(a.out / f"{gene}.hmm"))
            if uniq < cb.MIN_UNIQUE_SEQUENCES or (spread and spread[1] > cb.MAX_MEDIAN_IDENTITY):
                warnings.append(f"{gene}: low training diversity ({uniq} unique, median identity "
                                f"{spread[1] if spread else None}%)")
        loo = cb.leave_one_genus_out(rows, GENES, work)
    testable = [r for r in loo if not r["untestable_gene"]]
    correct = [r for r in testable if r["correct"]]
    manifest = dict(
        family=FAMILY, kind="gene-identity profiles (alpha box vs HMG box); not referenced by order.yml",
        builder="scripts/build_asco_core_hmms.py", deterministic=True,
        pyhmmer_version=pyhmmer.__version__, tool_versions=cb.tool_versions(),
        mafft_options=list(cb.MAFFT_OPTIONS), hmm_builder_seed=cb.BUILDER_SEED,
        training_sha256=cb.training_sha256(rows), notes=a.notes, warnings=warnings, genes=per_gene,
        leave_one_genus_out=dict(
            n_testable=len(testable), n_correct=len(correct),
            worst_correct_margin=min((r["margin"] for r in correct), default=None),
            failures=[r for r in testable if not r["correct"]],
            all_margins_sorted=sorted(r["margin"] for r in testable)[:10]),
        training=[dict(id=r["id"], gene=r["gene"], genus=r["genus"], source=r["source"]) for r in
                  sorted(rows, key=lambda r: (r["gene"], r["id"]))])
    (a.out / "manifest.yaml").write_text(yaml.safe_dump(manifest, sort_keys=False))
    print(f"alpha box: {per_gene['MAT1-1-1']['n_sequences']} sequences, HMG box: {per_gene['MAT1-2-1']['n_sequences']}; "
          f"leave-one-genus-out {len(correct)}/{len(testable)} correct; warnings: {warnings or 'none'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
