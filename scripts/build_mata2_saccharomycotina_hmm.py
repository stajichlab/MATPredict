#!/usr/bin/env python
"""Build the Saccharomycotina MATa2 (HMG-box-like) profile HMM: SUPPORTING EVIDENCE ONLY.

    python scripts/build_mata2_saccharomycotina_hmm.py --out db/Ascomycota/classifiers/MATa2_saccharomycotina \
        [--pool pool.faa --exclude-ids ID,ID,...]

Curator decision 2026-10-18: build the seven-protein MATa2 profile into the database as supporting
evidence. It is NOT referenced by any order.yml, vote, gate or classification step and changes no call.
Training proteins are listed in `training.faa` in the output directory (annotated MATa2 proteins of
Komagataella, Kluyveromyces, Candida, Debaryomyces, Clavispora, Yarrowia, Geotrichum; the set may be
incomplete). The build is deterministic (MATPredict.detect.classifier_build: single-threaded MAFFT L-INS-i,
inputs sorted by id, fixed Builder seed and HMM header). With `--pool` (a protein FASTA used as the
background, e.g. the UniProt proteomes of 13 yeasts) the script also runs leave-one-lineage-out and writes
the held-out scores and margins into the manifest; without it those fields are copied from the existing
manifest if present.
"""
import argparse
import hashlib
import sys
import tempfile
from pathlib import Path

import yaml

from MATPredict.detect import classifier_build as cb

NAME = "MATa2_saccharomycotina"


def read_training(path: Path) -> list[dict]:
    rows = []
    for header, seq in cb.read_fasta(path).items():
        parts = header.split("|")
        kv = dict(p.split("=", 1) for p in parts[2:] if "=" in p)
        rows.append(dict(id=parts[0], gene=parts[1], genus=kv["genus"], lineage=kv["lineage"],
                         species=kv["species"], strain=kv["strain"], accession=kv["accession"],
                         record=kv["record"], source="training.faa", sequence=seq))
    return rows


def lineage_scores(rows, pool: dict[str, str], exclude: set[str], work: Path) -> dict:
    import pyhmmer

    aa = pyhmmer.easel.Alphabet.amino()
    pool = dict(pool)
    present = set(pool.values())
    for r in rows:  # training proteins that are not in the pool (e.g. curated records) are added to it
        if r["sequence"] not in present:
            pool["training|" + r["id"]] = r["sequence"]
    pool_seqs = [pyhmmer.easel.TextSequence(name=k.encode(), sequence=v).digitize(aa) for k, v in pool.items()]
    train_seqs = {r["sequence"] for r in rows}
    out = {}
    for lineage in sorted({r["lineage"] for r in rows}):
        train = {r["id"]: r["sequence"] for r in rows if r["lineage"] != lineage}
        held = [r for r in rows if r["lineage"] == lineage]
        hmm, _ = cb.build_hmm(train, f"LO_{lineage}", work)
        hits = list(pyhmmer.hmmsearch(hmm, pool_seqs, E=1000))[0]
        ranked = [(h.name if isinstance(h.name, str) else h.name.decode(), h.score, h.evalue) for h in hits]
        # background = pool proteins that are neither a training protein (exact sequence) nor a set-aside id
        bg = [(n, s) for n, s, e in ranked if pool[n] not in train_seqs and not (set(n.split("|")) & exclude)]
        best_bg = bg[0] if bg else (None, None)
        for r in held:
            mine = [(i, s, e) for i, (n, s, e) in enumerate(ranked) if pool[n] == r["sequence"]]
            if not mine:
                out[r["id"]] = dict(lineage=lineage, held_out_score=None)
                continue
            i, s, e = max(mine, key=lambda x: x[1])
            out[r["id"]] = dict(lineage=lineage, held_out_score=round(s, 1), evalue=f"{e:.0e}",
                                best_background_hmg_score=round(best_bg[1], 1) if best_bg[1] is not None else None,
                                best_background_id=best_bg[0],
                                margin=round(s - best_bg[1], 1) if best_bg[1] is not None else None)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--pool", type=Path)
    ap.add_argument("--exclude-ids", default="")
    a = ap.parse_args(argv)
    import pyhmmer

    rows = read_training(a.out / "training.faa")
    seqs = {r["id"]: r["sequence"] for r in rows}
    old = {}
    mpath = a.out / "manifest.yaml"
    if mpath.exists():
        old = yaml.safe_load(mpath.read_text()) or {}
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        hmm, aln = cb.build_hmm(seqs, NAME, work)
        with open(a.out / f"{NAME}.hmm", "wb") as fh:
            hmm.write(fh)
        uniq, spread = cb.identity_spread(aln)
        lolo = old.get("leave_one_lineage_out")
        if a.pool:
            pool = cb.read_fasta(a.pool)
            exclude = {x for x in a.exclude_ids.split(",") if x}
            lolo = dict(pool=f"{a.pool.name}: {len(pool)} proteins, sha256 "
                             f"{hashlib.sha256(a.pool.read_bytes()).hexdigest()}",
                        set_aside_background_ids=sorted(exclude),
                        per_protein=lineage_scores(rows, pool, exclude, work))
    manifest = dict(
        profile=NAME, kind="Saccharomycotina MATa2 (HMG-box-like) profile; SUPPORTING EVIDENCE ONLY",
        status_statement=(
            "This profile is NOT referenced by any order.yml, vote, gate or classification step and changes no call. "
            "The set of seven training proteins may be incomplete. Margins over the best background HMG-box "
            "protein are thin for single divergent lineages (see leave_one_lineage_out). MATa2 identity of the "
            "training proteins rests on annotation and the homology evidence in "
            "docs/notes/2026-10-18_komagataella-mata-analysis.md, not on a published domain annotation."),
        builder="scripts/build_mata2_saccharomycotina_hmm.py", deterministic=True,
        pyhmmer_version=pyhmmer.__version__, tool_versions=cb.tool_versions(),
        mafft_options=list(cb.MAFFT_OPTIONS), hmm_builder_seed=cb.BUILDER_SEED,
        training_sha256=cb.training_sha256(rows), hmm_sha256=cb._sha256(a.out / f"{NAME}.hmm"),
        n_sequences=len(rows), n_unique=uniq, pairwise_identity_min_median_max=spread,
        leave_one_lineage_out=lolo,
        leave_one_lineage_out_definition=(
            "For each lineage the profile is rebuilt without it and searched against the background pool; "
            "margin = held-out score minus the best-scoring pool protein that is not a training protein "
            "(exact sequence) and not in set_aside_background_ids (probable unlabelled MATa2 orthologues: "
            "Z. rouxii C5DUN6, L. thermotolerans C5DJY3, P. kudriavzevii A0A1Z8JKM8; locus-tag adjacency to the "
            "MATa1-type gene only, not checked on the genomes)."),
        training=[dict(id=r["id"], species=r["species"].replace("_", " "), strain=r["strain"].replace("_", " "),
                       lineage=r["lineage"], source_accession=r["accession"], record=r["record"])
                  for r in sorted(rows, key=lambda r: r["id"])])
    mpath.write_text(yaml.safe_dump(manifest, sort_keys=False))
    print(f"{len(rows)} proteins; hmm sha256 {manifest['hmm_sha256'][:12]}; "
          f"lineage scores: {'recomputed' if a.pool else 'kept' if lolo else 'absent'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
