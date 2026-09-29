"""Build a family's profile-HMM idiomorph classifier, and measure it.

The ONLY way classifier files are made (curator's ruling 2026-09-26): nobody
edits a `.hmm` by hand, so the files in `db/<Phylum>/classifiers/<family>/`
are always reproducible from the curated database plus the one recorded extra
training file, and every build re-measures itself.

Training set, per idiomorph-specific core gene (`classifier.classifier_genes`):

* every accepted curated record's protein for that gene
  (`reference_fasta.build_reference_fasta`, restricted to the family), and
* `training_extra.faa` in the classifier directory, if present -- sequences
  admitted by a documented rule, each header `>{id}|{gene}|genus={genus}`.
  For Mucoromycota:MAT these are members of the UFBoot-99 sexP clade of the
  2026-09-26 HMG-box tree; the Zygo 23 loci are deliberately NOT in it, so
  the Zygo regression stays an independent test (curator's option (a)).

Each gene's sequences are aligned with MAFFT and built into an HMM with
pyhmmer. The build also runs leave-one-genus-out: for every genus, both HMMs
are rebuilt without it and its sequences are scored, so the manifest records
how often the right idiomorph wins and by how much.

Non-MAT paralog classes (R4; curator's ruling 2026-09-29): every
`paralogs/<name>.faa` in the classifier directory, with its provenance in
`paralogs/<name>.yaml` (`source`, `reason`, `evidence` -- refused without),
is built into `paralogs/<name>.hmm` and listed in the manifest under
`paralog_classes`. A one-sequence class is built from the sequence itself
(`Builder.build`, as the 2026-09-29 replay did), so it has no
leave-one-genus-out; instead the build checks that no MAT training protein
would be classed as the paralog and that the source sequences themselves are.

Training diversity is reported and checked: the 2026-09-21 a1/alpha1 HMM
failed partly because 778 peptides collapsed to 56 near-identical sequences
(78-99% identity). A gene whose training set has fewer than
`MIN_UNIQUE_SEQUENCES` unique sequences, or whose median pairwise identity is
above `MAX_MEDIAN_IDENTITY`, is flagged in the manifest and on stderr.
"""
from __future__ import annotations

import datetime
import hashlib
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

import yaml

from MATPredict.detect.classifier import classifier_genes
from MATPredict.detect.family_registry import FamilyKey, load_all_families
from MATPredict.detect.reference_fasta import build_reference_fasta

MIN_UNIQUE_SEQUENCES = 5
MAX_MEDIAN_IDENTITY = 80.0
MIN_PROTEIN_LENGTH = 50


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    cur = None
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            cur = line[1:].split()[0]
            seqs[cur] = []
        elif cur is not None:
            seqs[cur].append(line.strip())
    return {k: "".join(v).replace("*", "") for k, v in seqs.items()}


def write_fasta(path: Path, seqs: dict[str, str]) -> None:
    with open(path, "w") as fh:
        for k, v in seqs.items():
            fh.write(f">{k}\n{v}\n")


def _record_genera(db_root: Path) -> dict[str, str]:
    out = {}
    for meta in db_root.glob("*/*/*/metadata.yaml"):
        try:
            sp = yaml.safe_load(meta.read_text())["organism"]["species"]
            out[meta.parent.name] = sp.split()[0]
        except Exception:
            continue
    return out


def training_set(db_root: Path, family, classifier_dir: Path) -> list[dict]:
    """One dict per training protein: id, gene, genus, source, sequence."""
    genes = classifier_genes(family)
    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        ref = build_reference_fasta(db_root, Path(tmp) / "ref.faa", family_keys={family.key})
        genera = _record_genera(db_root)
        for header, seq in read_fasta(ref).items():
            record, _idx, gene = header.split("|")
            if gene in genes and len(seq) >= MIN_PROTEIN_LENGTH:
                rows.append(dict(id=f"REF|{header}", gene=gene, genus=genera.get(record, record),
                                 source="curated_record", sequence=seq))
    extra = classifier_dir / "training_extra.faa"
    if extra.exists():
        for header, seq in read_fasta(extra).items():
            parts = header.split("|")
            gene = parts[1]
            genus = next((p.split("=", 1)[1] for p in parts if p.startswith("genus=")), "?")
            if gene in genes and len(seq) >= MIN_PROTEIN_LENGTH:
                rows.append(dict(id=header, gene=gene, genus=genus, source="training_extra",
                                 sequence=seq))
    return rows


def _mafft(seqs: dict[str, str], workdir: Path, tag: str) -> Path:
    fa, aln = workdir / f"{tag}.faa", workdir / f"{tag}.afa"
    write_fasta(fa, seqs)
    mafft = shutil.which("mafft")
    if mafft is None:
        raise RuntimeError("mafft not found on PATH (it is in the pixi environment)")
    with open(aln, "w") as out:
        subprocess.run([mafft, "--auto", "--quiet", "--thread", "4", str(fa)],
                       stdout=out, check=True)
    return aln


def build_hmm(seqs: dict[str, str], name: str, workdir: Path):
    import pyhmmer

    alphabet = pyhmmer.easel.Alphabet.amino()
    if len(seqs) == 1:  # MAFFT needs two; a one-sequence "alignment" is the sequence
        (k, v), = seqs.items()
        seqs = {k: v, f"{k}_dup": v}
    aln = _mafft(seqs, workdir, name)
    with pyhmmer.easel.MSAFile(str(aln), digital=True, alphabet=alphabet) as fh:
        msa = fh.read()
    msa.name = name.encode()
    builder = pyhmmer.plan7.Builder(alphabet)
    hmm, _, _ = builder.build_msa(msa, pyhmmer.plan7.Background(alphabet))
    return hmm, aln


def identity_spread(aln: Path) -> tuple[int, list[float] | None]:
    a = read_fasta(aln)
    ks = list(a)
    vals = []
    for i in range(len(ks)):
        for j in range(i + 1, len(ks)):
            cols = [(x, y) for x, y in zip(a[ks[i]], a[ks[j]]) if x != "-" and y != "-"]
            if cols:
                vals.append(100 * sum(x == y for x, y in cols) / len(cols))
    uniq = len({v.replace("-", "") for v in a.values()})
    if not vals:
        return uniq, None
    vals.sort()
    return uniq, [round(vals[0], 1), round(vals[len(vals) // 2], 1), round(vals[-1], 1)]


def score(hmm, seqs: dict[str, str]) -> dict[str, float]:
    import pyhmmer

    alphabet = pyhmmer.easel.Alphabet.amino()
    block = pyhmmer.easel.DigitalSequenceBlock(alphabet, [
        pyhmmer.easel.TextSequence(name=k.encode(), sequence=v).digitize(alphabet)
        for k, v in seqs.items()
    ])
    pipeline = pyhmmer.plan7.Pipeline(alphabet, Z=1, E=1e9, domE=1e9, bias_filter=False,
                                      F1=1.0, F2=1.0, F3=1.0)
    out = {k: 0.0 for k in seqs}
    for hit in pipeline.search_hmm(hmm, block):
        name = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
        out[name] = max(out[name], float(hit.score))
    return out


def leave_one_genus_out(rows: list[dict], genes: dict[str, str], workdir: Path) -> list[dict]:
    """For each genus: rebuild every gene's HMM without it, score its proteins.
    A genus that is a gene's ONLY training source leaves that gene with no
    model; its proteins are then scored against the remaining HMMs only and
    marked `untestable_gene`."""
    out = []
    by_genus = defaultdict(list)
    for r in rows:
        by_genus[r["genus"]].append(r)
    for genus, held in sorted(by_genus.items()):
        train = [r for r in rows if r["genus"] != genus]
        hmms = {}
        for gene in genes:
            seqs = {r["id"]: r["sequence"] for r in train if r["gene"] == gene}
            if seqs:
                hmms[gene], _ = build_hmm(seqs, f"loo_{gene}", workdir)
        test = {r["id"]: r["sequence"] for r in held}
        scores = {gene: score(h, test) for gene, h in hmms.items()}
        for r in held:
            s = {genes[g]: scores[g][r["id"]] for g in hmms}
            own = genes[r["gene"]]
            other = max((v for k, v in s.items() if k != own), default=0.0)
            out.append(dict(id=r["id"], genus=genus, gene=r["gene"], truth=own,
                            own_score=round(s.get(own, 0.0), 1), other_score=round(other, 1),
                            margin=round(s.get(own, 0.0) - other, 1),
                            correct=own in s and s[own] > other,
                            untestable_gene=r["gene"] not in hmms))
    return out


#: Percentile of the HMG-paralog negative set's best scores that sets the
#: MAT-gene gate threshold. In the 2026-09-28 validation build the paralog
#: 95th percentile was 99.9 bits -- the 100 bits the curator adopted
#: (results/2026-09-28_validation_f3_f4/f3_absolute.txt: >=100 kept 96/108
#: true full proteins and 9/189 paralogs). Anchoring on the negatives makes the
#: threshold move WITH a rebuild's score calibration.
GATE_PARALOG_PERCENTILE = 95
NEGATIVES_FILE = "paralog_negatives.faa"


def gate_threshold(best_scores: list[float], percentile: int = GATE_PARALOG_PERCENTILE):
    """The nearest-rank `percentile` of `best_scores`, rounded to 0.1, or None
    when there are none."""
    if not best_scores:
        return None
    ordered = sorted(best_scores)
    rank = max(1, -(-percentile * len(ordered) // 100))  # ceil(p * n / 100)
    return round(ordered[rank - 1], 1)


def compute_gate(out_dir: Path, genes: dict[str, str], loo: list[dict]) -> dict | None:
    """The `mat_gene_gate` manifest block for the HMMs in `out_dir`, or None
    (with a warning) when the classifier has no negative set."""
    import pyhmmer

    negatives_path = out_dir / NEGATIVES_FILE
    if not negatives_path.exists():
        print(f"WARNING no {NEGATIVES_FILE} in {out_dir}: manifest gets no mat_gene_gate "
              "threshold; the gate falls back to the roster/default", file=sys.stderr)
        return None
    negatives = read_fasta(negatives_path)
    per = []
    for gene in genes:
        with pyhmmer.plan7.HMMFile(str(out_dir / f"{gene}.hmm")) as fh:
            per.append(score(fh.read(), negatives))
    best = [max(p[k] for p in per) for k in negatives]
    threshold = gate_threshold(best)
    own = [r["own_score"] for r in loo if r["correct"] and not r["untestable_gene"]]
    return dict(
        min_score=threshold,
        rule=(f"nearest-rank {GATE_PARALOG_PERCENTILE}th percentile of the best classifier "
              f"scores of the HMG-paralog negative set ({NEGATIVES_FILE}) against this "
              "build's HMMs; curator ruling 2026-09-29"),
        n_negatives=len(best),
        negatives_sha256=_sha256(negatives_path),
        negatives_at_or_above=sum(b >= threshold for b in best),
        loo_correct_at_or_above=sum(o >= threshold for o in own),
        loo_correct_tested=len(own),
    )


def update_gate(db_root: Path, family_key: str, out_dir: Path) -> dict | None:
    """Recompute `mat_gene_gate` for the EXISTING HMMs in `out_dir` and write it
    into their manifest, without rebuilding the HMMs. Rebuilds are not
    bit-reproducible (the same training set and tools moved held-out scores by
    up to 8 bits on 2026-09-29), so adding the threshold must not rebuild."""
    phylum, locus = family_key.split(":")
    family = next(f for f in load_all_families(db_root)
                  if f.key == FamilyKey(phylum=phylum, locus_name=locus))
    manifest = yaml.safe_load((out_dir / "manifest.yaml").read_text())
    gate = compute_gate(out_dir, classifier_genes(family),
                        manifest["leave_one_genus_out"]["margins"])
    manifest["mat_gene_gate"] = gate
    with open(out_dir / "manifest.yaml", "w") as fh:
        yaml.safe_dump(manifest, fh, sort_keys=False, width=100)
    return gate


PARALOG_DIR = "paralogs"
PARALOG_PROVENANCE_FIELDS = ("source", "reason", "evidence")


def _paralog_hmm(seqs: dict[str, str], name: str, workdir: Path):
    """An HMM for a paralog class: from the single sequence itself when there
    is one (the 2026-09-29 replay's construction), else from a MAFFT MSA."""
    import pyhmmer

    if len(seqs) > 1:
        hmm, _ = build_hmm(seqs, name, workdir)
        return hmm
    alphabet = pyhmmer.easel.Alphabet.amino()
    (k, v), = seqs.items()
    seq = pyhmmer.easel.TextSequence(name=k.encode(), sequence=v).digitize(alphabet)
    hmm, _, _ = pyhmmer.plan7.Builder(alphabet).build(seq, pyhmmer.plan7.Background(alphabet))
    hmm.name = name.encode()
    return hmm


def build_paralog_classes(out_dir: Path, genes: dict[str, str], min_margin: float,
                          check_proteins: dict) -> list[dict]:
    """Build `paralogs/<name>.hmm` for every curated `paralogs/<name>.faa` in
    `out_dir` and return the manifest's `paralog_classes` block.

    `genes` maps each MAT gene to its idiomorph (its HMM is `<gene>.hmm` in
    `out_dir`); `check_proteins` maps an id to `(gene, sequence)` for the MAT
    training proteins, each of which must NOT be classed as a paralog."""
    import pyhmmer

    from MATPredict.detect.classifier import paralog_call

    pdir = out_dir / PARALOG_DIR
    if not pdir.is_dir():
        return []
    mat = {}
    for gene, idiomorph in genes.items():
        with pyhmmer.plan7.HMMFile(str(out_dir / f"{gene}.hmm")) as fh:
            mat[idiomorph] = fh.read()
    block = []
    with tempfile.TemporaryDirectory() as tmp:
        for faa in sorted(pdir.glob("*.faa")):
            name = faa.stem
            prov_path = pdir / f"{name}.yaml"
            provenance = yaml.safe_load(prov_path.read_text()) if prov_path.exists() else None
            if not provenance or any(not provenance.get(k) for k in PARALOG_PROVENANCE_FIELDS):
                raise RuntimeError(f"{faa}: paralog class needs {prov_path.name} with "
                                   f"{', '.join(PARALOG_PROVENANCE_FIELDS)}")
            seqs = read_fasta(faa)
            hmm = _paralog_hmm(seqs, name, Path(tmp))
            hmm_path = pdir / f"{name}.hmm"
            with open(hmm_path, "wb") as fh:
                hmm.write(fh)

            def classed(proteins: dict[str, str]) -> dict[str, bool]:
                per = {i: score(h, proteins) for i, h in mat.items()} if proteins else {}
                par = score(hmm, proteins) if proteins else {}
                return {k: paralog_call({i: per[i][k] for i in per}, {name: par[k]},
                                        min_margin) is not None for k in proteins}

            check = classed({k: v for k, (_g, v) in check_proteins.items()})
            own = classed(seqs)
            block.append(dict(
                name=name, hmm=f"{PARALOG_DIR}/{name}.hmm",
                n_sequences=len(seqs), sequence_ids=sorted(seqs),
                source_sha256=_sha256(faa), sha256=_sha256(hmm_path),
                provenance=provenance,
                rule=(f"a protein whose best {name} score beats its best MAT score by >= "
                      f"min_margin ({min_margin:g}) is classed {name} and its call withheld"),
                training_mat_proteins_checked=len(check),
                training_mat_proteins_classed_paralog=sum(check.values()),
                self_check_classed_paralog=sum(own.values()),
                leave_one_out=("not possible: one sequence" if len(seqs) == 1 else
                               "not run"),
            ))
    return block


def update_paralogs(db_root: Path, family_key: str, out_dir: Path) -> list[dict]:
    """Build the paralog classes for the EXISTING MAT HMMs in `out_dir` and
    write them into the manifest, without rebuilding the MAT HMMs (rebuilds
    are not bit-reproducible; see `update_gate`)."""
    phylum, locus = family_key.split(":")
    family = next(f for f in load_all_families(db_root)
                  if f.key == FamilyKey(phylum=phylum, locus_name=locus))
    genes = classifier_genes(family)
    rows = training_set(db_root, family, out_dir)
    min_margin = float((family.idiomorph_classifier or {}).get("min_margin", 25.0))
    block = build_paralog_classes(out_dir, genes, min_margin,
                                  {r["id"]: (r["gene"], r["sequence"]) for r in rows})
    manifest = yaml.safe_load((out_dir / "manifest.yaml").read_text())
    manifest["paralog_classes"] = block
    with open(out_dir / "manifest.yaml", "w") as fh:
        yaml.safe_dump(manifest, fh, sort_keys=False, width=100)
    return block


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build(db_root: Path, family_key: str, out_dir: Path, min_margin: float | None = None,
          notes: str = "") -> dict:
    """Build `<gene>.hmm` for every idiomorph-specific gene of `family_key`
    into `out_dir`, run leave-one-genus-out, and write `manifest.yaml`.
    Returns the manifest."""
    import pyhmmer

    phylum, locus = family_key.split(":")
    family = next(f for f in load_all_families(db_root)
                  if f.key == FamilyKey(phylum=phylum, locus_name=locus))
    genes = classifier_genes(family)
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = training_set(db_root, family, out_dir)
    warnings = []
    per_gene = {}
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        for gene, idiomorph in genes.items():
            seqs = {r["id"]: r["sequence"] for r in rows if r["gene"] == gene}
            if not seqs:
                raise RuntimeError(f"{family_key}: no training sequences for {gene}")
            hmm, aln = build_hmm(seqs, gene, work)
            with open(out_dir / f"{gene}.hmm", "wb") as fh:
                hmm.write(fh)
            uniq, spread = identity_spread(aln)
            genera = sorted({r["genus"] for r in rows if r["gene"] == gene})
            per_gene[gene] = dict(
                idiomorph=idiomorph, n_sequences=len(seqs), n_unique=uniq,
                pairwise_identity_min_median_max=spread, n_genera=len(genera),
                genera=genera, sha256=_sha256(out_dir / f"{gene}.hmm"),
                sources={s: sum(1 for r in rows if r["gene"] == gene and r["source"] == s)
                         for s in sorted({r["source"] for r in rows})},
            )
            if uniq < MIN_UNIQUE_SEQUENCES or (spread and spread[1] > MAX_MEDIAN_IDENTITY):
                msg = (f"{gene}: low training diversity ({uniq} unique, median identity "
                       f"{spread[1] if spread else None}%) -- see memory "
                       f"matpredict_hmm_typing_rejected")
                warnings.append(msg)
                print(f"WARNING {msg}", file=sys.stderr)
        loo = leave_one_genus_out(rows, genes, work)
        gate = compute_gate(out_dir, genes, loo)
        paralog_classes = build_paralog_classes(
            out_dir, genes, float((family.idiomorph_classifier or {}).get("min_margin", 25.0)),
            {r["id"]: (r["gene"], r["sequence"]) for r in rows})
    testable = [r for r in loo if not r["untestable_gene"]]
    correct = [r for r in testable if r["correct"]]
    worst = min((r["margin"] for r in correct), default=None)
    manifest = dict(
        family=family_key,
        built=datetime.date.today().isoformat(),
        builder="scripts/build_idiomorph_hmms.py (MATPredict.detect.classifier_build)",
        pyhmmer_version=pyhmmer.__version__,
        mafft=subprocess.run(["mafft", "--version"], capture_output=True, text=True).stderr.strip(),
        genes=per_gene,
        training=[{k: r[k] for k in ("id", "gene", "genus", "source")} for r in rows],
        leave_one_genus_out=dict(
            n_tested=len(testable), n_correct=len(correct),
            n_untestable=len(loo) - len(testable),
            worst_correct_margin=worst,
            margins=[{k: r[k] for k in ("id", "genus", "truth", "own_score", "other_score",
                                          "margin", "correct", "untestable_gene")} for r in loo],
        ),
        recommended_min_margin=(None if worst is None else round(worst / 2, 1)),
        min_margin_rule=("half the worst correct leave-one-genus-out margin, so every held-out "
                         "training protein is still called and a margin below half the weakest "
                         "separation seen on held-out genera is reported undetermined"),
        mat_gene_gate=gate,
        paralog_classes=paralog_classes,
        warnings=warnings,
        notes=notes,
    )
    if min_margin is not None:
        manifest["min_margin_in_roster"] = min_margin
    with open(out_dir / "manifest.yaml", "w") as fh:
        yaml.safe_dump(manifest, fh, sort_keys=False, width=100)
    return manifest
