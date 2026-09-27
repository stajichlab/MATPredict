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
        warnings=warnings,
        notes=notes,
    )
    if min_margin is not None:
        manifest["min_margin_in_roster"] = min_margin
    with open(out_dir / "manifest.yaml", "w") as fh:
        yaml.safe_dump(manifest, fh, sort_keys=False, width=100)
    return manifest
