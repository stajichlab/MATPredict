"""Profile-HMM idiomorph classifier: decide a locus's idiomorph on its MODELLED
proteins.

Curator's ruling, 2026-09-26. After gene modelling, each modelled core protein
of an idiomorph-specific gene (e.g. Mucoromycota sexM / sexP) is scored against
one profile HMM per idiomorph; the idiomorph whose HMM scores higher wins, and
the margin between the two best scores is reported. A margin below the
family's `min_margin` leaves the call `undetermined` with both scores reported
(the 2026-09-21 ruling: report both rather than guess). The existing
model-pair / bitscore decision stays the fallback where the family has no
classifier.

Where NO core protein was modelled, the classifier scores the tblastn HSP
translations instead (curator's ruling 2026-09-27), keeping `min_margin`; the
verdict says so (`classifier_input: hsp_fragment`). Before, an identity
comparison of short hits decided those calls and mislabelled Umbelopsis sp.
M5902 (sexP +75.3 on its fragment, labelled Minus) and Mucor hiemalis
gzMucHiem1 (sexM -46.3, labelled Plus)
(results/2026-09-27_umbelopsis_diagnosis/NOTE.md).

Why a classifier and not more references: the 2026-09-26 experiment
(results/2026-09-26_sexMP_hmm/NOTE.md) scored modelled proteins, and sexM vs
sexP separated 38/38 leave-one-genus-out with a worst-case margin of 36.6 bits
(full length) against 5.1 for blastp. The 2026-09-21 a1/alpha1 HMM failed for
two reasons -- near-identical training sequences, and 6-frame ORF search --
and scoring MODELLED proteins removes the second; the build script
(`scripts/build_idiomorph_hmms.py`) reports training diversity for the first.

Files, never hand-edited, live in `db/<Phylum>/classifiers/<family>/`: one
`<gene>.hmm` per idiomorph-specific gene plus `manifest.yaml` (training IDs,
tool versions, checksums, leave-one-genus-out accuracy and margins). The
family's roster points at them:

    idiomorph_classifier: {type: hmm, dir: classifiers/MAT, min_margin: 20}

Non-MAT paralog classes (R4; curator's ruling 2026-09-29). A classifier may
also carry HMMs of known non-MAT HMG paralogs, listed in the manifest under
`paralog_classes` (`{name, hmm}` with `hmm` relative to the classifier
directory) and built by the build script from a curated source file in
`<dir>/paralogs/`. Each verdict records the paralog scores; when the best
paralog score beats the best MAT score by at least `min_margin` the verdict
names `paralog_class`, and the MAT-gene gate step withholds the call
(`mat_gene_gate.WITHHELD_PARALOG_CLASS`). Evidence: the "P1" HMG gene of the
Mucor hiemalis/indicus group, present in both mating types, scored 77-84 bits
as sexM but 145-159 against a P1 HMM; every other called core protein scored
>= 30 bits LOWER on P1 than on its own idiomorph
(results/2026-09-29_sexM_like_paralog/NOTE.md, p1_replay.tsv). Typing
(`idiomorph`) is unchanged by paralog classes.
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from pathlib import Path

import yaml

UNDETERMINED = "undetermined"


class ClassifierError(RuntimeError):
    """A roster names a classifier whose files are missing or inconsistent."""


@dataclass(frozen=True)
class IdiomorphClassifier:
    directory: Path
    #: idiomorph -> (gene name, pyhmmer HMM)
    models: dict
    min_margin: float
    manifest_sha256: str
    #: paralog class name -> pyhmmer HMM (R4); empty when the manifest lists none
    paralogs: dict = field(default_factory=dict)


@dataclass(frozen=True)
class ClassifierVerdict:
    #: idiomorph -> best HMM bit score over the proteins scored
    scores: dict[str, float]
    margin: float
    idiomorph: str
    min_margin: float
    proteins_scored: int
    manifest_sha256: str = ""
    method: str = "hmm"
    genes_scored: list[str] = field(default_factory=list)
    #: "model" (polished models), "hsp_fragment" (tblastn HSP translations,
    #: used only where no core protein was modelled), or "mixed" (a call
    #: spanning clusters of both kinds).
    classifier_input: str = "model"
    #: paralog class -> best HMM bit score (R4); empty without paralog classes
    paralog_scores: dict[str, float] = field(default_factory=dict)
    #: the paralog class that beats every MAT class by >= min_margin, else None
    paralog_class: str | None = None
    #: True when a scored model had exonerate frameshifts and was translated
    #: from its aligned blocks (possible assembly error; finding 026).
    frameshift_corrected: bool = False

    def as_report(self) -> dict:
        doc = {
            "method": self.method,
            "idiomorph": self.idiomorph,
            "scores": {k: round(v, 1) for k, v in sorted(self.scores.items())},
            "margin": round(self.margin, 1),
            "classifier_input": self.classifier_input,
            "min_margin": self.min_margin,
            "proteins_scored": self.proteins_scored,
            "genes_scored": sorted(self.genes_scored),
            "manifest_sha256": self.manifest_sha256,
        }
        if self.frameshift_corrected:
            doc["frameshift_corrected"] = True
        if self.paralog_scores:
            doc["paralog_scores"] = {k: round(v, 1) for k, v in sorted(self.paralog_scores.items())}
            doc["paralog_class"] = self.paralog_class
        return doc


_CACHE: dict[Path, IdiomorphClassifier] = {}


def load_classifier(spec: dict | None, family) -> IdiomorphClassifier | None:
    """The classifier a roster's `idiomorph_classifier` block names, or None.

    `spec["dir"]` is already resolved to an absolute path by the family
    loader. Every idiomorph-specific, idiomorph-informative core gene of the
    family must have a `<gene>.hmm`; a missing file is an error, not a silent
    fallback, because a roster that asks for a classifier and quietly runs
    without one would report calls on a basis nobody chose.
    """
    if not spec:
        return None
    if spec.get("type", "hmm") != "hmm":
        raise ClassifierError(f"unknown idiomorph_classifier type {spec.get('type')!r}")
    directory = Path(spec["dir"])
    if directory in _CACHE:
        return _CACHE[directory]
    import pyhmmer  # imported here so families without a classifier never need it

    manifest = directory / "manifest.yaml"
    if not manifest.exists():
        raise ClassifierError(f"no manifest.yaml in {directory}")
    models = {}
    for gene, idiomorph in classifier_genes(family).items():
        path = directory / f"{gene}.hmm"
        if not path.exists():
            raise ClassifierError(f"{path} missing for idiomorph {idiomorph}")
        with pyhmmer.plan7.HMMFile(str(path)) as handle:
            models[idiomorph] = (gene, handle.read())
    if len(models) < 2:
        raise ClassifierError(f"{directory}: need HMMs for at least two idiomorphs")
    paralogs = {}
    for entry in (yaml.safe_load(manifest.read_text()) or {}).get("paralog_classes") or []:
        path = directory / entry["hmm"]
        if not path.exists():
            raise ClassifierError(f"{path} missing for paralog class {entry['name']}")
        with pyhmmer.plan7.HMMFile(str(path)) as handle:
            paralogs[entry["name"]] = handle.read()
    clf = IdiomorphClassifier(
        directory=directory, models=models, min_margin=float(spec["min_margin"]),
        manifest_sha256=hashlib.sha256(manifest.read_bytes()).hexdigest(),
        paralogs=paralogs,
    )
    _CACHE[directory] = clf
    return clf


def classifier_genes(family) -> dict[str, str]:
    """gene name -> idiomorph, for the core genes a classifier models: core_MAT,
    restricted to exactly one idiomorph, and idiomorph-informative."""
    out = {}
    for g in family.genes:
        if not isinstance(g, dict) or g.get("role") != "core_MAT":
            continue
        if not g.get("idiomorph_informative", True):
            continue
        idioms = g.get("present_in_idiomorphs") or []
        if len(idioms) == 1:
            out[g["name"]] = idioms[0]
    return out


def _best_scores(hmms: dict, proteins: list[str]) -> dict[str, float]:
    """Best full-sequence bit score per named HMM over `proteins` (Z=1,
    permissive thresholds; a protein an HMM does not match scores 0.0)."""
    import pyhmmer

    alphabet = pyhmmer.easel.Alphabet.amino()
    seqs = [
        pyhmmer.easel.TextSequence(name=f"p{i}".encode(), sequence=p).digitize(alphabet)
        for i, p in enumerate(proteins) if p
    ]
    best = {name: 0.0 for name in hmms}
    if not seqs:
        return best
    block = pyhmmer.easel.DigitalSequenceBlock(alphabet, seqs)
    for name, hmm in hmms.items():
        pipeline = pyhmmer.plan7.Pipeline(alphabet, Z=1, E=1e9, domE=1e9, bias_filter=False,
                                          F1=1.0, F2=1.0, F3=1.0)
        for hit in pipeline.search_hmm(hmm, block):
            best[name] = max(best[name], float(hit.score))
    return best


def paralog_call(scores: dict[str, float], paralog_scores: dict[str, float],
                 min_margin: float) -> str | None:
    """The paralog class whose best score beats every MAT score by at least
    `min_margin`, else None."""
    if not paralog_scores:
        return None
    name, best = max(paralog_scores.items(), key=lambda kv: (kv[1], kv[0]))
    if best > 0 and best - max(scores.values(), default=0.0) >= min_margin:
        return name
    return None


def score_proteins(clf: IdiomorphClassifier, proteins: list[str]) -> dict[str, float]:
    """Best full-sequence HMM bit score per idiomorph over `proteins`.

    Scored with Z=1 and permissive thresholds so every protein gets a score
    (a protein the HMM does not match at all scores 0.0), exactly as the
    2026-09-26 experiment did with `hmmsearch -Z 1 -E 1000`.
    """
    return _best_scores({idiomorph: hmm for idiomorph, (_gene, hmm) in clf.models.items()},
                        proteins)


def classify(clf: IdiomorphClassifier, proteins: list[str], genes: list[str] | None = None,
             classifier_input: str = "model") -> ClassifierVerdict | None:
    """The idiomorph `proteins` belong to, or None when there is nothing to
    score. Below `clf.min_margin` the verdict is `undetermined`."""
    proteins = [p for p in proteins if p]
    if not proteins:
        return None
    scores = score_proteins(clf, proteins)
    ranked = sorted(scores.items(), key=lambda kv: (-kv[1], kv[0]))
    margin = ranked[0][1] - ranked[1][1]
    idiomorph = ranked[0][0] if margin >= clf.min_margin and ranked[0][1] > 0 else UNDETERMINED
    paralog_scores = _best_scores(clf.paralogs, proteins) if clf.paralogs else {}
    return ClassifierVerdict(
        scores=scores, margin=margin, idiomorph=idiomorph, min_margin=clf.min_margin,
        proteins_scored=len(proteins), manifest_sha256=clf.manifest_sha256,
        genes_scored=list(genes or []), classifier_input=classifier_input,
        paralog_scores=paralog_scores,
        paralog_class=paralog_call(scores, paralog_scores, clf.min_margin),
    )


def combine_verdicts(verdicts: list[ClassifierVerdict]) -> ClassifierVerdict | None:
    """One verdict for a result that spans several clusters (a fragmented
    call): the best score per idiomorph over all of them, re-decided against
    the same `min_margin`."""
    verdicts = [v for v in verdicts if v is not None]
    if not verdicts:
        return None
    if len(verdicts) == 1:
        return verdicts[0]
    scores: dict[str, float] = {}
    for v in verdicts:
        for k, x in v.scores.items():
            scores[k] = max(scores.get(k, 0.0), x)
    ranked = sorted(scores.items(), key=lambda kv: (-kv[1], kv[0]))
    margin = ranked[0][1] - ranked[1][1]
    first = verdicts[0]
    idiomorph = ranked[0][0] if margin >= first.min_margin and ranked[0][1] > 0 else UNDETERMINED
    paralog_scores: dict[str, float] = {}
    for v in verdicts:
        for k, x in v.paralog_scores.items():
            paralog_scores[k] = max(paralog_scores.get(k, 0.0), x)
    return ClassifierVerdict(
        paralog_scores=paralog_scores,
        paralog_class=paralog_call(scores, paralog_scores, first.min_margin),
        scores=scores, margin=margin, idiomorph=idiomorph, min_margin=first.min_margin,
        proteins_scored=sum(v.proteins_scored for v in verdicts),
        manifest_sha256=first.manifest_sha256,
        genes_scored=sorted({g for v in verdicts for g in v.genes_scored}),
        classifier_input=(inputs.pop() if len(inputs := {v.classifier_input for v in verdicts}) == 1
                          else "mixed"),
        frameshift_corrected=any(v.frameshift_corrected for v in verdicts),
    )


def read_manifest(directory: Path) -> dict:
    return yaml.safe_load((directory / "manifest.yaml").read_text())
