# tests/detect/test_polish_crash_isolation.py
"""One bad polish window must cost one gene, not an entire genome.

`exonerate --refine` segfaults on some windows. Measured on the real crash
from the 44-genus sweep (Gilbertella persicaria, scaffold_71:24582-37340, a
12,759 bp window against 19 reference proteins, exonerate 2.4.0):

* the crash is DETERMINISTIC, not transient -- three runs, three SIGSEGVs,
  byte-identical 184-line partial output each time;
* `--refine full` segfaults too, so it is the refinement step itself and not
  the `region` mode;
* dropping `--refine` succeeds, rc=0 with 21 gene records;
* miniprot on the same window succeeds, rc=0 with 12 mRNA records.

So a signal death is recoverable and must not be fatal. A NON-ZERO EXIT is a
different animal -- a missing binary or a bad path -- and stays fatal, or a
misconfiguration would silently mark every gene in every genome unpolished,
which is far worse than the 1-in-44 crash it would be papering over.

Returning None costs less than it sounds: `polish.classify` already turns a
missing exonerate model plus a miniprot model into `polished_single`, which
does NOT cap the confidence tier.
"""
from __future__ import annotations

import pytest

from MATPredict.detect.family_registry import Family, FamilyKey
from MATPredict.detect.search import SearchToolError, polish_with_exonerate

FAMILY = Family(
    FamilyKey("Basidiomycota", "aLocus"), "pattern", None, r"^a[0-9]+$",
    [{"name": "mfa1", "role": "core_MAT"}], [5270],
)

_GFF = (
    "c1\texonerate\tgene\t1\t400\t.\t+\t.\t"
    "gene_id 1 ; sequence rec1|gene0|mfa1 ; gene_orientation . ; identity 95.00 ; similarity 96.00\n"
    "c1\texonerate\texon\t1\t150\t.\t+\t.\tinsertions 0 ; deletions 0\n"
    "c1\texonerate\texon\t200\t400\t.\t+\t.\tinsertions 0 ; deletions 0\n"
)


class _Result:
    def __init__(self, returncode, stdout="", stderr=""):
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr


def _polish(tmp_path, runner):
    (tmp_path / "genome.fa").write_text(">c1\n" + "N" * 500 + "\n")
    return polish_with_exonerate(
        genome_fasta=tmp_path / "genome.fa",
        family=FAMILY, gene_name="mfa1",
        reference_fasta=tmp_path / "reference.faa",
        record_families={"rec1": FAMILY.key},
        window=("c1", 1, 500), runner=runner,
    )


def test_a_segfault_is_retried_without_refine(tmp_path):
    calls = []

    def runner(cmd, **kwargs):
        calls.append(list(cmd))
        if "--refine" in cmd:
            return _Result(-11)  # SIGSEGV, as measured on the real window
        return _Result(0, _GFF)

    model = _polish(tmp_path, runner)
    assert len(calls) == 2
    assert "--refine" in calls[0]
    assert "--refine" not in calls[1]
    assert model is not None
    assert model.exons


def test_the_retried_model_is_labelled_unrefined(tmp_path):
    # An unrefined model has less precise boundaries than a refined one, and
    # boundary precision is the entire point of the polish stage. It must not
    # be reported as though exonerate had refined it.
    def runner(cmd, **kwargs):
        return _Result(-11) if "--refine" in cmd else _Result(0, _GFF)

    assert _polish(tmp_path, runner).method == "exonerate_unrefined"


def test_a_successful_refine_is_not_retried_and_keeps_its_method(tmp_path):
    calls = []

    def runner(cmd, **kwargs):
        calls.append(list(cmd))
        return _Result(0, _GFF)

    assert _polish(tmp_path, calls and runner or runner).method == "exonerate_refine"
    assert len(calls) == 1


def test_a_segfault_on_both_attempts_returns_none_rather_than_raising(tmp_path):
    # classify() then yields polished_single from miniprot's model, with no
    # tier penalty. The locus survives; only exonerate's contribution is lost.
    def runner(cmd, **kwargs):
        return _Result(-11)

    assert _polish(tmp_path, runner) is None


def test_a_miniprot_segfault_returns_none_rather_than_killing_the_genome(tmp_path):
    # miniprot has no --refine to drop, so there is nothing to retry -- but a
    # signal death here would abort the whole genome exactly as the exonerate
    # one did, and classify() handles a missing miniprot model fine
    # (polished_single on exonerate alone, no tier penalty). No miniprot crash
    # has been observed; this is the same failure CLASS, closed symmetrically
    # rather than waiting to lose a genome to it.
    from MATPredict.detect.search import polish_with_miniprot

    (tmp_path / "genome.fa").write_text(">c1\n" + "N" * 500 + "\n")
    model = polish_with_miniprot(
        genome_fasta=tmp_path / "genome.fa",
        family=FAMILY, gene_name="mfa1",
        reference_fasta=tmp_path / "reference.faa",
        record_families={"rec1": FAMILY.key},
        window=("c1", 1, 500), runner=lambda cmd, **kw: _Result(-11),
    )
    assert model is None


def test_a_miniprot_nonzero_exit_still_raises(tmp_path):
    from MATPredict.detect.search import polish_with_miniprot

    (tmp_path / "genome.fa").write_text(">c1\n" + "N" * 500 + "\n")
    with pytest.raises(SearchToolError, match="exited 1"):
        polish_with_miniprot(
            genome_fasta=tmp_path / "genome.fa",
            family=FAMILY, gene_name="mfa1",
            reference_fasta=tmp_path / "reference.faa",
            record_families={"rec1": FAMILY.key},
            window=("c1", 1, 500), runner=lambda cmd, **kw: _Result(1),
        )


def test_a_nonzero_exit_still_raises_and_is_not_retried(tmp_path):
    # A bad path or a missing binary must stay loud. Degrading it would turn
    # a systematic misconfiguration into a silent per-gene downgrade.
    calls = []

    def runner(cmd, **kwargs):
        calls.append(list(cmd))
        return _Result(1, stderr="exonerate: could not open query file")

    with pytest.raises(SearchToolError, match="exited 1"):
        _polish(tmp_path, runner)
    assert len(calls) == 1
