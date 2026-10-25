"""A polish alignment may not span an intron longer than the cap.

Xylaria bambusicola (GCA_033171785.1): with no limit, miniprot scored an 8.3 kb MAT1-2-1 model (two real HMG exons plus
two exons next to and inside SLA2, joined by a 5.5 kb and a 2.3 kb "intron" that hold APN2 and COX13) above the 331 bp
real fragment. The fixture is that locus's window (genome positions 69,400-78,400) and the reference protein that
produced the chimera."""
from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from MATPredict.detect.family_registry import Family, FamilyKey
from MATPredict.detect.search import MAX_INTRON_BP, _max_intron_bp, polish_with_exonerate, polish_with_miniprot

DATA = Path(__file__).parent / "data"
KEY = FamilyKey("Ascomycota", "MAT")
FAMILY = Family(KEY, "enum", ["MAT1-2"], None, [{"name": "MAT1-2-1", "role": "core_MAT"}], [4890])
OFFSET = 69_400 - 1          # window coordinate 1 is genome position 69,400
APN2_START, COX13_START = 71_464, 74_197


class _Result:
    def __init__(self, returncode=0, stdout=""):
        self.returncode, self.stdout, self.stderr = returncode, stdout, ""


def test_the_default_cap_is_5kb_and_the_env_var_overrides_it(monkeypatch):
    monkeypatch.delenv("MATPREDICT_MAX_INTRON_BP", raising=False)
    assert MAX_INTRON_BP == 5_000 == _max_intron_bp()
    monkeypatch.setenv("MATPREDICT_MAX_INTRON_BP", "2000")
    assert _max_intron_bp() == 2_000
    assert _max_intron_bp(1234) == 1_234        # an explicit argument wins over the environment


def _run(tool, tmp_path, **kwargs):
    cmds = []

    def runner(cmd, **kw):
        cmds.append(list(cmd))
        return _Result(0, "")

    (tmp_path / "g.fa").write_text(">c1\n" + "ACGT" * 150 + "\n")
    (tmp_path / "r.faa").write_text(">rec|gene0|MAT1-2-1\nMKLV\n")
    tool(tmp_path / "g.fa", FAMILY, "MAT1-2-1", tmp_path / "r.faa", {"rec": KEY}, ("c1", 1, 600), runner=runner, **kwargs)
    return cmds[0]


def test_miniprot_is_given_the_cap(tmp_path, monkeypatch):
    monkeypatch.delenv("MATPREDICT_MAX_INTRON_BP", raising=False)
    cmd = _run(polish_with_miniprot, tmp_path)
    assert cmd[cmd.index("-G") + 1] == "5000"
    assert cmd.index("-G") < cmd.index("--gff")      # options before the two positional files
    cmd = _run(polish_with_miniprot, tmp_path, max_intron=2000)
    assert cmd[cmd.index("-G") + 1] == "2000"


def test_exonerate_is_given_the_cap(tmp_path, monkeypatch):
    monkeypatch.delenv("MATPREDICT_MAX_INTRON_BP", raising=False)
    cmd = _run(polish_with_exonerate, tmp_path)
    assert cmd[cmd.index("--maxintron") + 1] == "5000"
    monkeypatch.setenv("MATPREDICT_MAX_INTRON_BP", "3000")
    cmd = _run(polish_with_exonerate, tmp_path)
    assert cmd[cmd.index("--maxintron") + 1] == "3000"


@pytest.mark.skipif(shutil.which("miniprot") is None, reason="miniprot is not installed")
def test_miniprot_without_a_cap_reproduces_the_chimera_and_with_one_it_does_not(tmp_path):
    seq = "".join((DATA / "xylaria_bambusicola_mat121_window.fa").read_text().splitlines()[1:])
    fa = tmp_path / "w.fa"
    fa.write_text(">c\n" + seq + "\n")
    window = ("c", 1, len(seq))
    kw = dict(genome_fasta=fa, family=FAMILY, gene_name="MAT1-2-1", reference_fasta=DATA / "xylaria_bambusicola_mat121_ref.faa",
              record_families={"5518_3639_MAT_combined": KEY}, window=window)
    chimera = polish_with_miniprot(**kw, max_intron=200_000)
    assert chimera is not None and chimera.end + OFFSET > COX13_START          # the model runs across APN2 and COX13
    for cap in (5_000, 2_000):
        model = polish_with_miniprot(**kw, max_intron=cap)
        assert model is not None
        assert model.end + OFFSET < APN2_START, f"cap {cap}: model still reaches {model.end + OFFSET}"
        assert all(b.start - a.end - 1 <= cap for a, b in zip(model.exons, model.exons[1:]))

def test_a_family_cap_sits_between_the_env_override_and_the_default(monkeypatch):
    fam = Family(KEY, "enum", ["MAT1-2"], None, [{"name": "MAT1-2-1", "role": "core_MAT"}], [4890], max_intron_bp=10_000)
    monkeypatch.delenv("MATPREDICT_MAX_INTRON_BP", raising=False)
    assert _max_intron_bp(None, fam) == 10_000
    assert _max_intron_bp(None, FAMILY) == MAX_INTRON_BP          # a family with no cap keeps the 5 kb default
    assert _max_intron_bp(None, None) == MAX_INTRON_BP
    monkeypatch.setenv("MATPREDICT_MAX_INTRON_BP", "3000")
    assert _max_intron_bp(None, fam) == 3_000                      # the run-wide override still wins
    assert _max_intron_bp(1234, fam) == 1_234                      # and an explicit argument wins over both


def test_both_aligners_get_the_family_cap(tmp_path, monkeypatch):
    global FAMILY
    monkeypatch.delenv("MATPREDICT_MAX_INTRON_BP", raising=False)
    saved = FAMILY
    FAMILY = Family(KEY, "enum", ["MAT1-2"], None, [{"name": "MAT1-2-1", "role": "core_MAT"}], [4890], max_intron_bp=10_000)
    try:
        cmd = _run(polish_with_miniprot, tmp_path)
        assert cmd[cmd.index("-G") + 1] == "10000"
        cmd = _run(polish_with_exonerate, tmp_path)
        assert cmd[cmd.index("--maxintron") + 1] == "10000"
    finally:
        FAMILY = saved


def test_only_the_basidiomycota_blocus_has_a_family_cap_in_the_curated_order_files():
    from MATPredict.detect.family_registry import load_all_families

    db = Path(__file__).resolve().parents[2] / "db"
    caps = {(f.key.phylum, f.key.locus_name): f.max_intron_bp for f in load_all_families(db) if f.max_intron_bp}
    assert caps == {("Basidiomycota", "bLocus"): 10_000}


def test_the_order_schema_accepts_a_positive_integer_cap_and_rejects_zero():
    from MATPredict.db.schema import validate_order
    import yaml

    doc = yaml.safe_load((Path(__file__).resolve().parents[2] / "db" / "Basidiomycota" / "order.yml").read_text())
    assert validate_order(doc) == []
    next(l for l in doc["loci"] if l["locus_name"] == "bLocus")["max_intron_bp"] = 0
    assert validate_order(doc) != []
