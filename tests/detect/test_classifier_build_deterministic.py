"""Classifier builds are deterministic.

Curator's ruling, J. Stajich, 2026-09-29 (option c): make classifier builds
deterministic and keep --gate-only / --paralogs-only as fast paths.

Evidence (results/2026-09-29_gate_threshold/NOTE.md): rebuilding the PR #9
Mucoromycota classifier from the same training set and tool versions moved
held-out scores by up to 8.0 bits. The cause, measured 2026-09-29
(results/2026-09-29_deterministic_build/): MAFFT run with `--thread 4` returns
a different alignment on every run of the same input; `--thread 1` returns the
same one, and input order changes it too. The HMM files also carried the build
time (DATE) and the command line (COM).

Contract:
* two builds from the same inputs write byte-identical HMM files and identical
  leave-one-genus-out scores, in one process or in fresh processes;
* the input order of the training proteins does not change the HMM;
* MAFFT runs single-threaded;
* the manifest records the tool versions, the alignment options and a checksum
  of the training inputs.
"""
import hashlib
import random
import subprocess
import sys
from pathlib import Path

import pytest

from tests.detect.test_idiomorph_classifier import DB, SEXM, SEXP

REPO = Path(__file__).resolve().parents[2]


def _mutate(seq, n, seed):
    rng = random.Random(seed)
    s = list(seq)
    for i in rng.sample(range(len(s)), n):
        s[i] = rng.choice("ACDEFGHIKLMNPQRSTVWY")
    return "".join(s)


def _rows():
    """Six sexP and six sexM proteins in three genera (point-mutated curated ones)."""
    rows = []
    for gene, base in (("sexP", SEXP), ("sexM", SEXM)):
        for g, genus in enumerate(("Aaa", "Bbb", "Ccc")):
            for k in range(2):
                seq = _mutate(base, 25 + 10 * g + 3 * k, seed=hash((gene, g, k)) & 0xFFFF)
                rows.append(dict(id=f"{gene}_{genus}_{k}", gene=gene, genus=genus,
                                 source="fixture", sequence=seq))
    return rows


def _hmm_bytes(seqs, name, work):
    import io
    from MATPredict.detect.classifier_build import build_hmm
    hmm, _ = build_hmm(seqs, name, work)
    buf = io.BytesIO()
    hmm.write(buf)
    return buf.getvalue()


def test_build_hmm_is_byte_identical_across_calls(tmp_path):
    seqs = {r["id"]: r["sequence"] for r in _rows() if r["gene"] == "sexP"}
    a = _hmm_bytes(seqs, "sexP", tmp_path)
    b = _hmm_bytes(seqs, "sexP", tmp_path)
    assert a == b


def test_input_order_does_not_change_the_hmm(tmp_path):
    seqs = {r["id"]: r["sequence"] for r in _rows() if r["gene"] == "sexM"}
    rev = dict(reversed(list(seqs.items())))
    assert _hmm_bytes(seqs, "sexM", tmp_path) == _hmm_bytes(rev, "sexM", tmp_path)


def test_the_hmm_carries_no_build_time_or_command_line(tmp_path):
    seqs = {r["id"]: r["sequence"] for r in _rows() if r["gene"] == "sexP"}
    text = _hmm_bytes(seqs, "sexP", tmp_path).decode()
    assert "\nCOM " not in text
    date = next(l for l in text.splitlines() if l.startswith("DATE"))
    assert date.split(None, 1)[1].strip() == "Thu Jan  1 00:00:00 1970"


def test_mafft_runs_single_threaded(tmp_path, monkeypatch):
    from MATPredict.detect import classifier_build as cb
    seen = []
    real = subprocess.run

    def spy(cmd, *a, **kw):
        seen.append(list(cmd))
        return real(cmd, *a, **kw)

    monkeypatch.setattr(cb.subprocess, "run", spy)
    cb._mafft({"x": SEXP, "y": SEXM}, tmp_path, "t")
    (cmd,) = [c for c in seen if c and c[0].endswith("mafft")]
    assert cmd[cmd.index("--thread") + 1] == "1"
    for opt in cb.MAFFT_OPTIONS:
        assert opt in cmd


def _build(tmp_path, monkeypatch, rows):
    from MATPredict.detect import classifier_build as cb
    monkeypatch.setattr(cb, "training_set", lambda db_root, family, d: [dict(r) for r in rows])
    out = tmp_path
    out.mkdir(parents=True, exist_ok=True)
    return cb.build(DB, "Mucoromycota:MAT", out), out


def test_two_builds_are_byte_identical_and_score_the_same(tmp_path, monkeypatch):
    rows = _rows()
    m1, o1 = _build(tmp_path / "a", monkeypatch, rows)
    m2, o2 = _build(tmp_path / "b", monkeypatch, list(reversed(rows)))
    for gene in ("sexP", "sexM"):
        assert (o1 / f"{gene}.hmm").read_bytes() == (o2 / f"{gene}.hmm").read_bytes()
        assert m1["genes"][gene]["sha256"] == m2["genes"][gene]["sha256"]
    key = lambda r: r["id"]
    assert sorted(m1["leave_one_genus_out"]["margins"], key=key) == \
        sorted(m2["leave_one_genus_out"]["margins"], key=key)
    assert m1["training_sha256"] == m2["training_sha256"]


def test_the_manifest_records_versions_options_and_input_checksum(tmp_path, monkeypatch):
    from MATPredict.detect import classifier_build as cb
    m, _ = _build(tmp_path, monkeypatch, _rows())
    assert m["deterministic"] is True
    tools = m["tool_versions"]
    assert tools["pyhmmer"] and tools["hmmer"] and tools["mafft"]
    assert m["mafft_options"] == list(cb.MAFFT_OPTIONS)
    assert m["hmm_builder_seed"] == cb.BUILDER_SEED
    expect = hashlib.sha256("".join(
        f"{r['id']}\t{r['gene']}\t{r['genus']}\t{r['sequence']}\n"
        for r in sorted(_rows(), key=lambda r: (r["gene"], r["id"]))).encode()).hexdigest()
    assert m["training_sha256"] == expect


def test_builds_in_fresh_processes_are_byte_identical(tmp_path):
    """Two separate Python processes build the same HMM bytes."""
    seqs = {r["id"]: r["sequence"] for r in _rows() if r["gene"] == "sexP"}
    fa = tmp_path / "in.faa"
    fa.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
    code = (
        "import sys, io, pathlib, hashlib\n"
        "from MATPredict.detect.classifier_build import read_fasta, build_hmm\n"
        "p = pathlib.Path(sys.argv[1])\n"
        "hmm, _ = build_hmm(read_fasta(p), 'sexP', p.parent)\n"
        "b = io.BytesIO(); hmm.write(b)\n"
        "print(hashlib.sha256(b.getvalue()).hexdigest())\n")
    env_path = str(REPO / "src")
    import os
    env = dict(os.environ, PYTHONPATH=env_path + os.pathsep + os.environ.get("PYTHONPATH", ""))
    outs = [subprocess.run([sys.executable, "-c", code, str(fa)], capture_output=True, text=True,
                           env=env, check=True).stdout.strip() for _ in range(2)]
    assert outs[0] and outs[0] == outs[1]
