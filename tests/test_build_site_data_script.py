"""scripts/build_site_data.py and build_locus_shards.py: release folders, both table generations, the note columns."""
import csv
import importlib.util
import io
import json
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

import pytest
import yaml

_ROOT = Path(__file__).resolve().parents[1]


def _load(name):
    sys.path.insert(0, str(_ROOT / "scripts"))
    spec = importlib.util.spec_from_file_location(name, _ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


bsd = _load("build_site_data")

ASCO_G = ["genome", "species", "class_", "order", "family", "status", "routing", "n_loci", "families_called", "wall_s"]
ASCO_L = ["genome", "class_", "order", "family_called", "contig", "start", "end", "idiomorph", "confidence", "locus_class",
          "detection_pass", "genes_found", "verification"]
NOTES = ["tier_note", "label_unsupported", "high_basis"]


def _write(path, cols, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n", restval="", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def _campaigns(tmp_path, notes=False, mucoro_tables=False):
    asco, basidio, mucoro = tmp_path / "asco", tmp_path / "basidio", tmp_path / "mucoro"
    gl = ASCO_L + (NOTES if notes else [])
    _write(asco / "genomes.tsv", ASCO_G, [
        dict(genome="A1", species="Aureobasidium pullulans", class_="Dothideomycetes", order="Dothideales", family="Saccotheciaceae",
             status="called", routing="lineage", n_loci=2),
        dict(genome="A2", species="Fusarium sp.", class_="Sordariomycetes", order="Hypocreales", family="Nectriaceae",
             status="uncalled", routing="lineage", n_loci=0)])
    _write(asco / "loci.tsv", gl, [
        dict(genome="A1", family_called="Ascomycota:MAT", contig="c1", start=1, end=9, idiomorph="MAT1-1", confidence="medium",
             locus_class="mat_locus", tier_note="vote_collapse", label_unsupported="bit_margin<20;genus:Aureobasidium"),
        dict(genome="A1", family_called="Ascomycota:MAT", contig="c2", start=1, end=9, idiomorph="MAT1-2", confidence="high",
             locus_class="mat_locus", high_basis="supported")])
    _write(basidio / "genomes.tsv", ASCO_G, [dict(genome="B1", species="Coprinopsis cinerea", class_="Agaricomycetes",
                                                  order="Agaricales", family="Psathyrellaceae", status="called", routing="lineage")])
    _write(basidio / "loci.tsv", gl, [dict(genome="B1", family_called="Basidiomycota:HD", contig="c1", start=1, end=9,
                                           idiomorph="undetermined", confidence="high", locus_class="mat_locus")])
    if mucoro_tables:
        _write(mucoro / "genomes.tsv", ["source", "genome", "species", "class_", "order", "family", "status", "routing"], [
            dict(source="BFD", genome="M1", species="Mucor lusitanicus", order="Mucorales", status="called"),
            dict(source="LCG", genome="L1", species="Absidia sp.", order="", status="uncalled")])
        _write(mucoro / "loci.tsv", ["source", "genome", "contig", "start", "end", "idiomorph", "confidence", "locus_class"] + NOTES, [
            dict(source="BFD", genome="M1", contig="c9", start=5, end=50, idiomorph="Plus", confidence="high", locus_class="mat_locus")])
    else:
        _write(mucoro / "calls.tsv", ["source", "genome", "name", "bfd_order", "status", "idiomorph", "confidence", "locus_class",
                                      "contig", "locus_start", "locus_end"], [
            dict(source="BFD", genome="M1", name="Mucor lusitanicus", bfd_order="Mucorales", status="called", idiomorph="Plus",
                 confidence="high", locus_class="mat_locus", contig="c9", locus_start=5, locus_end=50),
            dict(source="LCG", genome="L1", name="Absidia sp.", bfd_order="", status="uncalled")])
    return ["--asco", str(asco), "--basidio", str(basidio), "--mucoro", str(mucoro)]


def _genomes(path):
    with open(path, newline="") as fh:
        return {r["genome"]: r for r in csv.DictReader(fh, delimiter="\t")}


def test_v060_tables_get_blank_note_columns(tmp_path):
    out = tmp_path / "out"
    bsd.main(["--release", "v060", "--out", str(out)] + _campaigns(tmp_path))
    taxa = json.loads((out / "taxa.json").read_text())
    assert taxa["fields"][:14] == ["dataset", "phylum", "class", "order", "family", "genus", "species", "cohort", "n", "called",
                                   "undet", "two", "alleles", "types"] and taxa["fields"][14] == "unsupported"
    assert all(r[14] == 0 for r in taxa["rows"])
    g = _genomes(out / "genomes_bfd_asco.tsv")["A1"]
    assert g["loci"] == "c1:1-9|MAT1-1|medium|mat_locus;c2:1-9|MAT1-2|high|mat_locus"
    assert (g["label_unsupported"], g["tier_notes"], g["high_basis"], g["locus_notes"]) == ("", "", "", "")
    s = json.loads((out / "summary.json").read_text())
    assert s["release"]["name"] == "v060" and s["datasets"]["bfd_asco"]["code"] == "v0.6.0 (2026-10-03)"
    assert s["datasets"]["bfd_asco"]["label_unsupported"] == 0 and s["datasets"]["bfd_asco"]["tier_note"] == {}


def test_cleanrun_tables_pass_the_notes_on(tmp_path):
    out = tmp_path / "out"
    bsd.main(["--release", "cleanrun", "--out", str(out)] + _campaigns(tmp_path, notes=True, mucoro_tables=True))
    g = _genomes(out / "genomes_bfd_asco.tsv")["A1"]
    assert g["label_unsupported"] == "1" and g["tier_notes"] == "vote_collapse" and g["high_basis"] == "supported"
    assert g["locus_notes"] == "vote_collapse,bit_margin<20+genus:Aureobasidium;supported"
    s = json.loads((out / "summary.json").read_text())
    d = s["datasets"]["bfd_asco"]
    assert d["label_unsupported"] == 1 and d["tier_note"] == {"vote_collapse": 1} and d["high_basis"] == {"supported": 1}
    assert d["code"].startswith("cleanrun-2026-10b")
    taxa = json.loads((out / "taxa.json").read_text())
    assert sum(r[14] for r in taxa["rows"]) == 1


def test_mucoromycotina_new_tables_read_like_calls_tsv(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    bsd.main(["--release", "v060", "--out", str(a)] + _campaigns(tmp_path / "x"))
    bsd.main(["--release", "v060", "--out", str(b)] + _campaigns(tmp_path / "y", mucoro_tables=True))
    for ds in ("bfd_mucoro", "zygolife_lcg"):
        ga, gb = _genomes(a / f"genomes_{ds}.tsv"), _genomes(b / f"genomes_{ds}.tsv")
        notes = ("label_unsupported", "tier_notes", "high_basis", "locus_notes")
        assert {k: {c: v for c, v in r.items() if c not in notes} for k, r in ga.items()} == \
               {k: {c: v for c, v in r.items() if c not in notes} for k, r in gb.items()}
    m1 = _genomes(b / "genomes_bfd_mucoro.tsv")["M1"]
    assert (m1["label_unsupported"], m1["tier_notes"], m1["high_basis"], m1["locus_notes"]) == ("0", "", "", "-")


def test_missing_or_unset_folders_stop_the_build(tmp_path):
    args = _campaigns(tmp_path)
    with pytest.raises(SystemExit, match="does not exist"):
        bsd.main(["--release", "v060", "--out", str(tmp_path / "o"), "--asco", str(tmp_path / "nope")] + args[2:])
    with pytest.raises(SystemExit, match="no mucoro folder"):
        bsd.resolve_release("cleanrun", asco=args[1], basidio=args[3])


def test_both_releases_name_their_folders():
    for name, rel in bsd.RELEASES.items():
        assert {"code", "asco", "basidio", "mucoro"} <= set(rel)
    assert bsd.DEFAULT_RELEASE == "cleanrun"


@pytest.mark.skipif(shutil.which("zstd") is None, reason="needs the zstd command")
def test_locus_shards_append_the_notes_in_report_order(tmp_path):
    bls = _load("build_locus_shards")
    camp = tmp_path / "asco"
    _write(camp / "loci.tsv", ASCO_L + NOTES, [
        dict(genome="G1", contig="c1", start=1, end=9, tier_note="margin_fallback"),
        dict(genome="G1", contig="c2", start=1, end=9, label_unsupported="genus:Xylaria")])
    gene = dict(gene="MAT1-1-1", start=2, end=8, strand="+", role="core_MAT", status="polished_agree")
    rep = dict(detected=[dict(contig="c1", start=1, end=9, idiomorph="MAT1-1", confidence="medium", locus_class="mat_locus",
                              gene_evidence=[gene]),
                         dict(contig="c2", start=1, end=9, idiomorph="MAT1-1", confidence="high", locus_class="mat_locus")])
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w") as tf:
        data = yaml.safe_dump(rep).encode()
        ti = tarfile.TarInfo("./wave_0/small_000/runs/G1/detection_report.yaml")
        ti.size = len(data)
        tf.addfile(ti, io.BytesIO(data))
    (camp / "reports_all.tar.zst").write_bytes(subprocess.run(["zstd", "-q", "-c"], input=buf.getvalue(), capture_output=True,
                                                              check=True).stdout)
    shards, noted = bls.campaign_shards("bfd_asco", str(camp))
    assert noted
    assert shards["bfd_asco"]["G1"] == [
        ["c1", 1, 9, "MAT1-1", "medium", "mat_locus", [["MAT1-1-1", 2, 8, "+", "core_MAT", "polished_agree"]], "margin_fallback", "", ""],
        ["c2", 1, 9, "MAT1-1", "high", "mat_locus", [], "", "genus:Xylaria", ""]]
