"""scripts/build_site_data.py and build_locus_shards.py: release folders and their codes, both table generations,
the reporting columns, the held-out sets, the shard notes."""
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
NOTES = ["baseline_confidence", "tier_note", "label_unsupported", "label_unsupported_reason", "high_basis", "high_basis_outside_scope",
         "defective_genes"]
MANIFEST = "clean full-BFD run\ntag: cleanrun-2026-10b  commit: 56e11cd74c8e86bd43045efa543ace8b96abc745\n"


def _write(path, cols, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n", restval="", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def _campaigns(tmp_path, notes=False, mucoro_tables=False, manifest=None, defective_col=True):
    asco, basidio, mucoro = tmp_path / "asco", tmp_path / "basidio", tmp_path / "mucoro"
    note_cols = [c for c in NOTES if defective_col or c != "defective_genes"]    # a clean-run table from before #151 lacks the last
    gl = ASCO_L + (note_cols if notes else [])
    _write(asco / "genomes.tsv", ASCO_G, [
        dict(genome="A1", species="Aureobasidium pullulans", class_="Dothideomycetes", order="Dothideales", family="Saccotheciaceae",
             status="called", routing="lineage", n_loci=2),
        dict(genome="A2", species="Fusarium sp.", class_="Sordariomycetes", order="Hypocreales", family="Nectriaceae",
             status="uncalled", routing="lineage", n_loci=0)])
    _write(asco / "loci.tsv", gl, [
        dict(genome="A1", family_called="Ascomycota:MAT", contig="c1", start=1, end=9, idiomorph="MAT1-1", confidence="medium",
             locus_class="mat_locus", baseline_confidence="high", tier_note="vote_collapse", label_unsupported="yes",
             label_unsupported_reason="bit_margin<20;genus:Aureobasidium", defective_genes="MAT1-1-1"),
        dict(genome="A1", family_called="Ascomycota:MAT", contig="c2", start=1, end=9, idiomorph="MAT1-2", confidence="medium",
             locus_class="mat_locus", baseline_confidence="high", label_unsupported="no", high_basis=""),
        dict(genome="A1", family_called="Ascomycota:MAT", contig="c3", start=1, end=9, idiomorph="MAT1-2", confidence="high",
             locus_class="mat_locus", label_unsupported="no", high_basis_outside_scope="grey")])
    _write(basidio / "genomes.tsv", ASCO_G, [dict(genome="B1", species="Coprinopsis cinerea", class_="Agaricomycetes",
                                                  order="Agaricales", family="Psathyrellaceae", status="called", routing="lineage")])
    _write(basidio / "loci.tsv", gl, [dict(genome="B1", family_called="Basidiomycota:HD", contig="c1", start=1, end=9,
                                           idiomorph="undetermined", confidence="high", locus_class="mat_locus", label_unsupported="no")])
    if mucoro_tables:
        _write(mucoro / "genomes.tsv", ["source", "genome", "species", "class_", "order", "family", "status", "routing"], [
            dict(source="BFD", genome="M1", species="Mucor lusitanicus", order="Mucorales", status="called"),
            dict(source="LCG", genome="L1", species="Absidia sp.", order="", status="uncalled")])
        _write(mucoro / "loci.tsv", ["source", "genome", "contig", "start", "end", "idiomorph", "confidence", "locus_class"] + note_cols, [
            dict(source="BFD", genome="M1", contig="c9", start=5, end=50, idiomorph="Plus", confidence="high", locus_class="mat_locus",
                 label_unsupported="no", defective_genes="sexP")])
    else:
        _write(mucoro / "calls.tsv", ["source", "genome", "name", "bfd_order", "status", "idiomorph", "confidence", "locus_class",
                                      "contig", "locus_start", "locus_end"], [
            dict(source="BFD", genome="M1", name="Mucor lusitanicus", bfd_order="Mucorales", status="called", idiomorph="Plus",
                 confidence="high", locus_class="mat_locus", contig="c9", locus_start=5, locus_end=50),
            dict(source="LCG", genome="L1", name="Absidia sp.", bfd_order="", status="uncalled")])
    if manifest:
        for d in (asco, basidio, mucoro):
            (d / "MANIFEST.txt").write_text(manifest)
    return ["--asco", str(asco), "--basidio", str(basidio), "--mucoro", str(mucoro)]


def _genomes(path):
    with open(path, newline="") as fh:
        return {r["genome"]: r for r in csv.DictReader(fh, delimiter="\t")}


def test_default_release_is_the_v070_definitive_run():
    assert bsd.DEFAULT_RELEASE == "v070"
    assert bsd.RELEASES["v070"]["code"].startswith("cleanrun-2026-10c")
    for rel in bsd.RELEASES.values():
        assert {"code", "asco", "basidio", "mucoro"} <= set(rel)


def test_v060_tables_get_blank_note_columns(tmp_path):
    out = tmp_path / "out"
    bsd.main(["--release", "v060", "--out", str(out)] + _campaigns(tmp_path))
    taxa = json.loads((out / "taxa.json").read_text())
    assert taxa["fields"][:14] == ["dataset", "phylum", "class", "order", "family", "genus", "species", "cohort", "n", "called",
                                   "undet", "two", "alleles", "types"] and taxa["fields"][14] == "unsupported"
    assert all(r[14] == 0 for r in taxa["rows"])
    for g in _genomes(out / "genomes_bfd_asco.tsv").values():
        assert [g[c] for c in ("n_label_unsupported", "tier_notes", "high_basis", "loci_baseline_confidence", "locus_notes",
                               "n_defective_loci")] == [""] * 6
    assert _genomes(out / "genomes_bfd_asco.tsv")["A1"]["loci"].startswith("c1:1-9|MAT1-1|medium|mat_locus;c2:1-9|MAT1-2|medium|mat_locus")
    s = json.loads((out / "summary.json").read_text())
    assert s["release"]["name"] == "v060" and s["datasets"]["bfd_asco"]["n_label_unsupported"] == 0


def test_cleanrun_tables_pass_the_notes_and_the_baseline_tier_on(tmp_path):
    out = tmp_path / "out"
    bsd.main(["--release", "cleanrun", "--out", str(out)] + _campaigns(tmp_path, notes=True, mucoro_tables=True, manifest=MANIFEST))
    G = _genomes(out / "genomes_bfd_asco.tsv")
    g = G["A1"]
    assert g["n_label_unsupported"] == "1" and g["tier_notes"] == "vote_collapse" and g["high_basis"] == ""
    assert g["loci_baseline_confidence"] == "high;high;-"
    assert g["locus_notes"] == "vote_collapse,bit_margin<20+genus:Aureobasidium;-;-"   # outside-scope basis: counts only
    assert G["A2"]["n_label_unsupported"] == "0"           # uncalled: a count, like a called genome
    assert g["n_defective_loci"] == "1" and G["A2"]["n_defective_loci"] == "0"     # #151: loci with a defective gene
    d = json.loads((out / "summary.json").read_text())["datasets"]["bfd_asco"]
    assert d["n_label_unsupported"] == 1 and d["tier_note"] == {"vote_collapse": 1}
    assert d["tier_drop"] == 2 and d["tier_drop_without_note"] == 1   # the second drop has no listed mechanism
    assert d["code"] == "cleanrun-2026-10b (56e11cd)"
    assert d["defective_loci"] == 1                                # only A1's first locus
    others = json.loads((out / "summary.json").read_text())["datasets"]
    assert others["bfd_basidio"]["defective_loci"] == 0            # no defective locus
    assert others["bfd_mucoro"]["defective_loci"] == 1             # M1's locus (its sexP is defective in the new-table fixture)
    assert d["high_basis"] == {} and d["high_basis_outside_scope"] == {"grey": 1}
    taxa = json.loads((out / "taxa.json").read_text())
    assert sum(r[14] for r in taxa["rows"]) == 1
    assert taxa["fields"][15] == "high_basis_outside_scope" and [r[15] for r in taxa["rows"] if r[15]] == [{"grey": 1}]


def test_each_folder_carries_its_own_code_and_mixed_generations_are_refused(tmp_path):
    args = _campaigns(tmp_path, notes=True, manifest=MANIFEST)
    (tmp_path / "mucoro" / "MANIFEST.txt").unlink()          # a v0.6.0-style Mucoromycotina folder without a manifest
    with pytest.raises(SystemExit, match="mixes result generations"):
        bsd.main(["--release", "cleanrun", "--out", str(tmp_path / "o")] + args)
    out = tmp_path / "o2"
    bsd.main(["--release", "cleanrun", "--out", str(out), "--allow-mixed"] + args)
    s = json.loads((out / "summary.json").read_text())
    assert s["datasets"]["bfd_asco"]["code"] == "cleanrun-2026-10b (56e11cd)"
    assert s["datasets"]["zygolife_lcg"]["code"] == s["datasets"]["bfd_mucoro"]["code"] == "unrecorded"
    assert s["release"]["code"].startswith("mixed: ")


def test_folder_code_falls_back_to_the_release_that_names_the_folder(tmp_path):
    d = tmp_path / bsd.RELEASES["v060"]["mucoro"]
    d.mkdir()
    assert bsd.folder_code(str(d)) == bsd.RELEASES["v060"]["code"]
    assert bsd.folder_code(str(tmp_path)) == "unrecorded"


def test_held_out_sets_are_separate_and_out_of_the_headline(tmp_path):
    out = tmp_path / "out"
    bsd.main(["--out", str(out)] + _campaigns(tmp_path))
    s = json.loads((out / "summary.json").read_text())
    assert s["datasets"]["zygolife_lcg"]["held_out"] and s["datasets"]["zygolife_jena"]["held_out"]
    assert not s["datasets"]["bfd_mucoro"]["held_out"]
    assert s["headline"]["datasets"] == 3 and s["headline"]["genomes"] == 2 + 1 + 1   # A1, A2, B1, M1; not L1
    taxa = json.loads((out / "taxa.json").read_text())
    S = taxa["strings"]
    assert ("ZyGoLife", "LCG") in {(S[r[0]], S[r[7]]) for r in taxa["rows"]}


def test_summary_records_outside_folders_as_absolute_paths(tmp_path):
    out = tmp_path / "out"
    bsd.main(["--out", str(out)] + _campaigns(tmp_path))
    folders = json.loads((out / "summary.json").read_text())["release"]["folders"]
    assert not folders["asco"].startswith("..") and Path(folders["asco"]).is_absolute()
    assert bsd.show_path(str(_ROOT / "results" / "x")) == "results/x"


def test_mucoromycotina_new_tables_read_like_calls_tsv(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    bsd.main(["--out", str(a)] + _campaigns(tmp_path / "x"))
    bsd.main(["--out", str(b)] + _campaigns(tmp_path / "y", mucoro_tables=True))
    notes = ("n_label_unsupported", "tier_notes", "high_basis", "loci_baseline_confidence", "locus_notes", "n_defective_loci")
    for ds in ("bfd_mucoro", "zygolife_lcg"):
        ga, gb = _genomes(a / f"genomes_{ds}.tsv"), _genomes(b / f"genomes_{ds}.tsv")
        assert {k: {c: v for c, v in r.items() if c not in notes} for k, r in ga.items()} == \
               {k: {c: v for c, v in r.items() if c not in notes} for k, r in gb.items()}
    m1 = _genomes(b / "genomes_bfd_mucoro.tsv")["M1"]
    assert (m1["n_label_unsupported"], m1["loci_baseline_confidence"], m1["locus_notes"], m1["n_defective_loci"]) == ("0", "-", "-", "1")
    assert _genomes(b / "genomes_zygolife_lcg.tsv")["L1"]["n_label_unsupported"] == "0"


def test_missing_or_unset_folders_stop_the_build(tmp_path):
    args = _campaigns(tmp_path)
    with pytest.raises(SystemExit, match="does not exist"):
        bsd.main(["--out", str(tmp_path / "o"), "--asco", str(tmp_path / "nope")] + args[2:])
    with pytest.raises(SystemExit, match="no mucoro folder"):
        bsd.resolve_release("cleanrun", asco=args[1], basidio=args[3])


def test_a_cleanrun_table_from_before_the_defective_column_gives_zero_not_an_error(tmp_path):
    out = tmp_path / "out"
    bsd.main(["--release", "cleanrun", "--out", str(out)] + _campaigns(tmp_path, notes=True, mucoro_tables=True, manifest=MANIFEST,
                                                                      defective_col=False))
    assert {g["n_defective_loci"] for g in _genomes(out / "genomes_bfd_asco.tsv").values()} == {"0"}
    s = json.loads((out / "summary.json").read_text())["datasets"]
    assert all(d["defective_loci"] == 0 for d in s.values())
    assert _genomes(out / "genomes_bfd_asco.tsv")["A1"]["n_label_unsupported"] == "1"      # the other note columns are unaffected


# ------------------------------------------------------------------------------------------ shards

def _tar_zst(path, members):
    path.parent.mkdir(parents=True, exist_ok=True)
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w") as tf:
        for name, rep in members:
            data = yaml.safe_dump(rep).encode()
            ti = tarfile.TarInfo(name)
            ti.size = len(data)
            tf.addfile(ti, io.BytesIO(data))
    path.write_bytes(subprocess.run(["zstd", "-q", "-c"], input=buf.getvalue(), capture_output=True, check=True).stdout)


GENE = dict(gene="MAT1-1-1", start=2, end=8, strand="+", role="core_MAT", status="polished_agree")
REP = dict(detected=[dict(contig="c1", start=1, end=9, idiomorph="MAT1-1", confidence="medium", locus_class="mat_locus",
                          gene_evidence=[GENE]),
                     dict(contig="c2", start=1, end=9, idiomorph="MAT1-1", confidence="high", locus_class="mat_locus")])
needs_zstd = pytest.mark.skipif(shutil.which("zstd") is None, reason="needs the zstd command")


@needs_zstd
def test_locus_shards_append_the_notes_in_report_order(tmp_path):
    bls = _load("build_locus_shards")
    camp = tmp_path / "asco"
    _write(camp / "loci.tsv", ASCO_L + NOTES, [
        dict(genome="G1", contig="c1", start=1, end=9, baseline_confidence="high", tier_note="margin_fallback", label_unsupported="no",
             defective_genes="MAT1-1-1"),
        dict(genome="G1", contig="c2", start=1, end=9, label_unsupported="yes", label_unsupported_reason="genus:Xylaria")])
    _tar_zst(camp / "reports_all.tar.zst", [("./wave_0/small_000/runs/G1/detection_report.yaml", REP)])
    shards, noted = bls.campaign_shards("bfd_asco", str(camp))
    assert noted
    assert shards["bfd_asco"]["G1"] == [
        ["c1", 1, 9, "MAT1-1", "medium", "mat_locus", [["MAT1-1-1", 2, 8, "+", "core_MAT", "polished_agree"]],
         "high", "margin_fallback", "", "", "MAT1-1-1"],
        ["c2", 1, 9, "MAT1-1", "high", "mat_locus", [], "", "", "genus:Xylaria", "", ""]]
    assert bls.NOTE_FIELDS[-1] == "defective_genes"           # appended last: a viewer that reads fields 0-6 is unaffected


@needs_zstd
def test_locus_shards_from_a_table_that_predates_the_defective_column_get_a_blank_field(tmp_path):
    bls = _load("build_locus_shards")
    camp = tmp_path / "asco"
    _write(camp / "loci.tsv", ASCO_L + [c for c in NOTES if c != "defective_genes"], [
        dict(genome="G1", contig="c1", start=1, end=9, baseline_confidence="high", label_unsupported="no"),
        dict(genome="G1", contig="c2", start=1, end=9, label_unsupported="no")])
    _tar_zst(camp / "reports_all.tar.zst", [("./wave_0/small_000/runs/G1/detection_report.yaml", REP)])
    shards, noted = bls.campaign_shards("bfd_asco", str(camp))
    assert noted and [row[-1] for row in shards["bfd_asco"]["G1"]] == ["", ""]


@needs_zstd
def test_locus_shards_choose_the_same_run_as_the_tables_for_a_duplicated_genome(tmp_path):
    bls = _load("build_locus_shards")
    bct = _load("build_cleanrun_tables")
    camp = tmp_path / "basidio"
    _tar_zst(camp / "reports_all.tar.zst", [
        ("./runs/small_000/runs/G1/detection_report.yaml", REP),             # an earlier chunk called it
        ("./runs/zz_rerun/runs/G1/detection_report.yaml", dict(detected=[]))])  # the re-run is uncalled
    shards, _ = bls.campaign_shards("bfd_basidio", str(camp))
    runs, _ = bct.read_campaign(str(camp / "reports_all.tar.zst"))
    picked, _ = bct.pick_runs(runs)
    assert picked["G1"]["summary"]["loci"] == [] and "G1" not in shards.get("bfd_basidio", {})


@needs_zstd
def test_locus_shards_stop_when_loci_tsv_disagrees_with_the_reports(tmp_path):
    bls = _load("build_locus_shards")
    camp = tmp_path / "asco"
    _write(camp / "loci.tsv", ASCO_L + NOTES, [dict(genome="G1", contig="c1", start=1, end=9, label_unsupported="no")])
    _tar_zst(camp / "reports_all.tar.zst", [("./wave_0/small_000/runs/G1/detection_report.yaml", REP)])
    with pytest.raises(SystemExit, match="disagree on the number of loci"):
        bls.campaign_shards("bfd_asco", str(camp))


def test_locus_rows_mark_the_winner_and_the_gene_set_aside_in_a_resolved_pair():
    bls = _load("build_locus_shards")
    locus = dict(contig="c1", start=1, end=99, idiomorph="Plus", confidence="high", locus_class="mat_locus", gene_evidence=[
        dict(gene="sexP", start=10, end=20, strand="-", role="core_MAT", status="polished_disagree", lost_to=None),
        dict(gene="sexM", start=12, end=18, strand="-", role="core_MAT", status="polished_single", lost_to="sexP"),
        dict(gene="rnhA", start=1, end=5, strand="-", role="flanking_conserved", status="polished_agree", lost_to=None)],
        idiomorph_resolutions=[dict(winner="sexP", loser="sexM", basis="hmm_bitscore", winner_hmm_score=364.33, loser_hmm_score=70.79,
                                    winner_identity=44.9, loser_identity=34.1, winner_span=[10, 20], loser_span=[12, 18])])
    genes = bls.locus_row(locus)[6]
    assert genes[0][6:] == [["W", "sexM", "hmm_bitscore", "364 vs 71"]]
    assert genes[1][6:] == [["L", "sexP", "hmm_bitscore", "364 vs 71"]]
    assert len(genes[2]) == 6                                 # a gene outside a resolved pair keeps the six-item form
    r = locus["idiomorph_resolutions"][0]
    r.update(basis="first_pass_identity", winner_hmm_score=None, loser_hmm_score=None)
    assert bls.locus_row(locus)[6][1][6] == ["L", "sexP", "first_pass_identity", "44.9% vs 34.1%"]
    r.update(winner_span=[10, 21])                            # a span that matches no gene leaves the pair unmarked
    assert len(bls.locus_row(locus)[6][0]) == 6


def test_locus_rows_without_resolutions_are_unchanged():
    bls = _load("build_locus_shards")
    locus = dict(contig="c1", start=1, end=9, idiomorph="MAT1-1", confidence="high", locus_class="mat_locus", gene_evidence=[
        dict(gene="MAT1-1-1", start=2, end=8, strand="+", role="core_MAT", status="polished_agree", lost_to=None)])
    assert bls.locus_row(locus)[6] == [["MAT1-1-1", 2, 8, "+", "core_MAT", "polished_agree"]]


def test_mucoromycotina_family_reaches_the_taxa_rows_and_the_genome_table(tmp_path):
    out = tmp_path / "out"
    camps = _campaigns(tmp_path, notes=True, mucoro_tables=True, manifest=MANIFEST)
    mu = Path(camps[camps.index("--mucoro") + 1])
    g = mu / "genomes.tsv"
    text = g.read_text().splitlines()
    cols = text[0].split("\t")
    assert "family" in cols
    fi = cols.index("family")
    text[1:] = ["\t".join(f if i != fi else "Mucoraceae" for i, f in enumerate(l.split("\t"))) for l in text[1:]]
    g.write_text("\n".join(text) + "\n")
    bsd.main(["--release", "cleanrun", "--out", str(out)] + camps)
    taxa = json.loads((out / "taxa.json").read_text())
    i = {k: n for n, k in enumerate(taxa["fields"])}
    fam = {taxa["strings"][r[i["family"]]] for r in taxa["rows"] if taxa["strings"][r[i["phylum"]]] == "Mucoromycota"}
    assert fam == {"Mucoraceae"}


def test_idiomorph_labels_sit_beside_the_raw_names(tmp_path):
    """Labels from the registry: taxa field `label_alleles`, genome columns `call_label` / `alleles_label`, and the tooltip file; the raw columns are unchanged."""
    out = tmp_path / "out"
    bsd.main(["--release", "v060", "--out", str(out)] + _campaigns(tmp_path))
    taxa = json.loads((out / "taxa.json").read_text())
    assert taxa["fields"][16] == "label_alleles"
    S = taxa["strings"]
    # the genome totals agree between the raw and the label counts in every taxon row
    for r in taxa["rows"]:
        assert sum(r[12].values()) == sum(r[16].values())
    labels = json.loads((out / "idiomorph_labels.json").read_text())["labels"]
    assert labels["Peziz__BipMay_MAT1-1"]["raw"] == ["MAT1-1 (Ascomycota:MAT)"]
    assert labels["Mucor__PhyBla_plus"]["reference"] == "Phycomyces blakesleeanus"
    G = _genomes(out / "genomes_bfd_asco.tsv")
    assert G["A1"]["call"] == "MAT1-1+MAT1-2" and G["A1"]["call_label"] == "Peziz__BipMay_MAT1-1+Peziz__BipMay_MAT1-2"
    assert G["A1"]["alleles_label"] == "Peziz__BipMay_MAT1-1,Peziz__BipMay_MAT1-2"
    assert G["A2"]["call"] == G["A2"]["call_label"] == "uncalled"
    assert _genomes(out / "genomes_bfd_basidio.tsv")["B1"]["call_label"] == "undetermined"     # undetermined idiomorph: no allele to label
    assert any(S[int(k)] == "two:Peziz__BipMay_MAT1-1+Peziz__BipMay_MAT1-2" for r in taxa["rows"] for k in r[16])      # A1 carries both idiomorphs
