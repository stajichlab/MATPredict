"""scripts/build_cleanrun_tables.py: post-run tables and the curator's reporting columns (tier_note,
label_unsupported, high_basis), on small synthetic report trees."""
import csv
import importlib.util
import os
import tarfile
from pathlib import Path

import yaml

_ROOT = Path(__file__).resolve().parents[1]


def _load(name):
    spec = importlib.util.spec_from_file_location(name, _ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


bct = _load("build_cleanrun_tables")

POL = "polished_agree"


def ge(gene, status=POL, identity=60.0, coverage=80.0, lost_to=None, role="core_MAT"):
    return dict(gene=gene, role=role, status=status, identity=identity, coverage=coverage, lost_to=lost_to,
                contig="c1", start=100, end=400, strand="+", bitscore=100.0)


def det(idio="MAT1-2", conf="high", cls="mat_locus", margin=None, cands=1, genes=(), polished=None, ignored=(),
        family="Ascomycota:MAT", contig="c1", start=1, end=5000, homothallic_vote=None, defective=()):
    genes = list(genes)
    d = dict(family=family, contig=contig, start=start, end=end, confidence=conf, idiomorph=idio,
             idiomorph_candidates=[dict(idiomorph=f"I{i}", score=50.0 - i) for i in range(cands)],
             detection_pass="strict", locus_class=cls, verification=None,
             polished_genes=sum(g["status"] == POL for g in genes) if polished is None else polished,
             idiomorph_margin=margin, genes_found=[g["gene"] for g in genes], gene_evidence=genes,
             confidence_ignored_genes=list(ignored))
    if homothallic_vote:
        d["homothallic_note"] = dict(vote_idiomorph=homothallic_vote)
    if defective:                       # the report writes the key only when a gene is defective (issue #151)
        d["defective_genes"] = sorted(defective)
    return d


def L(**kw):
    """The compact locus summary the rules work on."""
    return bct.summarise_locus(det(**kw))


# ------------------------------------------------------------------------------------------ margin basis


def test_margin_basis_bits_needs_two_candidates_identity_fallback_otherwise():
    assert bct.margin_basis(L(margin=12.0, cands=2)) == "bits"
    assert bct.margin_basis(L(margin=12.0, cands=1)) == "identity_fallback"
    assert bct.margin_basis(L(margin=12.0, cands=0)) == "identity_fallback"
    assert bct.margin_basis(L(margin=None, cands=2)) == ""


# ------------------------------------------------------------------------------------------ label_unsupported


def test_label_unsupported_marks_a_bit_margin_under_20_only():
    assert bct.label_unsupported(L(margin=19.9, cands=2), "Fusarium") == "bit_margin<20"
    assert bct.label_unsupported(L(margin=20.0, cands=2), "Fusarium") == ""


def test_label_unsupported_does_not_mark_the_identity_fallback():
    # #145, not #138: a single-candidate locus' margin is an identity gap, not bits
    assert bct.label_unsupported(L(margin=1.0, cands=1), "Beauveria") == ""
    assert bct.label_unsupported(L(margin=None, cands=0), "Beauveria") == ""


def test_label_unsupported_never_marks_an_undetermined_locus():
    # curator 2026-10-24: a homothallic candidate publishes no MAT1-1/MAT1-2 label, so neither rule applies
    homo = L(idio="undetermined", cls="homothallic_candidate", margin=3.0, cands=2, homothallic_vote="MAT1-1")
    assert bct.label_unsupported(homo, "Neurospora") == ""
    assert bct.label_unsupported(homo, "Aureobasidium") == ""
    assert bct.label_unsupported(L(idio="A", family="Ascomycota:MTL", margin=3.0, cands=2), "Candida") == "bit_margin<20"


def test_label_unsupported_genus_rule_marks_every_locus_of_the_three_genera():
    for g in ("Aureobasidium", "Xylaria", "Friedmanniomyces"):
        assert bct.label_unsupported(L(margin=300.0, cands=2), g) == f"genus:{g}"
        assert bct.label_unsupported(L(margin=None, cands=0), g) == f"genus:{g}"
    assert bct.label_unsupported(L(margin=3.0, cands=2), "Aureobasidium") == "bit_margin<20;genus:Aureobasidium"


# ------------------------------------------------------------------------------------------ tier_note


def test_tier_note_blank_when_the_tier_did_not_drop():
    old = L(conf="medium", margin=1.0)
    assert bct.tier_note(L(conf="medium", margin=1.0), old) == ("", False)
    assert bct.tier_note(L(conf="high", margin=1.0), old) == ("", False)
    assert bct.tier_note(L(conf="low"), None) == ("", False)


def test_tier_note_margin_fallback():
    old = L(conf="high", margin=5.18, cands=1, genes=[ge("MAT1-2-1"), ge("APN2")])
    new = L(conf="medium", margin=4.96, cands=1, genes=[ge("MAT1-2-1"), ge("APN2")])
    assert bct.tier_note(new, old) == ("margin_fallback", True)


def test_tier_note_vote_collapse():
    old = L(conf="high", margin=12.0, cands=2)
    new = L(conf="medium", margin=2.0, cands=2)
    assert bct.tier_note(new, old) == ("vote_collapse", True)


def test_tier_note_polished_model_lost():
    old = L(idio="A", family="Ascomycota:MTL", conf="high", genes=[ge("MTLA1"), ge("MTLA2")])
    new = L(idio="A", family="Ascomycota:MTL", conf="medium", genes=[ge("MTLA1", status="unpolished"), ge("MTLA2")])
    assert bct.tier_note(new, old) == ("polished_model_lost", True)


def test_tier_note_distant_fragment_is_an_unpolished_gene_no_longer_ignored():
    base = [ge("MAT1-2-1"), ge("APN2"), ge("COX13")]
    frag = ge("MAT1-1-2", status="unpolished", identity=40.0, coverage=11.0)
    old = L(conf="high", margin=33.0, cands=2, genes=base + [frag], ignored=["MAT1-1-2"])
    new = L(conf="medium", margin=36.6, cands=2, genes=base + [frag])
    assert bct.tier_note(new, old) == ("distant_fragment", True)


def test_tier_note_blank_when_the_locus_class_also_changed():
    old = L(conf="high", cls="mat_locus", margin=1.0)
    new = L(conf="medium", cls="partial_locus", margin=1.0)
    assert bct.tier_note(new, old) == ("", True)


def test_tier_note_blank_when_no_mechanism_applies():
    assert bct.tier_note(L(conf="medium", margin=50.0, cands=2), L(conf="high", margin=50.0, cands=2)) == ("", True)


# ------------------------------------------------------------------------------------------ high_basis

SCOPE = ("Pleosporaceae",)


def test_high_basis_supported_by_a_polished_extra_gene():
    loc = L(idio="MAT1-1", genes=[ge("MAT1-1-1"), ge("GAP1"), ge("APN2", role="flanking_conserved")])
    assert bct.high_basis(loc, "Pleosporaceae", SCOPE) == "supported"


def test_high_basis_supported_by_an_unpolished_gene_at_50_identity_30_coverage():
    loc = L(idio="MAT1-1", genes=[ge("MAT1-1-1"), ge("SLA2", status="unpolished", identity=50.0, coverage=30.0)])
    assert bct.high_basis(loc, "Pleosporaceae", SCOPE) == "supported"


def test_high_basis_tie_break_when_extras_are_weak_or_the_opposite_core():
    loc = L(idio="MAT1-1", genes=[ge("MAT1-1-1"), ge("GAP1"), ge("BGL1"),
                                   ge("SLA2", status="unpolished", identity=70.0, coverage=29.0),
                                   ge("MAT1-2-1")])   # polished, but the opposite idiomorph's core gene
    assert bct.high_basis(loc, "Pleosporaceae", SCOPE) == "tie_break"
    assert bct.high_basis(L(idio="MAT1-2", genes=[ge("MAT1-2-1"), ge("GAP1")]), "Pleosporaceae", SCOPE) == "tie_break"


def test_high_basis_lost_to_the_called_core_counts_as_opposite():
    loc = L(idio="MAT1-2", genes=[ge("MAT1-2-1"), ge("MAT1-1-3", status="not_polish_candidate", identity=45.0,
                                                      coverage=40.0, lost_to="MAT1-2-1")])
    assert bct.high_basis(loc, "Pleosporaceae", SCOPE) == "tie_break"


def test_high_basis_grey():
    loc = L(idio="MAT1-2", genes=[ge("MAT1-2-1"), ge("MAT1-1-3", status="unpolished", identity=40.0, coverage=40.0)])
    assert bct.high_basis(loc, "Pleosporaceae", SCOPE) == "grey"


def test_high_basis_blank_outside_scope_below_high_or_undetermined():
    genes = [ge("MAT1-1-1"), ge("APN2")]
    assert bct.high_basis(L(idio="MAT1-1", genes=genes), "Nectriaceae", SCOPE) == ""
    assert bct.high_basis(L(idio="MAT1-1", genes=genes), "Nectriaceae", None) == "supported"   # scope 'all'
    assert bct.high_basis(L(idio="MAT1-1", conf="medium", genes=genes), "Pleosporaceae", SCOPE) == ""
    assert bct.high_basis(L(idio="undetermined", cls="homothallic_candidate", genes=genes, homothallic_vote="MAT1-1"),
                          "Pleosporaceae", SCOPE) == ""
    assert bct.high_basis(L(idio="MAT1-1", family="Ascomycota:MATtub", genes=genes), "Pleosporaceae", SCOPE) == ""


# ------------------------------------------------------------------------------------------ matching, categories


def test_match_baseline_by_contig_overlap():
    old = [L(contig="c1", start=1, end=100), L(contig="c1", start=50, end=5000), L(contig="c2", start=1, end=5000)]
    new = [L(contig="c1", start=10, end=4000), L(contig="c3", start=1, end=10)]
    assert bct.match_baseline(new, old) == [1, None]


def test_match_baseline_prefers_the_same_family_over_a_larger_overlap():
    old = [L(contig="c1", start=1, end=5000, family="Ascomycota:MATsc"), L(contig="c1", start=1000, end=2000)]
    new = [L(contig="c1", start=1, end=5000)]
    assert bct.match_baseline(new, old) == [1]


def test_match_baseline_skips_a_null_span():
    old = [L(contig="c1", start=None, end=None), L(contig="c1", start=1, end=50)]
    assert bct.match_baseline([L(contig="c1", start=None, end=None), L(contig="c1", start=1, end=50)], old) == [None, 1]


def test_choose_run_prefers_a_report_then_the_last_path():
    assert bct.choose_run([("a/runs/G", True), ("b/runs/G", False)]) == "a/runs/G"
    assert bct.choose_run([("a/runs/G", True), ("b/runs/G", True)]) == "b/runs/G"


def test_change_category():
    def summ(*loci):
        return {"loci": list(loci)}
    a = L(conf="high")
    assert bct.change_category(None, summ(a)) == "no new report"
    assert bct.change_category(summ(a), None) == "NEW call"
    assert bct.change_category(summ(), summ()) == "uncalled both"
    assert bct.change_category(summ(L(conf="medium")), summ(a)) == "confidence DOWN"
    assert bct.change_category(summ(L(idio="undetermined")), summ(a)) == "idiomorph -> undetermined"
    assert bct.change_category(summ(L(start=600)), summ(a)) == "identical call"
    assert bct.change_category(summ(L(start=1200)), summ(a)) == "span changed >1 kb"


# ------------------------------------------------------------------------------------------ end to end


def _report(root, rel, detected):
    d = root / rel
    d.mkdir(parents=True)
    (d / "detection_report.yaml").write_text(yaml.safe_dump(dict(routing_mode="lineage", detected=list(detected))))
    (d / "wall_seconds").write_text("42\n")


def _tree(tmp_path):
    new, old = tmp_path / "new", tmp_path / "old"
    (new / "lists").mkdir(parents=True)
    (new / "lists" / "wave_0.tsv").write_text("G1\t1\nG2\t2\nG3\t3\nG4\t4\nG5\t5\n")
    # G1 Macrophomina-like identity fallback drop; G2 Aureobasidium vote collapse; G3 no new report but a failure
    _report(new, "wave_0/small_000/runs/G1", [det(conf="medium", margin=4.96, cands=1, genes=[ge("MAT1-2-1")])])
    _report(old, "wave_0/runs/G1", [det(conf="high", margin=5.18, cands=1, genes=[ge("MAT1-2-1")])])
    _report(new, "wave_0/small_001/runs/G2", [det(idio="MAT1-1", conf="medium", margin=2.0, cands=2, genes=[ge("MAT1-1-1")])])
    _report(old, "wave_0/runs/G2", [det(idio="MAT1-2", conf="high", margin=9.0, cands=2, genes=[ge("MAT1-1-1")])])
    (new / "wave_0/small_001/runs/G3").mkdir(parents=True)
    (new / "wave_0/small_001/runs/G3/stderr.log").write_text("boom\n")
    _report(new, "wave_0/small_002/runs/G4", [det(idio="MAT1-1", genes=[ge("MAT1-1-1"), ge("GAP1")], defective=["GAP1", "MAT1-1-1"])])
    _report(new, "wave_0/small_002/runs/UNLISTED", [])
    # G5 outside the Pleosporaceae scope: a high call resting on the opposite core gene only
    _report(new, "wave_0/small_002/runs/G5", [det(idio="MAT1-2", genes=[ge("MAT1-2-1"), ge("MAT1-1-1")])])
    samples = tmp_path / "samples.csv"
    with open(samples, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["ASMID", "SPECIES", "CLASS", "ORDER", "FAMILY", "GENUS"])
        w.writerow(["G1", "Macrophomina phaseolina", "Dothideomycetes", "Botryosphaeriales", "Botryosphaeriaceae", "Macrophomina"])
        w.writerow(["G2", "Aureobasidium pullulans", "Dothideomycetes", "Dothideales", "Saccotheciaceae", "Aureobasidium"])
        w.writerow(["G4", "Curvularia lunata", "Dothideomycetes", "Pleosporales", "Pleosporaceae", "Curvularia"])
        w.writerow(["G5", "Cercospora beticola", "Dothideomycetes", "Mycosphaerellales", "Mycosphaerellaceae", "Cercospora"])
    return new, old, samples


def _rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_end_to_end_tables(tmp_path):
    new, old, samples = _tree(tmp_path)
    out = tmp_path / "out"
    bct.main(["--campaign", str(new), "--baseline", str(old), "--samples", str(samples), "--out", str(out)])
    G = {r["genome"]: r for r in _rows(out / "genomes.tsv")}
    assert list(G) == ["G1", "G2", "G3", "G4", "G5"]          # list order; the unlisted report is left out
    assert G["G1"]["change_vs_baseline"] == "confidence DOWN" and G["G1"]["n_tier_note"] == "1" and G["G1"]["n_tier_drop"] == "1"
    assert G["G3"]["n_label_unsupported"] == "0"
    assert G["G2"]["change_vs_baseline"] == "idiomorph flip"
    assert G["G3"]["status"] == "failed" and G["G3"]["baseline_status"] == "no_report"
    assert G["G4"]["change_vs_baseline"] == "NEW call" and G["G4"]["high_basis"] == "tie_break"
    loci = {r["genome"]: r for r in _rows(out / "loci.tsv")}
    assert loci["G1"]["tier_note"] == "margin_fallback" and loci["G1"]["margin_basis"] == "identity_fallback"
    assert loci["G1"]["label_unsupported"] == "no" and loci["G1"]["label_unsupported_reason"] == "" and loci["G1"]["baseline_confidence"] == "high"
    assert loci["G2"]["tier_note"] == "vote_collapse"
    assert loci["G2"]["label_unsupported"] == "yes"
    assert loci["G2"]["label_unsupported_reason"] == "bit_margin<20;genus:Aureobasidium"
    assert loci["G2"]["idiomorph"] == "MAT1-1" and loci["G2"]["confidence"] == "medium"   # the report's values are kept
    assert loci["G4"]["high_basis"] == "tie_break" and loci["G4"]["high_basis_outside_scope"] == ""
    assert loci["G5"]["high_basis"] == "" and loci["G5"]["high_basis_outside_scope"] == "tie_break"
    # defective genes (#151): the locus lists them `|`-joined, the genome counts loci, an intact locus stays blank/zero
    assert loci["G4"]["defective_genes"] == "GAP1|MAT1-1-1" and G["G4"]["n_defective_loci"] == "1"
    assert all(loci[g]["defective_genes"] == "" for g in ("G1", "G2", "G5"))
    assert [G[g]["n_defective_loci"] for g in ("G1", "G2", "G3", "G5")] == ["0", "0", "0", "0"]    # G3 failed: a count of 0
    by = _rows(out / "by_class.tsv")
    assert by[0]["defective_loci"] == "1"
    assert (by[0]["high_tie_break"], by[0]["outside_scope_high_tie_break"]) == ("1", "1")
    assert by[0]["class"] == "Dothideomycetes" and by[0]["margin_fallback"] == "1" and by[0]["vote_collapse"] == "1"


def test_by_rank_counts_loci_with_a_defective_gene_not_genes():
    base = dict(baseline_confidence="", confidence="high", tier_note="", high_basis="", high_basis_outside_scope="",
                label_unsupported="no")
    G = [dict(class_="C1", status="called", wall_s="5", routing="lineage"),
         dict(class_="C1", status="uncalled", wall_s="7", routing="lineage"),
         dict(class_="C2", status="called", wall_s="9", routing="lineage")]
    Lr = [dict(base, class_="C1", defective_genes="MAT1-1-1|GAP1"),      # one locus with two defective genes counts once
          dict(base, class_="C1", defective_genes="PIK1"),
          dict(base, class_="C1", defective_genes=""),
          dict(base, class_="C2", defective_genes="")]
    head, rows = bct.by_rank(G, Lr, "class_")
    assert head[-1] == "defective_loci"
    assert {r[0]: r[-1] for r in rows} == {"C1": 2, "C2": 0}


def test_output_is_deterministic_and_a_tarball_reads_like_the_folder(tmp_path):
    new, old, samples = _tree(tmp_path)
    tar = tmp_path / "old.tar.gz"
    with tarfile.open(tar, "w:gz") as tf:
        tf.add(old, arcname=".")
    a, b = tmp_path / "a", tmp_path / "b"
    bct.main(["--campaign", str(new), "--baseline", str(old), "--samples", str(samples), "--out", str(a)])
    bct.main(["--campaign", str(new), "--baseline", str(tar), "--samples", str(samples), "--out", str(b)])
    for f in ("genomes.tsv", "loci.tsv", "by_class.tsv"):
        assert (a / f).read_bytes() == (b / f).read_bytes()


def test_basidiomycota_layout_and_duplicate_run_folders(tmp_path):
    new = tmp_path / "basidio"
    (new / "lists").mkdir(parents=True)
    (new / "lists" / "small_000.tsv").write_text("B1\t5\n")
    _report(new, "runs/small_000/runs/B1", [])
    _report(new, "runs/zz_rerun/runs/B1", [det(idio="undetermined", family="Basidiomycota:HD", genes=[ge("HD1")])])
    out = tmp_path / "out"
    bct.main(["--campaign", str(new), "--out", str(out), "--rank", "order"])
    G = _rows(out / "genomes.tsv")
    assert [(r["genome"], r["status"], r["baseline_status"], r["change_vs_baseline"]) for r in G] == [("B1", "called", "", "")]
    assert (out / "by_order.tsv").exists()


def test_mucoromycotina_layout(tmp_path):
    new = tmp_path / "mucoro"
    new.mkdir()
    (new / "inputs.tsv").write_text("source\tid\tgenome\torder\tspecies\ttaxid\n"
                                    "BFD\tGCA_1\t/x.fa\tMucorales\tMucor lusitanicus\t1\nLCG\tAbsidia_x_1\t/y.fa\t\t\t\n")
    _report(new, "runs/BFD__GCA_1", [det(idio="Plus", family="Mucoromycota:MAT", margin=40.0, cands=2, genes=[ge("sexP")],
                                         defective=["sexP"])])
    _report(new, "runs/LCG__Absidia_x_1", [])
    calls = tmp_path / "calls.tsv"
    calls.write_text("source\tgenome\tname\tbfd_order\nLCG\tAbsidia_x_1\tAbsidia x NRRL 1\t\n")
    out = tmp_path / "out"
    bct.main(["--campaign", str(new), "--out", str(out), "--layout", "mucoro", "--baseline-calls", str(calls)])
    G = _rows(out / "genomes.tsv")
    assert [(r["source"], r["genome"], r["species"], r["class_"], r["order"], r["status"]) for r in G] == [
        ("BFD", "GCA_1", "Mucor lusitanicus", "Mucoromycotina", "Mucorales", "called"),
        ("LCG", "Absidia_x_1", "Absidia x NRRL 1", "Mucoromycotina", "Mucorales", "uncalled")]    # order from the genus lineage
    assert [(r["genus"], r["family"]) for r in G] == [("Mucor", "Mucoraceae"), ("Absidia", "Cunninghamellaceae")]
    loci = _rows(out / "loci.tsv")
    assert loci[0]["source"] == "BFD" and loci[0]["margin_basis"] == "bits" and loci[0]["label_unsupported"] == "no"
    assert loci[0]["defective_genes"] == "sexP" and [r["n_defective_loci"] for r in G] == ["1", "0"]     # #151 in the mucoro layout


def test_mucoromycotina_genus_lineage_fills_blanks_and_never_overwrites(tmp_path):
    new = tmp_path / "mucoro"
    new.mkdir()
    (new / "inputs.tsv").write_text("source\tid\tgenome\torder\tspecies\ttaxid\n"
                                    "BFD\tU1\t/x.fa\tUmbelopsidales\tUmbelopsis ramanniana\t1\n"
                                    "LCG\tM1\t/y.fa\t\tMucor circinelloides\t\nLCG\tQ1\t/z.fa\t\tUnknownus sp\t\n")
    for g in ("BFD__U1", "LCG__M1", "LCG__Q1"):
        _report(new, f"runs/{g}", [])
    table = tmp_path / "lin.tsv"
    table.write_text("genus\ttaxid\tclass\torder\tfamily\nUmbelopsis\t1\tMucoromycetes\tMucorales\tWrongaceae\n"
                     "Mucor\t2\tMucoromycetes\tMucorales\tMucoraceae\n")
    out = tmp_path / "out"
    bct.main(["--campaign", str(new), "--out", str(out), "--layout", "mucoro", "--genus-lineage", str(table)])
    got = {r["genome"]: (r["order"], r["family"], r["genus"]) for r in _rows(out / "genomes.tsv")}
    assert got["U1"] == ("Umbelopsidales", "Wrongaceae", "Umbelopsis")     # a given order is kept; a blank family is filled
    assert got["M1"] == ("Mucorales", "Mucoraceae", "Mucor")
    assert got["Q1"] == ("", "", "Unknownus")                               # a genus outside the table stays blank
    bct.main(["--campaign", str(new), "--out", str(tmp_path / "out2"), "--layout", "mucoro", "--genus-lineage", str(tmp_path / "none.tsv")])
    assert {r["genome"]: (r["order"], r["family"]) for r in _rows(tmp_path / "out2" / "genomes.tsv")}["M1"] == ("", "")   # no table: as before


def test_the_committed_genus_lineage_table_covers_orders_and_families():
    rows = bct.load_genus_lineage(os.path.join(os.path.dirname(bct.__file__), "mucoromycotina_genus_lineage.tsv"))
    assert len(rows) == 51 and all(o and f for o, f in rows.values())
    assert rows["Thamnidium"] == ("Mucorales", "Mucoraceae") and rows["Umbelopsis"][0] == "Umbelopsidales"
