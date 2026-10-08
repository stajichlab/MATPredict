"""`matpredict report genome`: the HTML page renders every report shape and escapes its input."""
from __future__ import annotations

from html.parser import HTMLParser
from pathlib import Path

import pytest
import yaml

from MATPredict.__main__ import main
from MATPredict.report import pdf
from MATPredict.report.html import _css_string, render_genome_report
from MATPredict.report.cli import sample_name
from MATPredict.report.svg import assign_lanes, contig_end

FIXTURES = Path(__file__).parent / "fixtures"
CASES = sorted(p.stem for p in FIXTURES.glob("*.yaml"))


def _doc(name: str) -> dict:
    return yaml.safe_load((FIXTURES / f"{name}.yaml").read_text())


class _Collector(HTMLParser):
    """Element names and the document's text, so tests check what a browser would show."""

    def __init__(self) -> None:
        super().__init__()
        self.tags: list[str] = []
        self.text: list[str] = []
        self._skip = 0

    def handle_starttag(self, tag, attrs):
        self.tags.append(tag)
        if tag in ("style", "script"):
            self._skip += 1

    def handle_endtag(self, tag):
        if tag in ("style", "script"):
            self._skip -= 1

    def handle_data(self, data):
        if not self._skip:
            self.text.append(data)


def _parse(html: str) -> _Collector:
    c = _Collector()
    c.feed(html)
    c.close()
    return c


@pytest.mark.parametrize("case", CASES)
def test_every_fixture_renders_and_shows_every_gene(case):
    doc = _doc(case)
    html = render_genome_report(doc, generated="2026-10-07 00:00 UTC")
    parsed = _parse(html)
    text = " ".join(parsed.text)
    assert html.startswith("<!doctype html>")
    for call in doc.get("detected") or []:
        assert call["idiomorph"] in text
        for g in call.get("gene_evidence") or []:
            assert g["gene"] in text
        for name in call.get("genes_missing") or []:
            assert name in text
    # One figure per contig a call's genes lie on.
    n_contigs = sum(len({g["contig"] for g in c.get("gene_evidence") or []}) or 1 for c in doc.get("detected") or [])
    assert parsed.tags.count("svg") == n_contigs


def test_result_sentence_for_each_outcome():
    assert "Not searched" in render_genome_report(_doc("not_searched"))
    none = render_genome_report(_doc("none_called_gap"))
    assert "No MAT locus called" in none and "Assembly gap" in none
    two = render_genome_report(_doc("two_idiomorphs"))
    assert "Two idiomorphs" in two and "not a verdict" in two
    basid = render_genome_report(_doc("basidio_hd_pr"))
    assert "unconfirmed" in basid and "receptor array" in basid.lower() and "Subloci" in basid


def test_hostile_strings_are_escaped():
    html = render_genome_report(_doc("hostile_strings"))
    parsed = _parse(html)
    assert parsed.tags.count("script") == 1  # only the page's own print helper
    assert "img" not in parsed.tags
    assert "<script>alert" not in html
    assert "<img src=x" not in html


def test_css_string_cannot_close_the_style_element():
    s = _css_string('</style><script>alert("x")</script>')
    assert "<" not in s and ">" not in s
    assert s.startswith('"') and s.endswith('"')


def test_old_report_without_run_block_says_not_recorded():
    html = render_genome_report(_doc("fola_50a_old_format"), sample="Fola 50a")
    assert "not recorded" in html
    assert "Fola 50a" in html


def test_run_block_is_shown():
    html = render_genome_report(_doc("phycomyces_classifier"))
    assert "Phycomyces blakesleeanus" in html and "taxid 4837" in html
    assert "N50" in html and "3f2a9c41d07be5aa" in html


def test_print_mode_opens_collapsed_sections():
    doc = _doc("fola_50a_old_format")
    assert "<details>" in render_genome_report(doc)
    printed = render_genome_report(doc, print_mode=True)
    assert "<details>" not in printed and "<details open>" in printed


def test_overlapping_genes_and_labels_go_to_separate_lanes():
    extents = [(100, 300), (150, 200), (310, 400), (600, 700)]
    assert assign_lanes(extents, [True, False, True, True]) == [0, 1, 0, 0]


def test_contig_end_side_from_edge_distance():
    assert contig_end({"start": 116422, "end": 128495, "contig_edge_distance": 41212}) == ("right", 169707)
    assert contig_end({"start": 1200, "end": 1980, "contig_edge_distance": 1199}) == ("left", 1)
    assert contig_end({"start": 5, "end": 9}) is None


def test_sample_name_fallback(tmp_path):
    assert sample_name({"run": {"sample": "s1"}}, tmp_path, None) == ("s1", False)
    assert sample_name({"run": {"genome": {"file": "GCA_1.fna.gz"}}}, tmp_path, None) == ("GCA_1", False)
    assert sample_name({}, tmp_path, "given") == ("given", False)
    assert sample_name({}, tmp_path, None) == (tmp_path.name, True)


def test_two_idiomorphs_lead_with_review_verdict_and_do_not_repeat_causes():
    html = render_genome_report(_doc("two_idiomorphs"))
    result = html[html.index('id="result-h"'):html.index("Called loci")]
    assert result.index("needs review") < result.index('class="chip"')
    for cause in ("a duplication or paralog", "a mixed culture or heterokaryon"):
        assert result.count(cause) == 1


def test_locus_without_flanks_is_not_called_complete():
    html = render_genome_report(_doc("basidio_hd_pr"))
    pr = html[html.index('id="locus-1"'):]
    assert "no flanking gene identified" in pr
    assert "B43-like</b> (PR family): complete locus" not in html


def test_cli_writes_report_next_to_the_run(tmp_path):
    run = tmp_path / "run"
    run.mkdir()
    (run / "detection_report.yaml").write_text((FIXTURES / "phycomyces_classifier.yaml").read_text())
    (run / "detected_loci.gff3").write_text("##gff-version 3\n")
    assert main(["report", "genome", "--run", str(run)]) == 0
    html = (run / "report.html").read_text()
    assert "detected_loci.gff3" in html and "Phybl2" in html


def test_cli_rejects_a_file_that_is_not_a_detection_report(tmp_path):
    bad = tmp_path / "x.yaml"
    bad.write_text("a: 1\n")
    assert main(["report", "genome", "--run", str(bad)]) == 1


def test_no_pdf_engine_gives_a_clear_error(monkeypatch):
    import builtins
    real_import = builtins.__import__

    def no_weasy(name, *a, **k):
        if name == "weasyprint":
            raise ImportError(name)
        return real_import(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", no_weasy)
    monkeypatch.setattr(pdf, "find_chrome", lambda: None)
    with pytest.raises(pdf.NoPdfEngine, match="Save as PDF"):
        pdf.available_engine()


def test_pdf_render_with_installed_engine(tmp_path):
    try:
        pdf.available_engine()
    except pdf.NoPdfEngine:
        pytest.skip("no PDF engine installed")
    out = tmp_path / "r.pdf"
    pdf.html_to_pdf(render_genome_report(_doc("phycomyces_classifier"), print_mode=True), out)
    assert out.read_bytes()[:5] == b"%PDF-"


def test_element_ids_are_unique():
    import re
    html = render_genome_report(_doc("two_idiomorphs"))
    ids = re.findall(r'\bid="([^"]+)"', html)
    assert len(ids) == len(set(ids))


def test_figure_colours_do_not_depend_on_page_css():
    """WeasyPrint ignores page CSS inside SVG: role colours must be attributes."""
    html = render_genome_report(_doc("phycomyces_classifier"))
    assert 'fill="#b8471a"' in html and 'fill="#2a6fb0"' in html


def test_single_candidate_margin_names_the_resolution_it_came_from():
    """With one idiomorph candidate the pipeline's margin is the narrowest overlap resolution (identity points);
    the page must say so rather than claim no other idiomorph scored."""
    html = render_genome_report(_doc("fola_50a_old_format"))
    assert "46.7" in html and "identity points" in html
    assert "No other idiomorph scored" not in html
    assert "No other idiomorph scored" in render_genome_report(_doc("basidio_hd_pr"))


def test_candidate_scores_match_the_gloss():
    for case in CASES:
        for call in _doc(case).get("detected") or []:
            if call.get("idiomorph_classifier"):
                continue
            core = [g["bitscore"] for g in call["gene_evidence"] if g.get("role") == "core_MAT" and g.get("bitscore")]
            if core and call["idiomorph_candidates"]:
                assert call["idiomorph_candidates"][0]["score"] == max(core), case
