"""The reads screen rule: flags from a GmG record, the trigger for a deeper screen, and the command line."""
import importlib.util
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("screen_reads", REPO / "datasets/reads/sra_survey/pilot/screen_reads.py")
sr = importlib.util.module_from_spec(spec)
sys.modules["screen_reads"] = sr
spec.loader.exec_module(sr)

COCC = ["Ascomycota", "Eurotiomycetes", "Onygenales", "Onygenaceae", "Coccidioides", "Coccidioides posadasii"]
HISTO = ["Ascomycota", "Eurotiomycetes", "Onygenales", "Ajellomycetaceae", "Histoplasma", "Histoplasma capsulatum"]
LECAN = ["Ascomycota", "Sordariomycetes", "Hypocreales", "Cordycipitaceae", "Lecanicillium", "Lecanicillium sp."]


def rec(rid, genus, lineages, second=None, hashes=0, conc=0.0, suspected=False, status="ok", cov=12.0):
    return {"id": rid, "status": status, "qc_flags": ["mixture_suspected"] if suspected else [],
            "guess": {"calls": {"genus": {"taxon": genus}, "species": {"taxon": None, "vote": (lineages[0][-1] if lineages else None)}} if genus else {},
                      "hits": [{"lineage": lin} for lin in lineages]},
            "qc": {"coverage_est": cov, "mixture": {"primary": genus, "suspected": suspected, "explained_fraction": 0.9,
                   "second": {"taxon": second, "hashes": hashes, "explained_fraction": hashes / 2500, "concentration": conc, "ref_coverage": 0.1} if second else None}}}


CLEAN = rec("r", "Coccidioides", [COCC] * 3, second="Aspergillus", hashes=40, conc=0.17)      # a chance second genus: many hashes, scattered
WHOLE = rec("r", "Histoplasma", [HISTO] * 5)                                                   # NM_7898: the whole set is another organism
MIXED = rec("r", "Lecanicillium", [LECAN, LECAN, COCC], second="Coccidioides", hashes=2334, conc=0.94, suspected=True)   # TX7
BORDER = rec("r", "Coccidioides", [COCC] * 3, second="Gigaspora", hashes=13, conc=0.85)


def test_expected_genus_comes_from_the_list_label():
    assert sr.expected_genus("aspergillus_fumigatus") == "Aspergillus"
    assert sr.expected_genus("Candidozyma auris") == "Candidozyma"
    assert sr.expected_genus("") == ""


def test_a_clean_set_with_a_chance_second_genus_is_ok_and_does_not_trigger():
    assert sr.trigger(CLEAN) == "ok"
    out = sr.finalize(CLEAN, "coccidioides_posadasii", 2_000_000)
    assert out["status"] == "ok" and out["flags"] == [] and out["rescreened"] is False


def test_a_read_set_that_is_another_genus_is_a_genus_mismatch():
    out = sr.finalize(WHOLE, "coccidioides_posadasii", 2_000_000)
    assert out["status"] == "genus_mismatch" and out["flags"] == ["genus_mismatch"] and out["genus"] == "Histoplasma"


def test_a_mixture_is_mixed_and_not_a_genus_mismatch_because_the_expected_genus_is_among_the_hits():
    out = sr.finalize(MIXED, "coccidioides_posadasii", 2_000_000)
    assert out["status"] == "mixed" and out["flags"] == ["mixed"]
    assert (out["second_genus"], out["second_hashes"]) == ("Coccidioides", 2334)


def test_a_mixture_whose_expected_genus_is_not_among_the_hits_is_both():
    out = sr.finalize(rec("r", "Lecanicillium", [LECAN] * 3, second="Histoplasma", hashes=900, conc=0.9, suspected=True), "coccidioides_posadasii", 2_000_000)
    assert out["flags"] == ["genus_mismatch", "mixed"] and out["status"] == "genus_mismatch"      # the whole-set problem comes first


def test_a_scattered_but_named_second_genus_triggers_only_when_it_is_concentrated_enough():
    assert sr.trigger(BORDER) == "deeper"
    assert sr.trigger(rec("r", "Coccidioides", [COCC], second="Gigaspora", hashes=4, conc=0.9)) == "ok"      # too few hashes
    assert sr.trigger(rec("r", "Coccidioides", [COCC], second="Gigaspora", hashes=30, conc=0.69)) == "ok"    # not concentrated
    assert sr.trigger(rec("r", "Coccidioides", [COCC], second="Gigaspora", hashes=5, conc=0.7)) == "deeper"  # the limits themselves
    assert sr.trigger(MIXED) == "ok"                                                                          # already flagged: nothing to ask
    assert sr.trigger(None) == "ok"


def test_a_borderline_set_with_no_deeper_screen_stays_borderline_and_says_why():
    out = sr.finalize(BORDER, "coccidioides_posadasii", 2_000_000, deep_possible=False)
    assert out["status"] == "borderline" and out["deeper_unavailable"] is True


def test_the_deeper_screen_decides():
    flagged = rec("r", "Coccidioides", [COCC] * 3, second="Gigaspora", hashes=400, conc=0.9, suspected=True)
    out = sr.finalize(BORDER, "coccidioides_posadasii", 2_000_000, deep_rec=flagged, deep_spots=6_000_000)
    assert out["status"] == "mixed" and out["rescreened"] and out["rescreen_spots"] == 6_000_000
    assert out["second_hashes"] == 400 and out["at_screened_spots"]["second_hashes"] == 13      # both depths are kept
    resolved = rec("r", "Coccidioides", [COCC] * 3, second="Aspergillus", hashes=60, conc=0.15)
    out = sr.finalize(BORDER, "coccidioides_posadasii", 2_000_000, deep_rec=resolved, deep_spots=6_000_000)
    assert out["status"] == "ok" and out["flags"] == []                                            # chance, not confirmed deeper
    still = rec("r", "Coccidioides", [COCC] * 3, second="Gigaspora", hashes=30, conc=0.8)
    out = sr.finalize(BORDER, "coccidioides_posadasii", 2_000_000, deep_rec=still, deep_spots=6_000_000)
    assert out["status"] == "borderline"                                                           # still looks real, still not flagged


def test_no_genus_call_and_a_missing_record():
    out = sr.finalize(rec("r", None, [], status="no_match"), "coccidioides_posadasii", 2_000_000)
    assert "no_call" in out["flags"]
    assert sr.finalize(None, "coccidioides_posadasii", 2_000_000)["status"] == "no_result"


def test_a_species_call_without_a_genus_call_is_not_a_no_call():
    """Clonostachys rosea reads get a species call and no genus call from GmG: that is a classification (wave 1, 16 runs)."""
    r = rec("r", "Clonostachys", [["Fungi", "Clonostachys", "Clonostachys rosea"]] * 3)
    r["guess"]["calls"] = {"species": {"taxon": "Clonostachys rosea"}}
    out = sr.finalize(r, "clonostachys_rosea", 2_000_000)
    assert out["status"] == "ok" and "no_call" not in out["flags"]
    r["guess"]["calls"] = {}
    assert "no_call" in sr.finalize(r, "clonostachys_rosea", 2_000_000)["flags"]               # neither call: still flagged


def test_without_an_expected_species_the_genus_is_not_checked():
    assert sr.finalize(WHOLE, "", 2_000_000)["status"] == "ok"


def test_the_command_line_trigger_and_finalize(tmp_path, capsys):
    first, deep = tmp_path / "a.jsonl", tmp_path / "b.jsonl"
    first.write_text(json.dumps({**BORDER, "id": "SRR1"}) + "\n" + json.dumps({**CLEAN, "id": "SRR2"}) + "\n")
    deep.write_text(json.dumps({**rec("r", "Coccidioides", [COCC], second="Gigaspora", hashes=400, conc=0.9, suspected=True), "id": "SRR1"}) + "\n")
    assert sr.main(["trigger", "--results", str(first), "--run", "SRR1"]) == 0 and capsys.readouterr().out.strip() == "deeper"
    assert sr.main(["trigger", "--results", str(first), "--run", "SRR2"]) == 0 and capsys.readouterr().out.strip() == "ok"
    out = tmp_path / "SRR1.screen.json"
    sr.main(["finalize", "--results", str(first), "--run", "SRR1", "--group", "coccidioides_posadasii", "--spots", "2000000",
             "--deep-results", str(deep), "--deep-spots", "6000000", "--out", str(out)])
    assert json.loads(out.read_text())["status"] == "mixed"
