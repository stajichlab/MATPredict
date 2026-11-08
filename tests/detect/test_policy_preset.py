"""`detect --preset NAME`: a named policy set that reproduces the v0.7.0 site data's flags in one option.

The code defaults stay as they were (every policy off). A preset fills only the policy options that were NOT given on
the command line, so an explicit flag wins even when it names the default value. The run records the preset it used.
"""
import argparse

import pytest

from MATPredict.detect import cli
from MATPredict.detect.cli import PRESETS, apply_preset, register_subcommands

BASE = ["detect", "--genome", "g.fa", "--taxid", "1"]
POLICIES = ("partial_core_policy", "homothallic_policy", "confidence_cap_policy", "strong_single_gene_policy",
            "scoped_roster_policy")


def _parse(*extra):
    parser = argparse.ArgumentParser()
    register_subcommands(parser.add_subparsers(dest="cmd"))
    return parser.parse_args(BASE + list(extra))


def _resolved(args):
    return {k: getattr(args, k) for k in POLICIES}


def test_the_v07_preset_is_the_clean_run_setting():
    assert PRESETS["v0.7"] == {"partial_core_policy": "coverage_profile", "homothallic_policy": "both_genes",
                               "confidence_cap_policy": "off", "strong_single_gene_policy": "off",
                               "scoped_roster_policy": "on"}


def test_without_a_preset_every_policy_keeps_its_default():
    args = apply_preset(_parse())
    assert args.preset is None
    assert _resolved(args) == {"partial_core_policy": "off", "homothallic_policy": "off", "confidence_cap_policy": "off",
                               "strong_single_gene_policy": "off", "scoped_roster_policy": "off"}


def test_the_preset_sets_the_five_policies():
    args = apply_preset(_parse("--preset", "v0.7"))
    assert _resolved(args) == PRESETS["v0.7"]
    assert args.preset == "v0.7"


def test_an_explicit_flag_beats_the_preset_even_when_it_names_the_default():
    args = apply_preset(_parse("--preset", "v0.7", "--homothallic-policy", "off", "--scoped-roster-policy", "off"))
    assert args.homothallic_policy == "off" and args.scoped_roster_policy == "off"      # kept, though they are the defaults
    assert args.partial_core_policy == "coverage_profile"                              # the rest still come from the preset
    args = apply_preset(_parse("--preset", "v0.7", "--strong-single-gene-policy", "on"))
    assert args.strong_single_gene_policy == "on" and args.homothallic_policy == "both_genes"


def test_flag_order_on_the_command_line_does_not_matter():
    a = apply_preset(_parse("--partial-core-policy", "off", "--preset", "v0.7"))
    b = apply_preset(_parse("--preset", "v0.7", "--partial-core-policy", "off"))
    assert _resolved(a) == _resolved(b) and a.partial_core_policy == "off"


def test_the_equals_and_abbreviated_forms_count_as_explicit():
    args = apply_preset(_parse("--preset=v0.7", "--homothallic-policy=off"))
    assert args.homothallic_policy == "off" and args.partial_core_policy == "coverage_profile"


def test_the_help_text_renders_with_the_preset_values(capsys):
    parser = argparse.ArgumentParser()
    register_subcommands(parser.add_subparsers(dest="cmd"))
    with pytest.raises(SystemExit):
        parser.parse_args(["detect", "--help"])
    out = capsys.readouterr().out
    assert "--preset" in out and "v0.7" in out and "homothallic-policy=both_genes" in out


def test_an_unknown_preset_is_rejected_by_the_parser():
    with pytest.raises(SystemExit):
        _parse("--preset", "v9.9")


def test_the_options_still_validate_their_values():
    with pytest.raises(SystemExit):
        _parse("--homothallic-policy", "nonsense")


def test_the_run_block_records_the_preset_and_the_resolved_policies(tmp_path):
    from types import SimpleNamespace

    from MATPredict.detect.provenance import RunClock

    db = tmp_path / "db"
    db.mkdir()
    for extra, expect in (([], None), (["--preset", "v0.7"], "v0.7")):
        args = apply_preset(_parse(*extra))
        doc = cli._run_doc(args, SimpleNamespace(db_root=db), SimpleNamespace(), RunClock())
        params = doc["parameters"]
        assert params["preset"] == expect
        assert params["partial_core_policy"] == ("coverage_profile" if expect else "off")
        assert params["scoped_roster_policy"] == ("on" if expect else "off")


def _flag_of(dest):
    """The option string of the `detect` option whose destination is `dest` (a flag need not be its dest spelled with dashes)."""
    parser = argparse.ArgumentParser()
    sub = register_subcommands(parser.add_subparsers(dest="cmd"))
    detect = parser._subparsers._group_actions[0].choices["detect"]
    return next(a.option_strings[0] for a in detect._actions if a.dest == dest and a.option_strings)


def test_every_preset_value_is_a_valid_choice_for_its_option():
    for name, values in PRESETS.items():
        extra = []
        for dest, value in values.items():
            extra += [_flag_of(dest), value]
        _parse(*extra)       # the parser would exit on an invalid choice


def test_the_cmd_detect_handler_applies_the_preset_first(monkeypatch):
    def stop(genome):
        raise RuntimeError("stop right after the preset was applied")

    monkeypatch.setattr(cli, "_refuse_compressed_genome", stop)
    args = _parse("--preset", "v0.7", "--out-dir", "o")
    assert _resolved(args)["homothallic_policy"] == "off"        # not applied yet: parsing alone changes nothing
    with pytest.raises(RuntimeError):
        cli._cmd_detect(args)
    assert _resolved(args) == PRESETS["v0.7"]                    # applied in place by the handler


def test_the_v09_preset_is_v081_plus_the_three_report_only_rules():
    v081, v09 = PRESETS["v0.8.1"], PRESETS["v0.9"]
    added = {k: v for k, v in v09.items() if v081.get(k) != v}
    assert added == {"cross_family_policy": "keep_labelled", "idiomorph_margin_policy": "separate", "out_of_scope_label_policy": "on"}
    assert all(v09[k] == v for k, v in v081.items() if k != "cross_family_policy")


def test_only_the_approved_tier_and_label_rules_are_in_a_preset():
    held_out = ("weak_fragment_cap_policy", "weak_coverage_noise_policy", "strong_single_gene_policy", "confidence_cap_policy")
    approved = {"weak_flank_cap_policy": "on", "short_fragment_vote_policy": "undetermined", "single_region_cap_policy": "on"}
    for name, vals in PRESETS.items():
        assert all(vals.get(k, "off") == "off" for k in held_out), name
        got = {k: vals[k] for k in approved if k in vals}
        assert got == (approved if name == "v0.9.1" else {}), name


def test_the_v091_preset_is_v09_plus_the_three_approved_rules():
    v09, v091 = PRESETS["v0.9"], PRESETS["v0.9.1"]
    assert {k: v for k, v in v091.items() if v09.get(k) != v} == {"weak_flank_cap_policy": "on", "short_fragment_vote_policy": "undetermined", "single_region_cap_policy": "on"}
    assert all(v091[k] == v for k, v in v09.items())
