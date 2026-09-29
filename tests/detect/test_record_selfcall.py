"""A curated record's source genome must be called at the record (review F2)."""
from MATPredict.detect.record_selfcall import RecordLocation, record_location, selfcall_verdict

LOC = RecordLocation("r1", "MAT", ("Plus",), "GCA_000000001.1", (("CTG1.1", 1000, 5000),))


def test_a_detected_locus_at_the_record_is_called():
    report = {"detected": [{"family": "Mucoromycota:MAT", "contig": "CTG1.1", "start": 900,
                            "end": 2000, "idiomorph": "Plus", "confidence": "high"}]}
    v = selfcall_verdict(LOC, report)
    assert v["verdict"] == "called" and v["idiomorph"] == "Plus"


def test_a_version_difference_in_the_contig_still_matches():
    report = {"detected": [{"family": "X:MAT", "contig": "CTG1.2", "start": 4000, "end": 6000}]}
    assert selfcall_verdict(LOC, report)["verdict"] == "called"


def test_a_suppressed_locus_at_the_record_is_withheld_with_its_reason():
    report = {"detected": [], "suppressed_loci": [
        {"family": "X:MAT", "contig": "CTG1.1", "start": 1000, "end": 5000,
         "withheld_reason": "below_fraction_floor"}]}
    v = selfcall_verdict(LOC, report)
    assert v["verdict"] == "withheld" and v["withheld_reason"] == "below_fraction_floor"


def test_a_call_elsewhere_or_of_another_family_is_missed():
    report = {"detected": [{"family": "X:MAT", "contig": "OTHER.1", "start": 1, "end": 9},
                           {"family": "X:HD", "contig": "CTG1.1", "start": 1000, "end": 5000}]}
    assert selfcall_verdict(LOC, report)["verdict"] == "missed"


def test_record_location_reads_segments_and_assembly(tmp_path):
    p = tmp_path / "metadata.yaml"
    p.write_text(
        "record_id: r1\nmating_type: {locus_name: MAT, idiomorphs: [Plus]}\n"
        "locus:\n  core:\n    definition_note: 'from GCA_002105135.1'\n    segments:\n"
        "    - sequence_source: {type: insdc_nucleotide, accession: MCGN01000004.1}\n"
        "      start: 1755515\n      end: 1760805\n")
    loc = record_location(p)
    assert loc.assembly == "GCA_002105135.1"
    assert loc.segments == (("MCGN01000004.1", 1755515, 1760805),)


def test_the_assembly_field_wins_over_text(tmp_path):
    """`locus.assembly_accession` (2026-09-28) is the record's own statement of
    its source assembly; the text search is only a fallback."""
    p = tmp_path / "metadata.yaml"
    p.write_text(
        "record_id: r1\nmating_type: {locus_name: MAT, idiomorphs: [Plus]}\n"
        "curation: {notes: 'compare GCA_002105135.1'}\n"
        "locus:\n  assembly_accession: GCA_000001405.1\n"
        "  core:\n    definition_note: 'see also GCA_002105135.1'\n    segments: []\n")
    assert record_location(p).assembly == "GCA_000001405.1"


def test_the_schema_accepts_the_assembly_field_and_rejects_a_bad_one():
    from MATPredict.db.schema import load_metadata_schema
    import jsonschema
    prop = load_metadata_schema()["properties"]["locus"]["properties"]["assembly_accession"]
    jsonschema.validate("GCF_025399195.1", prop)
    jsonschema.validate(None, prop)
    try:
        jsonschema.validate("AF029913.1", prop)
    except jsonschema.ValidationError:
        pass
    else:
        raise AssertionError("an INSDC nucleotide accession is not an assembly accession")
