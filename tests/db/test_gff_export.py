from __future__ import annotations

from Bio import SeqIO

from MATPredict.db.gff_export import write_genbank, write_gff3, write_proteins_fasta

RECORD = {
    "record_id": "4837_nrrl-1555_MAT_Plus",
    "locus": {
        "core": {
            "segments": [{"segment_index": 0, "sequence_source": {"seq_region": "scaffold_3"}, "start": 120345, "end": 128900}],
        },
    },
    "genes": [
        {"gene_index": 0, "name": "sexP", "role": "core_MAT", "present": True, "segment_index": 0,
         "start": 121002, "end": 122400, "strand": "+"},
        {"gene_index": 1, "name": "tptA", "role": "flanking_conserved", "present": True, "segment_index": 0,
         "start": 120345, "end": 121000, "strand": "+"},
        {"gene_index": 2, "name": "sexM", "role": "core_MAT", "present": False, "segment_index": 0,
         "start": None, "end": None, "strand": None},
    ],
}


def test_write_gff3_includes_only_present_genes(tmp_path):
    out_path = tmp_path / "locus.gff3"
    write_gff3(RECORD, out_path)
    content = out_path.read_text()
    assert "sexP" in content
    assert "tptA" in content
    assert "sexM" not in content
    assert content.startswith("##gff-version 3")


def test_write_proteins_fasta_header_convention(tmp_path):
    out_path = tmp_path / "proteins.faa"
    write_proteins_fasta(RECORD, sequences={0: "MKTAYIAKQRQ", 1: "GATTACAGATTACA"}, out_path=out_path)
    content = out_path.read_text()
    assert ">4837_nrrl-1555_MAT_Plus|gene_index=0|name=sexP|role=core_MAT" in content
    assert "MKTAYIAKQRQ" in content


def test_write_genbank_includes_only_present_genes_and_roundtrips(tmp_path):
    out_path = tmp_path / "locus.gbk"
    write_genbank(RECORD, sequences={0: "MKTAYIAKQRQ", 1: "GATTACAGATTACA"}, out_path=out_path)

    assert out_path.exists()
    content = out_path.read_text()
    assert "sexP" in content
    assert "tptA" in content
    assert "sexM" not in content

    records = list(SeqIO.parse(out_path, "genbank"))
    assert len(records) == 1
    gene_names = set()
    for rec in records:
        for feature in rec.features:
            if feature.type == "gene":
                gene_names.add(feature.qualifiers["gene"][0])
    assert gene_names == {"sexP", "tptA"}
