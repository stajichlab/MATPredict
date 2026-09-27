"""Assemble the final gene set: GFF-derived genes that translate cleanly on the
public contig, plus miniprot re-models for the genes whose GFF model did not."""
import json
from Bio import SeqIO
from Bio.Seq import Seq

N = "ncbi"
g = json.load(open("genes.json"))
keep = {  # (strain, gene) -> use GFF model
    ("CBS14", "STE3_A1"), ("CBS14", "HD1"),
    ("NBRC_0880", "STE3_A2"), ("NBRC_0880", "STE20"), ("NBRC_0880", "HD1"), ("NBRC_0880", "HD2"),
    ("CBS_20", "STE3_A1"), ("CBS_20", "HD1"), ("CBS_20", "HD2"),
    ("JY1105", "STE3_A1"), ("JY1105", "STE20"), ("JY1105", "HD1"),
    ("Y-7192", "STE3_A2"), ("Y-7192", "STE20"),
    ("LS11", "HD1"), ("LS11", "HD2"),
    ("JJ10.1", "STE3_A2"), ("JJ10.1", "STE20"),
}
out = [dict(o, model_source="group_gff") for o in g if (o["strain"], o["gene"]) in keep]
for o in out:
    assert o["internal_stops"] == 0, o["strain"] + o["gene"]

REMODEL = [  # name, strain, gene, contig, window_start, mRNA ID, query used
    ("CBS14_HD2", "CBS14", "HD2", "CAKLCE030000016.1", 68000, "MP000001", "R. toruloides NBRC 0880 HD2 (this set)"),
    ("JY1105_HD2", "JY1105", "HD2", "JANBVD010000009.1", 483000, "MP000001", "R. mucilaginosa HD2 KAG0665576"),
    ("CBS20_STE20", "CBS_20", "STE20", "CAKKSX020000013.1", 52000, "MP000001", "Rhodotorula STE20 panel (ctsai085 MAT_ref/PR/ste20.faa)"),
]
for name, strain, gene, contig, w0, mid, query in REMODEL:
    cds, strand = [], None
    for line in open(f"mp_{name}.gff"):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS" or f"Parent={mid}" not in f[8]:
            continue
        cds.append((int(f[3]) + w0 - 1, int(f[4]) + w0 - 1)); strand = f[6]
    cds.sort()
    seq = str(next(SeqIO.parse(f"{N}/{contig}.fa", "fasta")).seq).upper()
    nt = "".join(seq[a - 1:b] for a, b in cds)
    if strand == "-":
        nt = str(Seq(nt).reverse_complement())
    # include the stop codon if present right after the last CDS
    prot = str(Seq(nt[: len(nt) - len(nt) % 3]).translate())
    out.append(dict(strain=strain, locus="HD" if gene.startswith("HD") else "PR", gene=gene,
                    id=f"miniprot:{name}", contig=contig, start=cds[0][0], end=cds[-1][1], strand=strand,
                    exons=[list(c) for c in cds], codon_start=1, aa_len=len(prot.rstrip("*")),
                    starts_M=prot.startswith("M"), internal_stops=prot.rstrip("*").count("*"),
                    ends_stop=prot.endswith("*"), protein=prot.rstrip("*"),
                    model_source=f"miniprot 0.x from {query}"))
json.dump(out, open("final_genes.json", "w"), indent=1)
for o in sorted(out, key=lambda o: (o["strain"], o["locus"], o["gene"])):
    print(o["strain"], o["locus"], o["gene"], o["contig"], o["start"], o["end"], o["strand"], len(o["exons"]),
          "aa", o["aa_len"], "M", o["starts_M"], "stops", o["internal_stops"], o["model_source"][:30])
