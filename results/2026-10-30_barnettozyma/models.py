"""Hand-curated gene models of the B. californica PB4207 MAT locus (JAMAIE010000005.1), from sequence alone (no RNA-seq exists for the genus).
Models are defined in window coordinates (data/GCA_023627795.1_JAMAIE010000005.1_589225-655853.fa; genome = window + 589224) and checked here:
start codon, stop codon, canonical GT..AG introns, translation == the ORF used for scoring. Writes models.gff3, models_proteins.faa, models.tsv.
Evidence behind each model: README.md."""
import sys
sys.path.insert(0, ".")
from gm import load, sub, protein, cds, introns
FA = "data/GCA_023627795.1_JAMAIE010000005.1_589225-655853.fa"
CONTIG, OFF = "JAMAIE010000005.1", 589224
S = load(FA)

def first_M(seg_a, seg_b, strand):
    """Single-exon model from a stop-to-stop segment: first ATG to the stop codon (stop included in the CDS)."""
    p = protein(S, [(seg_a, seg_b)], strand); k = p.find("M")
    if strand == "+":
        a = seg_a + 3 * k; b = seg_b + 3          # stop codon follows the segment
    else:
        a = seg_a - 3; b = seg_b - 3 * k          # stop codon precedes (lower coordinates)
    return [(a, b)]

MODELS = [
    # name, strand, exons (window coords), role, note
    ("dic1", "-", first_M(16017, 17015, "-"), "flanking_conserved", "mitochondrial dicarboxylate carrier; single exon; phmmer 362 bits to DIC1 references"),
    ("MATa2", "-", [(24944, 25231), (24366, 24887)], "core_MAT",
     "HMG-box-like MATa2 (no Pfam HMG_box hit, as for Komagataella); ONE-INTRON MODEL chosen by exhaustive search (82 bits to the Saccharomycotina MATa2 profile vs 64 and 50 for the next donors); 56 bp intron GTATGT..CTTAG"),
    ("MATa1", "+", first_M(25832, 26341, "+"), "core_MAT", "homeodomain (Pfam 53 bits; helix 3 ...WFTNKR...), single exon, divergent partner of MATa2"),
    ("MATalpha1", "+", first_M(29172, 29876, "+"), "core_MAT", "alpha-box (MATalpha_HMGbox 33 bits, 62 bits to curated MATalpha1); zero-intron ORF; the C-terminus is uncertain (alternative one-intron models score the same)"),
    ("sla2", "-", first_M(38510, 41662, "-"), "flanking_conserved", "single exon, ANTH + I_LWEQ, 1217 bits to curated SLA2"),
]
rows, gff, faa = [], ["##gff-version 3"], []
for name, strand, exons, role, note in MODELS:
    ex = sorted(exons)
    p = protein(S, exons, strand)
    ok_start = p.startswith("M"); ok_stop = p.endswith("*") and p.count("*") == 1
    spl = []
    for a, b in introns(exons, strand):
        t = sub(S, a, b, strand); spl.append(f"{a + OFF}-{b + OFF} {len(t)}bp {t[:2]}..{t[-2:]}")
    ga, gb = ex[0][0] + OFF, ex[-1][1] + OFF
    rows.append((name, role, strand, ga, gb, len(exons), len(p) - 1, "yes" if ok_start else "NO", "yes" if ok_stop else "NO", "; ".join(spl) or "-", note))
    gid = f"BC_{name}"
    gff.append(f"{CONTIG}\tcuration\tgene\t{ga}\t{gb}\t.\t{strand}\t.\tID={gid};Name={name}")
    gff.append(f"{CONTIG}\tcuration\tmRNA\t{ga}\t{gb}\t.\t{strand}\t.\tID={gid}.t1;Parent={gid}")
    order = ex if strand == "+" else ex[::-1]
    phase = 0
    for a, b in order:
        gff.append(f"{CONTIG}\tcuration\tCDS\t{a + OFF}\t{b + OFF}\t.\t{strand}\t{phase}\tID={gid}.cds;Parent={gid}.t1")
        phase = (3 - ((b - a + 1) - phase) % 3) % 3
    faa.append(f">{gid} {CONTIG}:{ga}-{gb}({strand}) {len(p)-1} aa\n{p.rstrip('*')}\n")
open("models.gff3", "w").write("\n".join(gff) + "\n")
open("models_proteins.faa", "w").write("".join(faa))
with open("models.tsv", "w") as f:
    f.write("gene\trole\tstrand\tstart\tend\texons\taa\tM_start\tstop\tintrons\tnote\n")
    for r in rows: f.write("\t".join(map(str, r)) + "\n")
for r in rows: print(*r[:10], sep="\t")
