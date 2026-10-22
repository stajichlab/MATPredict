"""Add GAP1/BGL1 (flanking_conserved) to the two Bipolaris maydis MAT records. Text edits only."""
import sys
from pathlib import Path
root = Path(sys.argv[1])
SPEC = {
 '5016_c4_MAT_MAT1-2': dict(acc='AF027687.1', mat_start=5137, mat_end=6223, end=11742,
    gap=('AAB84002.1', [(2629,3006),(1965,2574),(1153,1917),(794,1073),(633,744)]),
    bgl=('AAB84005.1', [(8898,9003),(9072,9449),(9505,9720),(9830,11742)]),
    old_orient='GAP1->ORF1->MAT1-2-1->BGL1', new_orient='GAP1(-)->ORF1->MAT1-2-1(+)->BGL1(+)',
    old_note="""Curated: MAT1-2-1 (AAB84004.1).
      Not curated: GAP1 (AAB84002.1), unknown ORF (AAB84003.1), BGL1 (AAB84005.1).
      Idiomorph 1,171 bp per the paper.'""",
    new_note="""Curated: MAT1-2-1 (AAB84004.1); GAP1 (AAB84002.1) and BGL1 (AAB84005.1) added
      2026-10-22 as Pleosporaceae-scoped roster flanks (prototype; db/Ascomycota/order.yml),
      CDS joins translated from the deposit equal the protein records. Not curated:
      unknown ORF1 (AAB84003.1). The segment spans GAP1 to BGL1. Idiomorph 1,171 bp
      per the paper.'"""),
 '5016_c5_MAT_MAT1-1': dict(acc='AF029913.1', mat_start=5147, mat_end=6350, end=11866,
    gap=('AAB82943.1', [(2629,3006),(1965,2574),(1153,1917),(794,1073),(633,744)]),
    bgl=('AAB82946.1', [(9022,9127),(9196,9573),(9629,9844),(9954,11866)]),
    old_orient='GAP1->ORF1->MAT1-1-1->BGL1', new_orient='GAP1(-)->ORF1->MAT1-1-1(+)->BGL1(+)',
    old_note="""      MAT1-1-1 (AAB82945.1). Not curated, as they are not roster genes: GAP1 (AAB82943.1),
      an unknown ORF (AAB82944.1) and BGL1 (AAB82946.1), the Pleosporales flanks.
      The deposit has no APN2 or SLA2. Idiomorph 1,297 bp per the paper.'""",
    new_note="""      MAT1-1-1 (AAB82945.1); GAP1 (AAB82943.1) and BGL1 (AAB82946.1) added 2026-10-22
      as Pleosporaceae-scoped roster flanks (prototype; db/Ascomycota/order.yml), CDS
      joins translated from the deposit equal the protein records. Not curated: the
      unknown ORF1 (AAB82944.1). The segment spans GAP1 to BGL1. The deposit has no
      APN2 or SLA2. Idiomorph 1,297 bp per the paper.'"""),
}
def gene_block(idx, name, acc, exons, strand, order):
    lo = min(min(e) for e in exons); hi = max(max(e) for e in exons)
    s = f"""- gene_index: {idx}
  name: {name}
  protein_accession: ncbi_protein:{acc}
  role: flanking_conserved
  present: true
  locus_tag: null
  segment_index: 0
  start: {lo}
  end: {hi}
  strand: '{strand}'
  order_in_locus: {order}
  exons:
"""
    for a, b in exons:
        s += f"  - start: {a}\n    end: {b}\n"
    return s
for rec, sp in SPEC.items():
    p = root / rec / 'metadata.yaml'
    t = p.read_text()
    def sub(old, new):
        global t
        assert t.count(old) == 1, (rec, old)
        t = t.replace(old, new)
    sub(f"reference_orientation: {sp['old_orient']}", f"reference_orientation: {sp['new_orient']}")
    sub(sp['old_note'], sp['new_note'])
    sub(f"      start: {sp['mat_start']}\n      end: {sp['mat_end']}\n      contig_edge_distance: null",
        f"      start: 633\n      end: {sp['end']}\n      contig_edge_distance: null")
    # MAT gene becomes order_in_locus 1 (GAP1 is upstream)
    sub("  order_in_locus: 0\n", "  order_in_locus: 1\n")
    add = gene_block(1, 'GAP1', sp['gap'][0], sp['gap'][1], '-', 0) + gene_block(2, 'BGL1', sp['bgl'][0], sp['bgl'][1], '+', 2)
    sub("evidence:\n  locus_existence:", add + "evidence:\n  locus_existence:")
    pg = "    - gene_index: 0\n      percent_identity: 100.0\n      coverage: 100.0\n      status: pass\n"
    sub(pg, pg + pg.replace("gene_index: 0", "gene_index: 1") + pg.replace("gene_index: 0", "gene_index: 2"))
    p.write_text(t)
    print('edited', p)
