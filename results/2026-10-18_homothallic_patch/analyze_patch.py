"""Outcome tables for the patch. python analyze_patch.py PREP_DIR AUDIT_DIR TRUTH_SET HOLDOUT_MANIFEST FLIP_TSV DIAG_TSV OUT_DIR"""
import sys
from pathlib import Path
import pandas as pd
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 200); pd.set_option('display.max_columns', 40)
P, A, TRUTH, HOLD, FLIP, DIAG, OUT = [Path(x) for x in sys.argv[1:8]]
OUT.mkdir(exist_ok=True)
loci = pd.read_csv(P / 'loci_patched.tsv', sep='\t', low_memory=False)
sets = pd.read_csv(P / 'sets.tsv', sep='\t', dtype=str)
G = pd.read_csv(A / 'genome_audit_final.tsv', sep='\t', low_memory=False)
R = loci[loci.reclassed]
rg = set(R.genome)
print('re-classed loci', len(R), 'genomes', len(rg))
print(R.merge(sets.drop_duplicates('genome')[['genome', 'set']], on='genome').set.value_counts().to_dict())
print('previous class:', R.class_before.value_counts().to_dict(), ' confidence:', R.confidence.value_counts().to_dict(), ' idiomorph:', R.idiomorph.value_counts().to_dict())
print('hmg lost_to:', R.hmg_lost_to.fillna('none').value_counts().to_dict())
print('separation bp: median', R.separation_bp.median(), 'max', R.separation_bp.max(), '; alpha_bits min/med', R.alpha_bits.min(), R.alpha_bits.median(), '; hmg_bits min/med', R.hmg_bits.min(), R.hmg_bits.median())
D = G[['genome', 'genus', 'species', 'order', 'verdict', 'verdict_class', 'second_gene_scan', 'idiomorphs_called', 'primary_idiomorph', 'gw_class']].copy()
D['reclassed'] = D.genome.isin(rg)
print('\nby audit verdict:'); print(pd.crosstab(D.verdict_class, D.reclassed))
conf = D[D.verdict_class == 'conflicting']
print('\nconflicting 298 -> genomes with a candidate:', int(conf.reclassed.sum()))
print(pd.crosstab(conf.genus.where(conf.genus.isin(['Aureobasidium', 'Pyrenophora', 'Didymella', 'Friedmanniomyces', 'Elasticomyces', 'Phyllosticta', 'Stemphylium']), 'other'), conf.reclassed))
dg = pd.read_csv(DIAG, sep='\t')
dg['reclassed'] = dg.genome.isin(rg)
def cause(r):
    f = r.fails
    if r.pass_ == 'relaxed':
        return 'relaxed_pass'
    if 'gene_not_live' in f:
        return 'gene_not_live'
    if 'cross_match' in f:
        return 'cross_match' + ('+cov<50' if 'coverage_lt_50' in f else '')
    return f
dg['cause'] = dg.apply(cause, axis=1)
print('\nby legacy cause:'); print(pd.crosstab(dg.cause, dg.reclassed))
base = ['Alternaria', 'Parastagonospora', 'Venturia', 'Zymoseptoria', 'Ascochyta', 'Corynespora', 'Hortaea', 'Cercospora', 'Fulvia', 'Dothistroma', 'Macrophomina']
b = D[D.genus.isin(base) | (D.species == 'Pyrenophora teres')]
print('\nbaseline heterothallic groups: genomes', len(b), 're-classed', int(b.reclassed.sum()))
print(b[b.reclassed][['genome', 'species']].to_string())
for sp in ['Botryosphaeria dothidea', 'Macrophomina phaseolina']:
    x = D[D.species == sp]; print(sp, len(x), 're-classed', int(x.reclassed.sum()))
two = D[D.idiomorphs_called.astype(str).str.contains(r'\+')]
tl = loci[loci.genome.isin(two.genome)]
print('\ntwo-locus genomes', len(two), 'with a candidate locus', int(two.reclassed.sum()), '; their loci', len(tl), 're-classed', int(tl.reclassed.sum()))
print(D[D.reclassed].groupby(['genus', 'verdict']).size().sort_values(ascending=False).head(30))
# unverifiable / no-call
print('unverifiable re-classed', int(D[(D.verdict_class == 'unverifiable')].reclassed.sum()), 'no_call', int(D[D.verdict_class == 'no_call'].reclassed.sum()))
# idiomorph / confidence unchanged (by construction); record counts
R.to_csv(OUT / 'reclassed_loci.tsv', sep='\t', index=False)
D[D.reclassed].to_csv(OUT / 'reclassed_genomes_with_audit_verdict.tsv', sep='\t', index=False)
# truth, flip, holdout, extra
ts = pd.read_csv(TRUTH, sep='\t')
tg = sets[sets.set == 'truth'].genome
t = ts[ts.genome.isin(tg)].copy()
t['reclassed'] = t.genome.isin(rg)
print('\ntruth set', len(t), 're-classed', int(t.reclassed.sum())); print(t[t.reclassed][['genome', 'organism', 'stated']].to_string())
print(t.stated.value_counts().to_dict())
for s in ('flip', 'holdout', 'extra_afum', 'extra_neurospora_sordaria'):
    g = set(sets[sets.set == s].genome)
    print(s, len(g), 're-classed', len(g & rg), sorted(g & rg)[:6])
fl = pd.read_csv(FLIP, sep='\t')[['genome', 'organism', 'verdict']]
fl['reclassed'] = fl.genome.isin(rg)
print(fl[fl.verdict.str.contains('BOTH|both|fused', regex=True)].to_string())
