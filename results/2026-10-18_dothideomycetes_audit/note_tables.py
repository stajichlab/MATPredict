"""Markdown fragments for the analysis note -> tables/note_*.md.  python note_tables.py DIR"""
import sys
from pathlib import Path
import pandas as pd
D = Path(sys.argv[1]); T = D / 'tables'
G = pd.read_csv(D / 'genome_audit_final.tsv', sep='\t', low_memory=False)
G['order'] = G['order'].fillna('(blank)')


def md(df):
    cols = list(df.columns)
    out = ['| ' + ' | '.join(cols) + ' |', '|' + '|'.join(['---'] * len(cols)) + '|']
    for _, r in df.iterrows():
        out.append('| ' + ' | '.join('' if pd.isna(v) else (str(int(v)) if isinstance(v, float) and v == int(v) else str(v)) for v in r.values) + ' |')
    return '\n'.join(out)


def agg(df, by):
    g = df.groupby(by)
    t = pd.DataFrame({
        'genomes': g.size(),
        'report': g.apply(lambda x: (x.status != 'no_report').sum()),
        'called': g.apply(lambda x: x.primary_idiomorph.notna().sum()),
        'MAT1-1': g.apply(lambda x: (x.primary_idiomorph == 'MAT1-1').sum()),
        'MAT1-2': g.apply(lambda x: (x.primary_idiomorph == 'MAT1-2').sum()),
        'two-idiomorph loci': g.apply(lambda x: x.idiomorphs_called.astype(str).str.contains(r'\+').sum()),
        'not detected': g.apply(lambda x: (x.status == 'not_detected').sum()),
        'v0.6.0 called': g.apply(lambda x: (x.v060_idiomorphs.astype(str).str.len() > 0).sum() - (x.v060_idiomorphs.astype(str) == 'nan').sum()),
        'confirmed': g.apply(lambda x: (x.verdict_class == 'confirmed').sum()),
        'conflicting': g.apply(lambda x: (x.verdict_class == 'conflicting').sum()),
        'unverifiable': g.apply(lambda x: (x.verdict_class == 'unverifiable').sum()),
    })
    return t


o = agg(G, 'order').sort_values('genomes', ascending=False).reset_index()
tot = pd.DataFrame([dict(order='**all**', **{c: int(o[c].sum()) for c in o.columns if c != 'order'})])
(T / 'note_orders.md').write_text(md(pd.concat([o, tot])))
g = agg(G, ['order', 'genus']).reset_index()
lab = pd.read_csv(T / 'group_genus.tsv', sep='\t')[['order', 'genus', 'label']]
g = g.merge(lab, on=['order', 'genus'], how='left').sort_values('genomes', ascending=False)
(T / 'note_genera_ge10.md').write_text(md(g[g.genomes >= 10]))
(T / 'note_genera_all.tsv').write_text(g.to_csv(sep='\t', index=False))
s = pd.read_csv(T / 'suspect_groups.tsv', sep='\t')
cols = ['group', 'genomes', 'called', 'MAT11', 'MAT12', 'vote_profile_agree', 'vote_profile_disagree', 'profile_undetermined', 'scanned',
        'scan_second_gene_strong', 'scan_second_gene_weak', 'scan_none', 'core_set_aside_for_optional', 'confirmed', 'conflicting', 'unverifiable', 'no_call']
(T / 'note_suspects.md').write_text(md(s[cols]))
print('ok', len(o), len(g))
