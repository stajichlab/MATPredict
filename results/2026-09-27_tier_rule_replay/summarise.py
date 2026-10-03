import csv, glob, collections
out=open('summary.txt','w')
def p(*a):
    s=' '.join(str(x) for x in a); print(s); out.write(s+'\n')
allrows=[]
for f in sorted(glob.glob('calls_*.tsv')):
    rows=[r for r in csv.DictReader(open(f),delimiter='\t') if r.get('reported')]
    allrows+=rows
    ok=sum(r['reproduced_ok']=='True' for r in rows)
    p(f"\n== {f[6:-4]}: calls {len(rows)}, reproduced {ok}/{len(rows)} ({100*ok/max(1,len(rows)):.1f}%)")
    mm=collections.Counter((r['reported'],r['reproduced'],r['why_current']) for r in rows if r['reproduced_ok']!='True')
    for k,v in mm.most_common(6): p("   mismatch", k, v)
    fam=collections.defaultdict(collections.Counter)
    for r in rows:
        if r['reproduced_ok']!='True': continue
        c=fam[r['family']]
        c['medium']+= r['reported']=='medium'
        c['A_rise']+= r['reported']=='medium' and r['variant_A']=='high'
        c['B_rise']+= r['reported']=='medium' and r['variant_B']=='high'
        c['Bp_rise']+= r['reported']=='medium' and r['variant_Bp']=='high'
        c['A_other_change']+= r['variant_A']!=r['reported'] and not (r['reported']=='medium' and r['variant_A']=='high')
        c['B_other_change']+= r['variant_B']!=r['reported'] and not (r['reported']=='medium' and r['variant_B']=='high')
    for k,c in sorted(fam.items()):
        if c['medium'] or c['A_rise'] or c['B_rise']:
            p(f"   {k:28s} medium {c['medium']:4d}  A->high {c['A_rise']:4d}  B->high {c['B_rise']:4d}  B'->high {c['Bp_rise']:4d}  other A/B {c['A_other_change']}/{c['B_other_change']}")
risers=[r for r in allrows if r['reproduced_ok']=='True' and r['reported']=='medium' and (r['variant_A']=='high' or r['variant_B']=='high' or r['variant_Bp']=='high')]
with open('risers.tsv','w',newline='') as fo:
    w=csv.DictWriter(fo,fieldnames=list(risers[0]) if risers else ['run'],delimiter='\t'); w.writeheader(); w.writerows(risers)
p(f"\nrisers total A {sum(r['variant_A']=='high' for r in risers)}  B {sum(r['variant_B']=='high' for r in risers)}  B' {sum(r['variant_Bp']=='high' for r in risers)}")
strong=[r for r in risers if r['max_absent_identity'] and float(r['max_absent_identity'])>60]
p("risers with an ignored gene >60% identity:", len(strong))
for r in strong[:0]: p("   ", r['run'], r['genome'], r['family'], r['idiomorph'], r['absent_genes'], '| core', r['core_evidence'])

bystrong=collections.Counter((r['run'],r['family'],r['variant_A']=='high',r['variant_Bp']=='high') for r in strong)
p("strong-ignored risers by run/family (A_rises, Bp_rises):")
for k,v in sorted(bystrong.items()): p("   ",k,v)
