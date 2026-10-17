import glob,csv,collections,os
meta={r['genome']:r for r in csv.DictReader(open('dipo_meta.tsv'),delimiter='\t')}
loci=collections.defaultdict(list)
for r in csv.DictReader(open('loci.tsv'),delimiter='\t'):
    loci[r['genome']].append((r['contig'],int(r['start']),int(r['end']),r['family_called'],r['locus_class']))
core={'Ylip_AJ617306|mata1','Ylip_AJ617306|mata2','Ylip_AJ617307|matb1','Ylip_AJ617307|matb2','Gcand_HF558448|MATA','Gcand_HF558449|MATB','Ohpoly_AJ617305|matalpha1','Ohpoly_AJ617305|matalpha2','Ohpoly_AJ617305|mata1'}
rows=[]
for g,m in meta.items():
    f=f'all/out_{g}.tsv'
    if not os.path.exists(f): continue
    H=[]
    for l in open(f):
        q,s,pid,ln,ql,bs,a,b,ev=l.rstrip().split('\t')
        if float(ev)>1e-5: continue
        H.append(('|'.join(q.split('|')[:2]),s,float(pid),int(ln)/int(ql),float(bs),min(int(a),int(b)),max(int(a),int(b))))
    sla=[h for h in H if h[0] in('Ylip_AJ617306|sla2','Ohpoly_AJ617305|sla2') and h[3]>=0.5]
    sla.sort(key=lambda h:-h[4]); S=sla[0] if sla else None
    cn=[]
    if S:
        cn=[h for h in H if h[0] in core and h[1]==S[1] and abs(h[5]-S[5])<15000 and h[3]>=0.3 and not (h[5]>=S[5] and h[6]<=S[6])]
    typ='no_SLA2' if not S else 'SLA2_no_core_nearby'
    if cn:
        names={h[0].split('|')[1] for h in cn}
        best=max(cn,key=lambda h:h[4])
        if any(h[0]=='Ylip_AJ617307|matb2' and h[2]>=60 for h in cn) or any(h[0]=='Ylip_AJ617306|mata2' and h[2]>=60 for h in cn): typ='Yarrowia_clade_MAT'
        elif any(h[0]=='Gcand_HF558448|MATA' and h[2]>=40 for h in cn): typ='HMG_MATA_like'
        else: typ='alpha1_like_at_SLA2'
    # concordance for called
    conc=''
    if m['status']=='called' and S:
        conc='at_SLA2' if any(c==S[1] and min(e+30000,S[6])-max(s-30000,S[5])>0 for c,s,e,_,_ in loci[g]) else 'elsewhere'
    gene_only=sum(1 for x in loci[g] if x[4]=='idiomorph_gene_only')
    rows.append((m['genus'],m['status'],typ,conc,g,m['species'],m['strain'],m['busco'],m['n50'],S[1]+':'+str(S[5]) if S else '',max([h[2] for h in cn],default='')))
with open('classified.tsv','w') as o:
    o.write('genus\tstatus\ttype\tconc\tgenome\tspecies\tstrain\tbusco\tn50\tsla2\tbest_core_pid\n')
    for r in rows: o.write('\t'.join(map(str,r))+'\n')
c=collections.Counter((r[0],r[1],r[2]) for r in rows)
for k in sorted(c): print(k,c[k])
print(collections.Counter((r[1],r[3]) for r in rows))
