import sys,glob,csv,collections
meta={r['genome']:r for r in csv.DictReader(open('dipo_meta.tsv'),delimiter='\t')}
queries_core=['Ylip_AJ617306|mata1','Ylip_AJ617306|mata2','Ylip_AJ617307|matb1','Ylip_AJ617307|matb2','Gcand_HF558448|MATA','Gcand_HF558449|MATB','Ohpoly_AJ617305|matalpha1','Ohpoly_AJ617305|matalpha2','Ohpoly_AJ617305|mata1']
flank=['Ylip_AJ617306|sla2','Ylip_AJ617306|apn2','Ohpoly_AJ617305|sla2','Ohpoly_AJ617305|dic1']
for f in sorted(glob.glob('out/out_*.tsv')):
    g=f[8:-4]; m=meta[g]
    best={}
    for l in open(f):
        q,s,pid,ln,ql,bs,a,b,ev=l.rstrip().split('\t'); pid=float(pid);ln=int(ln);ql=int(ql);bs=float(bs)
        k='|'.join(q.split('|')[:2])
        if float(ev)>1e-5: continue
        h=(bs,pid,ln/ql,s,min(int(a),int(b)),max(int(a),int(b)))
        if k not in best or bs>best[k][0]: best[k]=h
    print('##',g,m['species'],m['strain'],m['status'],m['families_called'],'busco',m['busco'],'n50',m['n50'])
    fl={k:best[k] for k in flank if k in best}
    for k in queries_core+flank:
        if k in best:
            b=best[k]; near=''
            if k in queries_core:
                near=','.join(fk.split('|')[1] for fk,fv in fl.items() if fv[3]==b[3] and abs(fv[4]-b[4])<30000)
            print(f'   {k:28s} pid={b[1]:.0f} cov={b[2]:.2f} {b[3]}:{b[4]}-{b[5]} {"near:"+near if near else ""}')
