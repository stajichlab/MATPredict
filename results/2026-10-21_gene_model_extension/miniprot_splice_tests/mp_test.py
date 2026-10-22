import gzip,subprocess,tempfile,sys,re
G='/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes/GCA_003231325.2_CUR_Ptr134_2.0.fa.gz'
contig,ws,we='CM034231.1',1395500,1398500
seq=[];on=False
for l in gzip.open(G,'rt'):
    if l[0]=='>':
        if on: break
        on=l[1:].split()[0]==contig
    elif on: seq.append(l.strip())
seq=''.join(seq)[ws-1:we]
refs={};k=None
for l in open('/bigdata/stajichlab/jstajich/projects/hl_refs.faa'):
    if l[0]=='>': k=l[1:].split('|')[0]; refs[k]=''
    else: refs[k]+=l.strip()
comp=str.maketrans('ACGTacgt','TGCAtgca'); rc=lambda s:s.translate(comp)[::-1]
B=dict(zip('TCAG','TCAG'));aa='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
C={a+b+c:aa[i*16+j*4+k] for i,a in enumerate('TCAG') for j,b in enumerate('TCAG') for k,c in enumerate('TCAG')}
def tr(s): return ''.join(C.get(s[i:i+3].upper(),'X') for i in range(0,len(s)-2,3))
variants=[('default -j1',[]),('-j0',['-j','0']),('-j2',['-j','2']),('-C0 (ignore splice signals)',['-C','0']),('-C0.3',['-C','0.3']),('-J15',['-J','15']),('-J15 -C0.3',['-J','15','-C','0.3']),('-F50 (stop penalty up)',['-F','50']),('-j0 -F50',['-j','0','-F','50'])]
td=tempfile.mkdtemp(); open(td+'/t.fa','w').write('>w\n'+seq+'\n')
for qn in ('5016_c5_MAT_MAT1-1','13684_sn435pl98_MAT_MAT1-1'):
    open(td+'/q.fa','w').write('>q\n'+refs[qn]+'\n')
    print('QUERY',qn,len(refs[qn]))
    for name,extra in variants:
        o=subprocess.run(['miniprot','-t1','--gff','-G','400']+extra+[td+'/t.fa',td+'/q.fa'],capture_output=True,text=True).stdout
        ex=[];strd='+'
        for l in o.splitlines():
            f=l.split('\t')
            if len(f)>8 and f[2]=='CDS': ex.append((int(f[3]),int(f[4]))); strd=f[6]
            if len(f)>8 and f[2]=='mRNA': m=f[8]
        if not ex: print(f'  {name:32s} no model'); continue
        ex.sort(); order=ex[::-1] if strd=='-' else ex
        cds=''.join(rc(seq[a-1:b]) if strd=='-' else seq[a-1:b] for a,b in order)
        p=tr(cds)
        introns=[(ex[i][1]+1+ws-1,ex[i+1][0]-1+ws-1) for i in range(len(ex)-1)]
        sites=[]
        for a,b in [(ex[i][1]+1,ex[i+1][0]-1) for i in range(len(ex)-1)]:
            s=seq[a-1:b]; s=rc(s) if strd=='-' else s
            sites.append(s[:2]+'..'+s[-2:]+f'({len(s)}bp)')
        print(f'  {name:32s} {strd} exons={len(ex)} cds={len(cds)} aa={len(p)} internal_stops={p[:-1].count("*")} introns={sites} genomic={[(a+ws-1,b+ws-1) for a,b in ex]}')
