import gzip,subprocess,tempfile,itertools
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
aa='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
C={a+b+c:aa[i*16+j*4+k] for i,a in enumerate('TCAG') for j,b in enumerate('TCAG') for k,c in enumerate('TCAG')}
tr=lambda s:''.join(C.get(s[i:i+3].upper(),'X') for i in range(0,len(s)-2,3))
# exonerate's intron (forward genomic 1397174..1397226) in window coords, 1-based
a,b=1397174-ws+1,1397226-ws+1
print('exonerate intron forward',seq[a-1:b][:4],'...',seq[a-1:b][-4:],'len',b-a+1,' transcript-sense (minus):',rc(seq[a-1:b])[:4],'...',rc(seq[a-1:b])[-4:])
td=tempfile.mkdtemp(); open(td+'/t.fa','w').write('>w\n'+seq+'\n'); open(td+'/q.fa','w').write('>q\n'+refs['5016_c5_MAT_MAT1-1']+'\n')
def run(spsc=None,extra=()):
    cmd=['miniprot','-t1','--gff','-G','400']+list(extra)+(['--spsc='+spsc] if spsc else [])+[td+'/t.fa',td+'/q.fa']
    o=subprocess.run(cmd,capture_output=True,text=True)
    ex=[];strd='+'
    for l in o.stdout.splitlines():
        f=l.split('\t')
        if len(f)>8 and f[2]=='CDS': ex.append((int(f[3]),int(f[4]))); strd=f[6]
    if not ex: return 'no model '+o.stderr[-100:].replace('\n',' ')
    ex.sort(); order=ex[::-1] if strd=='-' else ex
    cds=''.join(rc(seq[x-1:y]) if strd=='-' else seq[x-1:y] for x,y in order); p=tr(cds)
    ints=[]
    for i in range(len(ex)-1):
        s=seq[ex[i][1]:ex[i+1][0]-1]; s=rc(s) if strd=='-' else s
        ints.append(f'{s[:2]}..{s[-2:]}({len(s)}bp,{ex[i][1]+ws}-{ex[i+1][0]-2+ws})')
    return f'{strd} exons={len(ex)} aa={len(p)} stops={p[:-1].count("*")} introns={ints}'
print('baseline', run())
# conventions for the minus-strand donor (transcript GC) / acceptor (AG): forward intron is a..b
for name,(donpos,accpos) in {'D at last intron base 0-based (b-1), A at first intron base 0-based (a-1)':(b-1,a-1),
                              'D at b, A at a-1':(b,a-1),'D at b-1, A at a':(b-1,a),'D at b, A at a':(b,a),
                              'D at b-2, A at a-1':(b-2,a-1),'D at b-1, A at a-2':(b-1,a-2)}.items():
    f=tempfile.mktemp()
    open(f,'w').write(f'w\t{donpos}\t-\tD\t10\nw\t{accpos}\t-\tA\t10\n')
    print(name,'->',run(f))
