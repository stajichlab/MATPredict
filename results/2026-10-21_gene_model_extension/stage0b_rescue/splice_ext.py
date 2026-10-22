"""Splice-aware 3' extension test (minimal). For each truncated alignment, from the alignment edge enumerate GT..AG introns (40-400 bp) that keep the frame open to a stop beyond the
unspliced stop; compare the chosen candidate's end with the deposited gene end. Deterministic, sequence + curated length range only. python splice_ext.py"""
import csv, numpy as np, pandas as pd
aa="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"; B="TCAG"
C={B[i//16]+B[(i//4)%4]+B[i%4]:aa[i] for i in range(64)}
comp=str.maketrans("ACGT","TGCA"); rc=lambda s:s.translate(comp)[::-1]
def rd(f):
    d,k={},None
    for l in open(f):
        l=l.strip()
        if l.startswith(">"): k=l[1:]; d[k]=""
        elif k: d[k]+=l
    return d
W=rd('truth_windows.fa'); meta={r['key']:r for r in csv.DictReader(open('truth_windows.tsv'),delimiter='\t')}
truth={f"{r['record']}|{r['gene']}":r for r in csv.DictReader(open('truth.tsv'),delimiter='\t')}
refs=rd('refs.faa'); L={}
for h,s in refs.items(): L.setdefault(h.split('|')[1],[]).append(len(s))
rng={g:(np.percentile(v,10),np.percentile(v,90),np.median(v)) for g,v in L.items()}
t=pd.read_csv('truncated_rows3.tsv',sep='\t')
t['outcome']=np.where(t.delta_aa.abs()<0.5,'exact',np.where(t.delta_aa.abs()<=5,'within5',np.where(t.delta_aa<-20,'undershoot','other')))
nex=pd.read_csv('truth.tsv',sep='\t').assign(key=lambda d:d.record+'|'+d.gene).set_index('key').nexons
t['nexons']=t.key.map(nex)
MININT,MAXINT,MINREOPEN=40,400,30
def translate_until_stop(s):
    n=0
    for i in range(0,len(s)-2,3):
        c=C.get(s[i:i+3],'X')
        if c=='*': return n,True
        n+=1
    return n,False
res=[]
for _,r in t.iterrows():
    key=r.key; seq=W[key]; off=int(meta[key]['seq_start'])-1; tr=truth[key]
    strand=r.strand; mend=int(r.model_end); L_=len(seq)
    exons=sorted((int(a)-off,int(b)-off) for a,b in (e.split('-') for e in tr['exons'].split(';')))
    tsum=sum(b-a+1 for a,b in exons); inc=(tsum==3*int(tr['prot_aa'])+3)
    tend=(min(a for a,_ in exons) if strand=='-' else max(b for _,b in exons))
    tend_ns=(tend+3 if strand=='-' else tend-3) if inc else tend
    if strand=='-':           # work in transcript orientation
        s=rc(seq); E=L_-mend+1; TE=L_-tend_ns+1
    else:
        s=seq; E=mend; TE=tend_ns
    # unspliced extension: codons from E+1
    n0,stop0=translate_until_stop(s[E:E+3000])
    un_end=E+3*n0
    cands=[]
    # donors: exon ends at p-1 (p = first intron base, 1-based); need GT at p,p+1 and p in (E, E+3*n0+3]
    for p in range(E+1,min(E+3*n0+3,len(s)-1)+1):
        if s[p-1:p+1]!='GT': continue
        for a in range(p+MININT-2,min(p+MAXINT-2,len(s)-2)+1):   # AG at a,a+1 (1-based); intron p..a+1; length a+2-p
            if s[a-1:a+1]!='AG': continue
            ex1=s[E:p-1]; rest=s[a+1:a+1+3000]
            cds=ex1+rest
            n,stop=translate_until_stop(cds)
            if not stop: continue
            # frame must reopen well beyond the unspliced stop and the first part must have no stop (guaranteed by p<=stop0 position)
            if '*' in ''.join(C.get(ex1[i:i+3],'X') for i in range(0,len(ex1)-len(ex1)%3,3)): continue
            if n < n0+MINREOPEN: continue
            # end of the CDS (last coding base before stop) in transcript coordinates
            used=n*3-len(ex1)         # bases consumed from the second exon
            cands.append(dict(p=p,a=a,intron=a+2-p,n=n,end=a+1+used))
    total_un=int(r.model_aa)+n0
    lo,hi,med=rng[r.gene]
    def total(c): return int(r.model_aa)+c['n']
    ok=[c for c in cands if lo<=total(c)<=hi]
    chosen=None
    if ok: chosen=min(ok,key=lambda c:abs(total(c)-med))     # closest to the curated median length
    d_un=(un_end-TE)/3
    d_ch=None if chosen is None else (chosen['end']-TE)/3
    res.append(dict(key=key,frac=r.frac,radius=r.radius,nexons=r.nexons,order=r.order,truth_aa=r.truth_aa,total_unspliced=total_un,delta_unspliced=round(d_un,1),
                    n_cands=len(cands),n_cands_in_range=len(ok),chosen_intron=None if chosen is None else chosen['intron'],delta_chosen=None if d_ch is None else round(d_ch,1),
                    truth_reachable=any(abs((c['end']-TE)/3)<0.5 for c in cands),truth_reachable_in_range=any(abs((c['end']-TE)/3)<0.5 for c in ok),outcome_unspliced=r.outcome))
o=pd.DataFrame(res); o.to_csv('splice_ext_out.tsv',sep='\t',index=False)
print('rows',len(o))
print('rows with >=1 candidate:',int((o.n_cands>0).sum()),' with >=1 in length range:',int((o.n_cands_in_range>0).sum()))
print('candidates per row (median, 90th pct):',o.n_cands.median(),o.n_cands.quantile(.9),'  in range:',o.n_cands_in_range.median(),o.n_cands_in_range.quantile(.9))
print('true end reachable by some candidate:',int(o.truth_reachable.sum()),' by one inside the length range:',int(o.truth_reachable_in_range.sum()))
bad=o[o.outcome_unspliced=='undershoot']; print('\nthe %d early-stop rows of the unspliced extension:'%len(bad))
print(' with a candidate in range:',int((bad.n_cands_in_range>0).sum()),' truth reachable in range:',int(bad.truth_reachable_in_range.sum()),' chosen candidate exact (|d|<0.5):',int((bad.delta_chosen.abs()<0.5).sum()),' within 5 aa:',int((bad.delta_chosen.abs()<=5).sum()))
good=o[o.outcome_unspliced.isin(['exact','within5'])]; print('\nthe %d rows the unspliced extension got right:'%len(good))
print(' a spliced candidate in range exists for',int((good.n_cands_in_range>0).sum()),'(a competing, wrong alternative); chosen candidate would REPLACE the right answer with a wrong one in',int(((good.delta_chosen.abs()>5)).sum()))
print('\nchosen vs unspliced over ALL rows (rule: use the spliced candidate only when exactly... see below)')
o['final']=np.where(o.delta_chosen.notna(),o.delta_chosen,o.delta_unspliced)
print(' always take the chosen spliced candidate when one exists: exact',int((o.final.abs()<0.5).sum()),'within5',int((o.final.abs()<=5).sum()),'of',len(o),' | unspliced alone: exact',int((o.delta_unspliced.abs()<0.5).sum()),'within5',int((o.delta_unspliced.abs()<=5).sum()))
