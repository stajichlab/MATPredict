"""Score every annotated protein of the in-scope LCG genomes with the
Mucoromycota sexM/sexP classifier HMMs (read-only use of the frozen db)."""
import csv, sys, pyhmmer
W='/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-076afe4/db/Mucoromycota/classifiers/MAT'
hmms={}
for g in ('sexM','sexP'):
    with pyhmmer.plan7.HMMFile(f'{W}/{g}.hmm') as f: hmms[g]=f.read()
alpha=pyhmmer.easel.Alphabet.amino()
rows=list(csv.DictReader(open('genomes.tsv'),delimiter='\t'))
out=open('annot_hmm_hits.tsv','w'); out.write('org\tprotein\tlength\tsexM_bits\tsexP_bits\tmargin_P_minus_M\tbest\n')
for x in rows:
    if x['group'] not in ('Mucorales','Umbelopsidales'): continue
    with pyhmmer.easel.SequenceFile(x['proteins'],digital=True,alphabet=alpha) as sf: seqs=sf.read_block()
    sc={}
    for g,h in hmms.items():
        for hits in pyhmmer.hmmsearch([h],seqs,cpus=8,E=10):
            for hit in hits:
                sc.setdefault(hit.name.decode() if isinstance(hit.name,bytes) else hit.name,{})[g]=hit.score
    lens={ (s.name.decode() if isinstance(s.name,bytes) else s.name):len(s) for s in seqs}
    for p,d in sc.items():
        m=d.get('sexM',0.0); pp=d.get('sexP',0.0)
        if max(m,pp)<30: continue
        out.write(f"{x['org']}\t{p}\t{lens.get(p,'')}\t{m:.1f}\t{pp:.1f}\t{pp-m:.1f}\t{'sexP' if pp>m else 'sexM'}\n")
out.close()
