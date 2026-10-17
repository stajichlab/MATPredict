import os,glob,yaml,collections,csv,sys
L=getattr(yaml,"CSafeLoader",yaml.SafeLoader)
want={l.strip() for l in open(sys.argv[1])}
P='/bigdata/stajichlab/jstajich/projects/core-hmms'
out=[]
for run,root in (('v060',P+'/asco'),('records',P+'/dothideo')):
    for sub in os.listdir(root):
        rd=f'{root}/{sub}/runs'
        if not os.path.isdir(rd): continue
        for g in os.listdir(rd):
            if g not in want: continue
            doc=yaml.load(open(f'{rd}/{g}/detection_report.yaml'),Loader=L)
            loci=[c for c in doc.get('detected') or [] if c.get('family')=='Ascomycota:MAT']
            if not loci: out.append(dict(run=run,genome=g,nloci=0)); continue
            for c in loci:
                res=c.get('idiomorph_resolutions') or []
                out.append(dict(run=run,genome=g,nloci=len(loci),contig=c['contig'],start=c['start'],end=c['end'],label=c['idiomorph'],conf=c['confidence'],lclass=c['locus_class'],
                  cands='|'.join(f"{x['idiomorph']}:{x['score']}" for x in c.get('idiomorph_candidates') or []),margin=c.get('idiomorph_margin'),
                  res='|'.join(f"{r['winner']}>{r['loser']}({r['basis']},{r['winner_identity']}v{r['loser_identity']})" for r in res),
                  genes='|'.join(c.get('genes_found') or []),missing='|'.join(c.get('genes_missing') or []),span=c['end']-c['start'],
                  ambiguous=','.join(map(str,c.get('ambiguous_with') or [])), cls=c.get('idiomorph_classifier'),pass_=c.get('detection_pass')))
keys=[]
for o in out:
    for k in o:
        if k not in keys: keys.append(k)
w=csv.DictWriter(open(sys.argv[2],'w'),fieldnames=keys,delimiter='\t');w.writeheader();w.writerows(out)
print(len(out))
