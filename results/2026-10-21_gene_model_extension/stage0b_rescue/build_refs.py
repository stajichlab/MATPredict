import glob, yaml, os, sys
db=sys.argv[1]; n=0
with open('refs.faa','w') as out:
    for md in sorted(glob.glob(f"{db}/*/*/*/metadata.yaml")):
        if "/candidates/" in md or "/classifiers/" in md: continue
        d=yaml.safe_load(open(md)); sp=((d.get('organism') or {}).get('species') or 'NA').replace(' ','_')
        fa=os.path.join(os.path.dirname(md),'proteins.faa'); cur=None; P={}
        for l in open(fa):
            if l[0]=='>':
                cur=[x.split('=')[1] for x in l[1:].split('|') if x.startswith('name=')][0]; P[cur]=''
            elif cur: P[cur]+=l.strip()
        for g,p in P.items():
            if g in ('MAT1-1-1','MAT1-2-1') and len(p)>=60:
                out.write(f">{d['record_id']}|{g}|{sp}\n{p}\n"); n+=1
print(n,'reference proteins')
