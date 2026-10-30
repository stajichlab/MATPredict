import csv,collections,subprocess,os,sys
W="/bigdata/stajichlab/jstajich/projects/MATPredict/results/sra_reads_typing/waves/wave2"; S="/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-26_v070_submit"
rd=lambda p:list(csv.DictReader(open(p,newline=""),delimiter="\t"))
rows=rd(W+"/final/summary/results_all.tsv"); by=collections.defaultdict(list)
for r in rows:
    if r["state"]=="done" and os.path.exists(f"{W}/out/{r['run']}.sig.zip"): by[r["group"]].append(r)
os.makedirs(W+"/clone",exist_ok=True)
py="/opt/linux/rocky/8.x/x86_64/pkgs/sourmash/4.9/bin/python"
for g,v in sorted(by.items()):
    if len(v)<10: continue
    p=subprocess.run([py,"-I",S+"/clone_clusters.py","--sigs",W+"/out","--runs",W+"/final/summary/results_all.tsv","--group",g,"--out",f"{W}/clone/{g}.clusters.tsv","--pairs",f"{W}/clone/{g}.pairs.tsv"],capture_output=True,text=True)
    print(g,len(v),"runs:",(p.stdout.strip().splitlines() or [p.stderr.strip()[-100:]])[-1])
