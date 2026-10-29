import json,os,csv
W="/bigdata/stajichlab/jstajich/projects/MATPredict/results/sra_reads_typing/waves/wave2/plan"
for s in "aspergillus_westerdijkiae bipolaris_zeicola fusarium_ambrosium fusarium_sambucinum fusarium_virguliforme morchella_sextelata candida_ipomoeae candidozyma_pseudohaemuli meyerozyma_caribbica sclerotinia_sclerotiorum ceratocystis_fimbriata".split():
    d=f"{W}/{s}"; st=open(d+"/status").read().strip()
    log=open(d+"/plan.log").read().splitlines()
    head=[l for l in log if "genomes," in l or "borrowed" in l or "loci with a core" in l]
    qc=json.load(open(d+"/qc_panel.json")) if os.path.exists(d+"/qc_panel.json") else {}
    lt=[f for f in os.listdir(d+"/panel") if f.endswith("_lineages.tsv")] if os.path.isdir(d+"/panel") else []
    rows=list(csv.DictReader(open(d+"/panel/"+lt[0]),delimiter="\t")) if lt else []
    names=qc.get("idiomorphs") or (open(d+"/panel/"+[f for f in os.listdir(d+"/panel") if f.endswith("_idiomorphs.txt")][0]).read().split() if os.path.isdir(d+"/panel") and any(f.endswith("_idiomorphs.txt") for f in os.listdir(d+"/panel")) else [])
    per={n:sum(int(r[n]) for r in rows) for n in names} if rows else {}
    print(f"{s}: {st} | {' ; '.join(head)[:170]} | loci per idiomorph {per} | lineages {len(rows)} | L1 {json.dumps(qc.get('lineages',{}).get('L1',''))[:120]}")
