import csv,collections
R="/bigdata/stajichlab/jstajich/projects/MATPredict/results/sra_reads_typing/wave_inputs"
rd=lambda p:list(csv.DictReader(open(p,newline=""),delimiter="\t"))
runs=[r for r in rd(R+"/fungi_wgs_illumina.tsv") if r["scientific_name"].startswith("Botrytis cinerea")]
print("runs:",len(runs),"| samples:",len({r["sample_accession"] for r in runs}),"| studies:",len({r["study_accession"] for r in runs}),"| names:",collections.Counter(r["scientific_name"] for r in runs).most_common(4))
S={r["study"]:r for r in rd(R+"/studies.tsv") if r["species"]=="Botrytis cinerea"}
print("studies table rows for the species:",len(S))
by=collections.defaultdict(list)
for r in runs: by[r["study_accession"]].append(r)
flag={k:v.get("excluded_as","") for k,v in S.items()}
tot_ok=0
for st,v in sorted(by.items(),key=lambda kv:-len(kv[1])):
    gb=[int(x["base_count"] or 0)/1e9 for x in v]
    print(f"  {st:14s} runs {len(v):4d} samples {len({x['sample_accession'] for x in v}):4d} flagged={flag.get(st,'(not in studies table)')!r:14} Gb min/median/max {min(gb):.1f}/{sorted(gb)[len(gb)//2]:.1f}/{max(gb):.1f}")
