import glob,yaml,collections
O="/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-11-04_148_hyphopichia/"
RANK={"low":0,"medium":1,"high":2}
def load(run):
    d={}
    for f in glob.glob(O+run+"/runs/*/detection_report.yaml"):
        r=yaml.safe_load(open(f)); d[f.split("/")[-2]]=sorted((l["idiomorph"],l["contig"],l["start"],l["end"],l["confidence"]) for l in r["detected"])
    return d
S={n:load(n) for n in ("base","cand","base1b","cand1b")}
def cmp(a,b,label):
    up=down=same=other=0; downs=[]
    for g in S[a]:
        x=S[a][g]; y=S[b].get(g,[])
        if x==y: same+=1; continue
        cx=sum(RANK[l[4]] for l in x); cy=sum(RANK[l[4]] for l in y)
        if cy>cx: up+=1
        elif cy<cx: down+=1; downs.append(g[:15])
        else: other+=1
    print(f"{label:34s} same {same:3d}  up {up:2d}  down {down:2d}  other {other:2d}", downs[:12] if downs else "")
cmp("base","cand","record added (floor as is)")
cmp("base","base1b","1b alone (no record)")
cmp("base1b","cand1b","record added under 1b")
cmp("base","cand1b","record + 1b vs today")
