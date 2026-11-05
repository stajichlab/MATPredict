import glob,yaml,statistics
O="/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-11-05_timing/"
def load(n):
    d={}
    for f in glob.glob(O+n+"/runs/*/detection_report.yaml"):
        r=yaml.safe_load(open(f)); d[f.split("/")[-2]]=(r["run"]["wall_seconds"],r["routing_mode"],r["run"]["genome"]["length_bp"]//10**6,r["run"]["genome"]["contigs"])
    return d
a=load("v070"); b=load("v080"); c=[g for g in a if g in b]; print(len(a),len(b),len(c))
for g in sorted(c,key=lambda g:b[g][1]): print("%-34s %-14s %4dMb %6d contigs  v070 %5.0fs  v080 %5.0fs  x%.1f"%(g[:34],b[g][1],b[g][2],b[g][3],a[g][0],b[g][0],b[g][0]/max(a[g][0],1)))
print("median v070 %.0f v080 %.0f ratio %.2f"%(statistics.median(a[g][0] for g in c),statistics.median(b[g][0] for g in c),statistics.median(b[g][0]/max(a[g][0],1) for g in c)))
