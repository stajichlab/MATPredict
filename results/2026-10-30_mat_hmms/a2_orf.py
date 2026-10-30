import sys
sys.argv=['x']
exec(open('score_models.py').read().split('MODELS =')[0])
def stop_down(k):
    while codon[sub(k,k+2)]!='*': k+=3
    return k+2
def stop_up(k):
    while True:
        k-=3
        if codon[sub(k,k+2)]=='*': return k
for start in (537367,):
    e=stop_down(start); u=stop_up(start)
    j=u+3
    while sub(j,j+2)!='ATG': j+=3
    print('a2 plus: upstream stop',u,'first ATG',j,'stop ends',e,'aa',(e-j+1)//3-1)
    print(tr(sub(j,e)))
