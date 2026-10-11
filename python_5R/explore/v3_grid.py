import pickle, json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore'); pd.set_option('display.width',250)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'ext2.py').read().split("rows=[]")[0]
exec(src)
D1=json.load(open(SC+'short_sites.json'))['D1']
ints0=[(D1[f'A{k}-F']['start']+D1[f'A{k}-F']['L'],D1[f'A{k}-R']['start']-1) for k in range(1,6)]
print(ints0)
base=ts.accuracy_direct(xall,ints0,gall,vall); print('base',round(base,4))
rows=[]
for s in list(range(300,358,4))+[357]:
    for e in list(range(504,540,4))+[504]:
        ints=list(ints0); ints[1]=(s,e); a=ts.accuracy_direct(xall,ints,gall,vall)
        rows.append(dict(s=s,e=e,len_interior=e-s+1,acc=a,gain=(a-base)*100))
R=pd.DataFrame(rows).drop_duplicates(['s','e']); R.to_csv(SC+'v3_grid.csv',index=False)
print(R.pivot(index='s',columns='e',values='gain').round(2).to_string())
