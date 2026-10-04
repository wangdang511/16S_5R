import pickle, json, itertools, numpy as np, pandas as pd, warnings
import tiling_search as ts
warnings.filterwarnings('ignore')
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
D=pickle.load(open(SC+'search_data.pkl','rb')); x=D['x']; gen=D['gen']; v_idx=D['v_idx']; d_idx=D['d_idx']
xv,gv=x[v_idx],gen[v_idx]; vv=(pd.Series(gv).map(pd.Series(gv).value_counts()).values>=2)
J=json.load(open(SC+'final5b_sets.json'))
iv={k:(J[f'A{k}-F']['st']+J[f'A{k}-F']['L'],J[f'A{k}-R']['st']-1) for k in range(1,6)}
names={1:'V1V2',2:'V3',3:'V4',4:'V6V7',5:'V8V9'}
oligos={k:len(J[f'A{k}-F']['multi'])+len(J[f'A{k}-R']['multi']) for k in range(1,6)}
single={k:2 for k in range(1,6)}
rows=[]
for r in range(1,6):
    for comb in itertools.combinations(range(1,6),r):
        ints=[iv[k] for k in comb]
        rows.append(dict(combo='+'.join(names[k] for k in comb),n=r,bases=sum(b-a+1 for a,b in ints),val_acc=ts.accuracy_direct(xv,ints,gv,vv),
                         spans=' '.join(f"{J[f'A{k}-F']['st']}-{J[f'A{k}-R']['st']+J[f'A{k}-R']['L']-1}" for k in comb),
                         sites=2*r,oligos_multi=sum(oligos[k] for k in comb)))
df=pd.DataFrame(rows); df.to_csv(SC+'subset_val.csv',index=False)
print(df.sort_values(['n','val_acc'],ascending=[True,False]).groupby('n').head(3).round(4).to_string(index=False))
