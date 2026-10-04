import pickle, json, itertools, numpy as np, pandas as pd, warnings
import tiling_search as ts
warnings.filterwarnings('ignore')
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
D=pickle.load(open(SC+'search_data.pkl','rb')); x=D['x']; gen=D['gen']; d_idx=D['d_idx']; v_idx=D['v_idx']
xd,gd=x[d_idx],gen[d_idx]; xv,gv=x[v_idx],gen[v_idx]
vd=(pd.Series(gd).map(pd.Series(gd).value_counts()).values>=2); vv=(pd.Series(gv).map(pd.Series(gv).value_counts()).values>=2)
J=json.load(open(SC+'final5b_sets.json'))
iv={f'A{k}':(J[f'A{k}-F']['st']+J[f'A{k}-F']['L'],J[f'A{k}-R']['st']-1) for k in range(1,6)}
names={'A1':'V1V2','A2':'V3','A3':'V4','A4':'V6V7','A5':'V8V9'}
print(iv)
# 对所有验证+设计序列合并后评估（样本更多，噪声更小）：用全部 5899 条序列
xa=x; ga=gen; va=(pd.Series(ga).map(pd.Series(ga).value_counts()).values>=2)
rows=[]
for r in range(1,6):
    for comb in itertools.combinations(iv,r):
        ints=[iv[c] for c in comb]
        acc_all=ts.accuracy_direct(xa,ints,ga,va)
        rows.append(dict(combo='+'.join(names[c] for c in comb),n=r,bases=sum(b-a+1 for a,b in ints),acc_all=acc_all))
df=pd.DataFrame(rows).sort_values(['n','acc_all'],ascending=[True,False]); df.to_csv(SC+'subset_acc.csv',index=False)
for n in range(1,6):
    print(df[df.n==n].head(4).round(4).to_string(index=False))
