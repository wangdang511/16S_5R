import pickle, numpy as np, pandas as pd, time, warnings
import tiling_search as ts
warnings.filterwarnings('ignore')
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
D=pickle.load(open(SC+'search_data.pkl','rb')); x=D['x']; gen=D['gen']; d_idx=D['d_idx']; v_idx=D['v_idx']
xd,gd=x[d_idx],gen[d_idx]; xv,gv=x[v_idx],gen[v_idx]
vd=(pd.Series(gd).map(pd.Series(gd).value_counts()).values>=2); vv=(pd.Series(gv).map(pd.Series(gv).value_counts()).values>=2)
C=ts.cumulative_matches(xd); NNo=ts.NN(gd,vd)
BS=20; edges=list(range(20,1510,BS)); blocks=[(a,min(a+BS-1,1510)) for a in edges]
Bm=[C[b].astype(np.int16)-C[a-1] for a,b in blocks]; Bl=[b-a+1 for a,b in blocks]
nb=len(blocks); keep=set(range(nb))
S=sum(m.astype(np.int32) for m in Bm); L=sum(Bl)
rows=[dict(n_blocks=nb,bases=L,acc=NNo.acc(L-S),removed=-1)]
print('all blocks',rows[0],flush=True)
while len(keep)>1:
    best=(-1,None)
    for k in keep:
        a=NNo.acc((L-Bl[k])-(S-Bm[k]))
        if a>best[0]: best=(a,k)
    k=best[1]; keep.remove(k); S=S-Bm[k]; L-=Bl[k]
    rows.append(dict(n_blocks=len(keep),bases=L,acc=best[0],removed=k))
df=pd.DataFrame(rows); df.to_csv(SC+'block_backward.csv',index=False)
for nbk in [75,60,50,40,30,25,20,15,12,10,8,6,4]:
    r=df[df.n_blocks==nbk].iloc[0]; print(nbk,int(r.bases),round(r.acc,4))
# 保留 40 个块（800 nt）时的位置
keep40=None
S=sum(m.astype(np.int32) for m in Bm); L=sum(Bl); keep=set(range(nb))
for _,r in df.iloc[1:].iterrows():
    if r.n_blocks<=60 and len(keep)>60: pass
rem=df.removed.values[1:]
def kept_at(n):
    rs=set(rem[:nb-n]); return [i for i in range(nb) if i not in rs]
out={}
for n in [60,40,30,20]:
    ks=kept_at(n); out[n]=[blocks[i] for i in ks]
    # merge adjacent
    seg=[]; 
    for a,b in out[n]:
        if seg and a==seg[-1][1]+1: seg[-1]=(seg[-1][0],b)
        else: seg.append((a,b))
    iv=seg; va=ts.accuracy_direct(xv,iv,gv,vv)
    print(f'keep {n} blocks ({sum(b-a+1 for a,b in iv)} nt) validation acc {va:.4f}; segments:',iv)
