import pickle, json, numpy as np, pandas as pd, time, warnings
import tiling_search as ts
warnings.filterwarnings('ignore')
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
B=pickle.load(open(SC+'big_labeled.pkl','rb')); tax=B['tax']; ids=B['ids']; A=B['aln']
gen=np.array([tax[i][5] for i in ids]); fam=np.array([tax[i][4] if tax[i][4] else '' for i in ids])
full=(A[:,19:1510]!=0).all(1); A,gen,fam=A[full],gen[full],fam[full]
cnt=pd.Series(gen).value_counts(); ok=np.isin(gen,cnt[cnt>=2].index); A,gen,fam=A[ok],gen[ok],fam[ok]
print('labelled full-length',len(A),'genera',len(set(gen)))
rng=np.random.default_rng(42); perm=rng.permutation(len(A))
nd=1100
d_idx=perm[:nd]; v_idx=perm[nd:]
x=ts.encode(A)
pickle.dump(dict(x=x,gen=gen,fam=fam,d_idx=d_idx,v_idx=v_idx),open(SC+'search_data.pkl','wb'))
xd,gd=x[d_idx],gen[d_idx]; xv,gv=x[v_idx],gen[v_idx]
vd=(pd.Series(gd).map(pd.Series(gd).value_counts()).values>=2); vv=(pd.Series(gv).map(pd.Series(gv).value_counts()).values>=2)
print('design',len(xd),'valid-queries',vd.sum(),'| validation',len(xv),vv.sum())
t=time.time(); C=ts.cumulative_matches(xd); print('cumulative built',C.shape,C.nbytes/1e9,'GB',time.time()-t)
def acc_c(iv): 
    S=sum(C[b].astype(np.int32)-C[a-1] for a,b in iv); L=sum(b-a+1 for a,b in iv); return ts.nn_accuracy(L-S,gd,vd)
def acc_v(iv): return ts.accuracy_direct(xv,iv,gv,vv)
# baselines (interior intervals)
J=json.load(open(SC+'final5b_sets.json'))
prev=[(J[f'A{k}-F']['st']+J[f'A{k}-F']['L'],J[f'A{k}-R']['st']-1) for k in range(1,6)]
FIVE=[(103,18,314,19),(338,18,519,18),(685,18,908,20),(944,21,1087,18),(1175,19,1374,18)]
r5=[(fs+fl,rs-1) for fs,fl,rs,rl in FIVE]
base={'5R 现有':r5,'上一轮 5 扩增子':prev,'仅 V4 (515F–806R)':[(534,786)],'整条 16S(20–1510)':[(20,1510)]}
rows=[]
for nm,iv in base.items():
    t=time.time(); a=acc_c(iv); b=acc_v(iv); rows.append(dict(design=nm,bases=sum(y-x_+1 for x_,y in iv),design_acc=a,val_acc=b)); print(rows[-1],round(time.time()-t,1),flush=True)
pd.DataFrame(rows).to_csv(SC+'search_baselines.csv',index=False)
