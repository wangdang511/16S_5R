import pickle, json, sys, numpy as np, pandas as pd, time, warnings
import tiling_search as ts
warnings.filterwarnings('ignore')
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
D=pickle.load(open(SC+'search_data.pkl','rb')); x=D['x']; gen=D['gen']; d_idx=D['d_idx']; v_idx=D['v_idx']
xd,gd=x[d_idx],gen[d_idx]; xv,gv=x[v_idx],gen[v_idx]
vd=(pd.Series(gd).map(pd.Series(gd).value_counts()).values>=2); vv=(pd.Series(gv).map(pd.Series(gv).value_counts()).values>=2)
C=ts.cumulative_matches(xd)
m=pd.read_csv(SC+'sites_joint.csv')
KMIN,KMEAN,CAP,tag,KMAX=float(sys.argv[1]),float(sys.argv[2]),int(sys.argv[3]),sys.argv[4],int(sys.argv[5])
MINLEN=int(sys.argv[6]) if len(sys.argv)>6 else 150
s=m[(m.key_min>=KMIN)&(m.key_mean>=KMEAN)].sort_values('key_min',ascending=False).drop_duplicates(['start','orient'])
F=s[s.orient=='F'].sort_values('start').reset_index(drop=True); R=s[s.orient=='R'].sort_values('start').reset_index(drop=True)
Fs,Fl,Rs,Rl=F.start.values,F.len.values,R.start.values,R.len.values
GAP=20
pairs=[]
for f in range(len(Fs)):
    for r in range(len(Rs)):
        L=Rs[r]+Rl[r]-1-Fs[f]+1
        if MINLEN<=L<=CAP and Rs[r]>Fs[f]+Fl[f]+60: pairs.append((f,r))
pairs=np.array(pairs); print(tag,'sites',len(Fs),len(Rs),'pairs',len(pairs),flush=True)
sp=np.array([(Fs[f],Rs[r]+Rl[r]-1) for f,r in pairs]); iv=np.array([(Fs[f]+Fl[f],Rs[r]-1) for f,r in pairs])
def ok_with(others):
    m_=np.ones(len(pairs),bool)
    for a,b in others: m_&=(sp[:,1]+GAP<a)|(sp[:,0]-GAP>b)
    return m_
NNo=ts.NN(gd,vd)
def score(S,L): return NNo.acc(L-S)
def best_for_slot(others_idx):
    others=[sp[i] for i in others_idx]
    S0=0;L0=0
    for i in others_idx: a,b=iv[i]; S0=S0+(C[b].astype(np.int32)-C[a-1]); L0+=b-a+1
    cand=np.flatnonzero(ok_with(others)); best=(-1,None)
    for i in cand:
        a,b=iv[i]; S=(S0+(C[b].astype(np.int32)-C[a-1])) if len(others_idx) else (C[b].astype(np.int32)-C[a-1])
        sc=score(S,L0+b-a+1)
        if sc>best[0]: best=(sc,i)
    return best
def fmt(st): return ' '.join(f'{sp[i][0]}-{sp[i][1]}({sp[i][1]-sp[i][0]+1})' for i in sorted(st,key=lambda i:sp[i][0]))

rng=np.random.default_rng(int(sys.argv[7]))
K=int(sys.argv[5]); NSEED=int(sys.argv[8])
def seed_state(mode):
    if mode=='prev':
        prev=[(11,260),(337,530),(558,802),(910,1195),(1222,1509)][:K] if K<=5 else None
        if prev is None: return None
        st=[]
        for a,b in prev:
            d=np.abs(sp[:,0]-a)+np.abs(sp[:,1]-b); st.append(int(d.argmin()))
        return st if len(set(st))==K and all(ok_with([sp[j] for j in st if j!=i])[i] for i in st) else None
    for _ in range(500):
        edges=np.linspace(1,1542,K+1); st=[]
        for k in range(K):
            c=np.flatnonzero((sp[:,0]>=edges[k]-30)&(sp[:,1]<=edges[k+1]+30))
            if len(c)==0: break
            st.append(int(rng.choice(c)))
        else:
            S_=sorted(tuple(sp[i]) for i in st)
            if all(S_[j+1][0]-S_[j][1]-1>=GAP for j in range(len(S_)-1)): return st
    return None
results=[]
modes=['prev']+['rand']*NSEED
for mode in modes:
    state=seed_state(mode)
    if state is None: print('seed',mode,'failed',flush=True); continue
    ivs=[tuple(iv[i]) for i in state]
    S0=sum((C[b].astype(np.int32)-C[a-1]) for a,b in ivs); sc=score(S0,sum(b-a+1 for a,b in ivs))
    print(mode,'seed score',round(sc,4),fmt(state),flush=True)
    for sweep in range(4):
        changed=False
        for pos in range(len(state)):
            others=state[:pos]+state[pos+1:]
            sc2,j=best_for_slot(others)
            if j is not None and sc2>sc+1e-9: state[pos]=j; sc=sc2; changed=True
        print('  sweep',sweep,round(sc,4),fmt(state),flush=True)
        if not changed: break
    ivs=[tuple(iv[i]) for i in state]; va=ts.accuracy_direct(xv,ivs,gv,vv)
    minq=min(min(F.key_min[pairs[i][0]],R.key_min[pairs[i][1]]) for i in state)
    row=dict(tag=tag,K=K,mode=mode,design_acc=float(sc),val_acc=float(va),bases=int(sum(b-a+1 for a,b in ivs)),minq=float(minq),amplicons=fmt(state),
             sites=json.dumps([[int(Fs[pairs[i][0]]),int(Fl[pairs[i][0]]),F.primer[pairs[i][0]],int(Rs[pairs[i][1]]),int(Rl[pairs[i][1]]),R.primer[pairs[i][1]]] for i in sorted(state,key=lambda i:sp[i][0])]))
    results.append(row); print({k:v for k,v in row.items() if k!='sites'},flush=True)
    pd.DataFrame(results).to_csv(SC+f'seeded_{tag}_K{K}.csv',index=False)
