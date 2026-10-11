"""
realistic_acc.py -- 把扩增失败算进属准确率
每条序列只在它能扩出的扩增子上有数据；两条序列的距离 = 共同有数据位点上的错配率（共同位点 < MINSHARED 视为不可比）。
没有任何扩增子扩出的序列算“未鉴定”。准确率 = 正确数 / 全部查询序列（未鉴定算错），同时给出已鉴定序列中的准确率。
"""
import pickle, json, itertools, numpy as np, pandas as pd, warnings, sys
import primer_design as pdz, iterate_5R as it
warnings.filterwarnings('ignore')
SC=sys.argv[1]
S=pickle.load(open(SC+'S.pkl','rb')); B=pickle.load(open(SC+'big_labeled.pkl','rb')); sets=json.load(open(SC+'final5b_sets.json')); sets5=pickle.load(open(SC+'iter5R_sets.pkl','rb'))
tax=B['tax']; ids=B['ids']; A0=B['aln']
gen0=np.array([tax[i][5] for i in ids])
ph0=np.array([tax[i][1] if tax[i][1] else 'unclassified' for i in ids])
# 与 prep_search 相同的过滤
full=(A0[:,19:1510]!=0).all(1); idx=np.flatnonzero(full)
cnt=pd.Series(gen0[idx]).value_counts(); ok=np.isin(gen0[idx],cnt[cnt>=2].index); idx=idx[ok]
# 排除新设计用过的那一半（Gd）
gi=np.random.default_rng(1).permutation(len(ids)); design_ids=set(gi[:len(ids)//2])
held=np.array([i for i in idx if i not in design_ids])
cnt=pd.Series(gen0[held]).value_counts(); held=held[np.isin(gen0[held],cnt[cnt>=2].index)]
print('held-out labelled',len(held),'genera',len(set(gen0[held])),flush=True)
A=A0[held]; gen=gen0[held]
G=pdz.Reference(A,[ids[i] for i in held],np.array(['Bacteria']*len(held)),ph0[held],S.col_of_pos,S.ref_seq,{})
x=np.full(A.shape,4,np.int64)
for k,v in {ord('A'):0,ord('C'):1,ord('G'):2,ord('T'):3}.items(): x[A==k]=v
rows_all=np.arange(len(held))
def hit_site(st,L,tops,orient):
    M=pdz.site_matrix(G,st,L); ok=(M!=0).all(1)&(M!=45).all(1); h=np.zeros(len(M),bool)
    h[ok]=it.hit_any(M[ok],[p if orient=='F' else pdz.revcomp(p) for p in tops],'right' if orient=='F' else 'left'); return h
FIVE=[('R1-F','R1-R',103,314,18,19),('R2-F','R2-R',338,519,18,18),('R3-F','R3-R',685,908,20,20),('R4-F','R4-R',944,1087,21,18),('R5-F','R5-R',1175,1374,19,18)]
FIVE=[('R1-F','R1-R',103,314,18,19),('R2-F','R2-R',338,519,18,18),('R3-F','R3-R',685,908,18,20),('R4-F','R4-R',944,1087,21,18),('R5-F','R5-R',1175,1374,19,18)]
def design_5R(cfg):
    out=[]
    for fn,rn,fs,rs,fl,rl in FIVE:
        hf=hit_site(fs,fl,sets5[(fn,cfg)],'F'); hr=hit_site(rs,rl,sets5[(rn,cfg)],'R'); out.append(((fs+fl,rs-1),hf&hr))
    return out
def design_new(which,keep=(1,2,3,4,5)):
    out=[]
    for k in keep:
        f=sets[f'A{k}-F']; r=sets[f'A{k}-R']
        out.append(((f['st']+f['L'],r['st']-1),hit_site(f['st'],f['L'],f[which],'F')&hit_site(r['st'],r['L'],r[which],'R')))
    return out
def realistic(design,minshared=60,chunk=500):
    n=len(x); Pres=np.zeros((n,x.shape[1]),bool)
    for (a,b),h in design: Pres[h,a-1:b]=True
    ident=Pres.any(1)
    Sm=np.zeros((n,n),np.float32); Lm=np.zeros((n,n),np.float32)
    cols=np.flatnonzero(Pres.any(0))
    for s0 in range(0,len(cols),300):
        c=cols[s0:s0+300]; Pm=Pres[:,c].astype(np.float32)
        Lm+=Pm@Pm.T
        X=np.zeros((n,len(c)*4),np.float32)
        for k,p in enumerate(c):
            m=Pres[:,p]&(x[:,p]<4); X[np.where(m)[0],k*4+x[m,p]]=1
        Sm+=X@X.T
    rate=np.where(Lm>=minshared,(Lm-Sm)/np.maximum(Lm,1),np.inf)
    np.fill_diagonal(rate,np.inf)
    mn=rate.min(1,keepdims=True); tied=(rate==mn)&np.isfinite(rate)
    same=(gen[None,:]==gen[:,None])
    credit=np.where(tied.any(1),(tied&same).sum(1)/np.maximum(tied.sum(1),1),0.0)
    credit=credit*ident
    valid=pd.Series(gen).map(pd.Series(gen).value_counts()).values>=2
    q=valid
    return dict(identified=float(ident[q].mean()),acc_identified=float(credit[q&ident].mean()),acc_overall=float(credit[q].mean()),
                mean_amps=float(sum(h for _,h in design).mean()) if False else float(np.sum([h for _,h in design],axis=0).mean()))
rows=[]
cfgs=[('5R 现有引物',design_5R('原引物')),('5R 位点，每位点 3 条（迭代引物）',design_5R('3条 各≤4')),
      ('新 5 扩增子，每位点 1 条',design_new('single')),('新 5 扩增子，每位点 2–3 条',design_new('multi')),
      ('新 4 扩增子（去 V1V2），2–3 条',design_new('multi',(2,3,4,5))),('新 3 扩增子（V3+V4+V8V9），2–3 条',design_new('multi',(2,3,5)))]
for nm,d in cfgs:
    r=realistic(d); r['design']=nm; rows.append(r); print(r,flush=True)
pd.DataFrame(rows).to_csv(SC+'realistic_acc.csv',index=False)
