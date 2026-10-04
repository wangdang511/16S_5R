# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换
import pickle, json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'robust_realistic.py').read().split("D={'S1'")[0]
exec(src)
sj=json.load(open(SC+'v2_tm_sites.json')); c4=store[('4 扩增子',0.95)]['sites']
A={'A1-F':sj['A1-F'],'A1-R':sj['A1-R_balanced']}
d5=build({**A,**{f'A{int(k[1])+1}-{k[3]}':v for k,v in c4.items()}})   # 第 0 个是 V1V2
d4=build(c4)
def mats(design):
    n=len(xh); Pres=np.zeros((n,xh.shape[1]),bool)
    for (a,b),h in design: Pres[h,a-1:b]=True
    Sm=np.zeros((n,n),np.float32); Lm=np.zeros((n,n),np.float32); cols=np.flatnonzero(Pres.any(0))
    for s0 in range(0,len(cols),300):
        c=cols[s0:s0+300]; Pm=Pres[:,c].astype(np.float32); Lm+=Pm@Pm.T
        X=np.zeros((n,len(c)*4),np.float32)
        for k,p in enumerate(c):
            m=Pres[:,p]&(xh[:,p]<4); X[np.where(m)[0],k*4+xh[m,p]]=1
        Sm+=X@X.T
    return Pres,Sm,Lm
def score(main,extra,minabs=60,minfrac=0.0):
    Pm,Sm,Lm=mats(main); n=len(Pm); ident=Pm.any(1)
    if extra is not None:
        Pe,Se,Le=mats(extra); ident=ident|Pe.any(1)
    need=np.maximum(minabs,minfrac*Pm.sum(1,keepdims=True))
    Rm=np.where(Lm>=need,(Lm-Sm)/np.maximum(Lm,1),np.inf)
    if extra is not None:
        Re=np.where(Le>=minabs,(Le-Se)/np.maximum(Le,1),np.inf)
    same=(genh[None,:]==genh[:,None]); q=pd.Series(genh).map(pd.Series(genh).value_counts()).values>=2
    np.fill_diagonal(Rm,np.inf)
    if extra is not None: np.fill_diagonal(Re,np.inf)
    credit=np.zeros(n)
    for i in range(n):
        r=Rm[i]
        if np.isfinite(r).any():
            cand=np.flatnonzero(r==r.min())
            if extra is not None and len(cand)>1:
                e=Re[i,cand]
                if np.isfinite(e).any(): cand=cand[e==e.min()]
        elif extra is not None and np.isfinite(Re[i]).any():
            cand=np.flatnonzero(Re[i]==Re[i].min())
        else: continue
        credit[i]=same[i,cand].mean()*ident[i]
    return float(credit[q].mean())
V=d5[0:1]; rest=d5[1:]
for lab,a,f in (('abs60',60,0),('abs300',300,0),('frac0.8',0,0.8)):
    print(lab,'4扩增子',round(score(d4,None,a,f),4),'| 5扩增子（V1V2 只用于平局/兜底）',round(score(rest,V,a,f),4),'| 5扩增子（原规则，一起算）',round(score(d5,None,a,f),4),flush=True)
