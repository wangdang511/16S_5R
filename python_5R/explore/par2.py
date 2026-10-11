# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/par_*
import pickle, json, itertools, time, numpy as np, pandas as pd, warnings, sys
warnings.filterwarnings('ignore'); pd.set_option('display.width',250); pd.set_option('display.max_colwidth',50)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'ten1.py').read().split("rows=[]; masks={}")[0]
exec(src)
U=pickle.load(open(SC+'par_U.pkl','rb'))
byslot={}
for i,u in enumerate(U): byslot.setdefault(u['slot'],[]).append(i)
PAIRS=[('A1-F','A1-R'),('A2-F','A2-R'),('A3-F','A3-R'),('A3-F','V5-R'),('A4-F','A4-R'),('A5-F','A5-R')]
WINS=[(25,245),(357,515),(577,784),(788,906),(986,1176),(1242,1491)]
n=len(xh); K=len(WINS)
# 预计算每个窗口的匹配矩阵 S_w（假设两条序列都有该窗口）
t=time.time(); Sw=[]; lens=np.array([b-a+1 for a,b in WINS],np.float32)
for a,b in WINS:
    cols=range(a-1,b); X=np.zeros((n,len(cols)*4),np.float32)
    for k,p in enumerate(cols):
        m=xh[:,p]<4; X[np.where(m)[0],k*4+xh[m,p]]=1
    Sw.append(X@X.T)
print('S_w',round(time.time()-t,1),'s',flush=True)
same=(genh[None,:]==genh[:,None]); q=pd.Series(genh).map(pd.Series(genh).value_counts()).values>=2
RULES=[('abs60',60,0.0),('abs300',300,0.0),('frac0.8',0,0.8)]
def amp_masks(sel):
    """sel：保留的寡核苷酸下标集合 -> 各窗口的 held 命中；以及 GG/SILVA 的命中"""
    def OR(slot,key):
        idx=[i for i in byslot.get(slot,[]) if i in sel]
        if not idx: return None
        if key=='held': return np.any([U[i]['held'] for i in idx],axis=0)
        return None
    P=[]
    for f,r in PAIRS:
        a=OR(f,'held'); b=OR(r,'held'); P.append(np.zeros(n,bool) if a is None or b is None else (a&b))
    return np.array(P)
def acc_rules(P):
    Pm=P.T.astype(np.float32)                    # n×K
    Lm=(Pm*lens)@Pm.T
    Sm=np.zeros((n,n),np.float32)
    for w in range(K): Sm+=Sw[w]*np.outer(Pm[:,w],Pm[:,w])
    ident=P.any(0); Lq=(Pm*lens).sum(1,keepdims=True)
    out={}
    for nm,ma,mf in RULES:
        need=np.maximum(ma,mf*Lq)
        rate=np.where(Lm>=need,(Lm-Sm)/np.maximum(Lm,1),np.inf); np.fill_diagonal(rate,np.inf)
        mn=rate.min(1,keepdims=True); tied=(rate==mn)&np.isfinite(rate)
        credit=np.where(tied.any(1),(tied&same).sum(1)/np.maximum(tied.sum(1),1),0.0)*ident
        out[nm]=float(credit[q].mean())
    return out
def acc_full(P):
    Pm=P.T.astype(np.float32); Lq=(Pm*lens).sum(1)
    Sq=np.zeros((n,n),np.float32)
    for w in range(K): Sq+=Pm[:,w][:,None]*Sw[w]
    rate=np.where(Lq[:,None]>0,(Lq[:,None]-Sq)/np.maximum(Lq[:,None],1),np.inf); np.fill_diagonal(rate,np.inf)
    mn=rate.min(1,keepdims=True); tied=(rate==mn)&np.isfinite(rate)
    credit=np.where(tied.any(1),(tied&same).sum(1)/np.maximum(tied.sum(1),1),0.0)*(Lq>0)
    return float(credit[q].mean())
def cov_metrics(sel):
    r={}
    for vn,(ref,rows_) in VAL.items():
        phs=ref.phylum[rows_]; M=[]
        for f,rr in PAIRS:
            a=[i for i in byslot.get(f,[]) if i in sel]; b=[i for i in byslot.get(rr,[]) if i in sel]
            M.append(np.zeros(len(rows_),bool) if not a or not b else (np.any([U[i]['val'][vn] for i in a],0)&np.any([U[i]['val'][vn] for i in b],0)))
        M=np.array(M).T
        per=[M[phs==p].mean(0) for p in pdz.KEY_PHYLA if (phs==p).sum()>=30]
        r[vn+'_amp_mean']=float(np.mean(per)); r[vn+'_all']=float(M.all(1).mean())
    for db,(ref,rr) in refs.items():
        M=[]
        for f,rr_ in PAIRS:
            a=[i for i in byslot.get(f,[]) if i in sel]; b=[i for i in byslot.get(rr_,[]) if i in sel]
            if not a or not b: M.append(np.zeros(len(U[0]['ten'][db][0]),bool)); continue
            ha=np.any([U[i]['ten'][db][0] for i in a],0); hb=np.any([U[i]['ten'][db][0] for i in b],0); M.append(ha&hb)
        M=np.array(M).T; r['Ten_'+db+'_ge3']=float((M.sum(1)>=3).mean())
    return r
def evaluate(sel):
    P=amp_masks(sel); a=acc_rules(P); r=dict(oligos=len(sel),expansions=int(sum(U[i]['nexp'] for i in sel)),**a); r['acc3']=float(np.mean(list(a.values()))); r['accf']=acc_full(P); r.update(cov_metrics(sel)); r['cov']=(r['GG_amp_mean']+r['SILVA_amp_mean'])/2; r['ten']=(r['Ten_GG_ge3']+r['Ten_SILVA_ge3'])/2; r['J']=(r['accf']+r['cov']+r['ten'])/3; return r
if __name__=='__main__':
    base=[i for i,u in enumerate(U) if 'SNAP' not in u['name']]
    t=time.time(); r=evaluate(set(base)); print('基线（22 条，无 SNAP 替代）',{k:round(v,4) for k,v in r.items()},round(time.time()-t,1),'s')
