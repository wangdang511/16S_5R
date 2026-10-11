# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/layout6_pareto.csv
import pickle, json, itertools, time, numpy as np, pandas as pd, warnings, sys
warnings.filterwarnings('ignore'); pd.set_option('display.width',250); pd.set_option('display.max_colwidth',50)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'ten1.py').read().split("rows=[]; masks={}")[0]
exec(src)
U=pickle.load(open(SC+'par_U.pkl','rb')); SNAP=pickle.load(open(SC+'snap_sites.pkl','rb'))
# 加入 515F（SNAP V4_f）
s=SNAP['V4_f']; u=dict(name='515F（SNAP V4_f）',slot='A3-F',s=dict(orient=s['orient'],start=s['start'],L=s['L'],prim=[s['prim'][0]],ext=[0]),seq=s['prim'][0],nexp=len(pdz.expand(s['prim'][0])))
rowsH=np.arange(len(held)); u['held']=site_hit(Gh,rowsH,u['s']); u['val']={vn:site_hit(ref,rows_,u['s']) for vn,(ref,rows_) in VAL.items()}; u['ten']={db:site_cov_ten(u['s'],ref,rr) for db,(ref,rr) in refs.items()}
U.append(u); I515=len(U)-1
byslot={}
for i,x in enumerate(U): byslot.setdefault(x['slot'],[]).append(i)
n=len(xh); same=(genh[None,:]==genh[:,None]); q=pd.Series(genh).map(pd.Series(genh).value_counts()).values>=2
def build_S(WINS):
    out=[]
    for a,b in WINS:
        cols=range(a-1,b); X=np.zeros((n,len(cols)*4),np.float32)
        for k,p in enumerate(cols):
            m=xh[:,p]<4; X[np.where(m)[0],k*4+xh[m,p]]=1
        out.append(X@X.T)
    return out
# 结构：窗口 + 每个窗口的 (F 位点, R 位点) 配对集合；扩增子对列表
STRUCT={
 'B':dict(excl={I515},wins=[(25,245),(357,515),(577,784),(788,906),(986,1176),(1242,1491)],
          wm=[[('A1-F','A1-R')],[('A2-F','A2-R')],[('A3-F','A3-R')],[('A3-F','A4-R')],[('A5-F','A5-R')],[('A6-F','A6-R')]],
          pairs=[('A1-F','A1-R'),('A2-F','A2-R'),('A3-F','A3-R'),('A3-F','A4-R'),('A5-F','A5-R'),('A6-F','A6-R')]),
 'C':dict(excl=set(byslot['A2-R']),wins=[(25,245),(357,477),(534,653),(665,784),(788,906),(986,1176),(1242,1491)],
          wm=[[('A1-F','A1-R')],[('A2-F','A3-R')],[('A3-F','A3-R'),('A3-F','A4-R')],[('A2-F','A3-R'),('A3-F','A3-R')],[('A3-F','A4-R')],[('A5-F','A5-R')],[('A6-F','A6-R')]],
          pairs=[('A1-F','A1-R'),('A2-F','A3-R'),('A3-F','A3-R'),('A3-F','A4-R'),('A5-F','A5-R'),('A6-F','A6-R')])}
for k,st in STRUCT.items(): st['S']=build_S(st['wins']); st['lens']=np.array([b-a+1 for a,b in st['wins']],np.float32)
def slot_hit(sel,slot,key,db=None):
    idx=[i for i in byslot.get(slot,[]) if i in sel]
    if not idx: return None
    if key=='held': return np.any([U[i]['held'] for i in idx],0)
    if key=='val': return np.any([U[i]['val'][db] for i in idx],0)
    if key=='ten': return np.any([U[i]['ten'][db][0] for i in idx],0)
def pair_hit(sel,f,r,key,db=None):
    a=slot_hit(sel,f,key,db); b=slot_hit(sel,r,key,db)
    if a is None or b is None:
        ln={'held':n}.get(key,None) or (len(U[0]['val'][db]) if key=='val' else len(U[0]['ten'][db][0]))
        return np.zeros(ln,bool)
    return a&b
def acc_full(P,st):
    Pm=P.T.astype(np.float32); Lq=(Pm*st['lens']).sum(1); Sq=np.zeros((n,n),np.float32)
    for w in range(len(st['wins'])): Sq+=Pm[:,w][:,None]*st['S'][w]
    rate=np.where(Lq[:,None]>0,(Lq[:,None]-Sq)/np.maximum(Lq[:,None],1),np.inf); np.fill_diagonal(rate,np.inf)
    mn=rate.min(1,keepdims=True); tied=(rate==mn)&np.isfinite(rate)
    credit=np.where(tied.any(1),(tied&same).sum(1)/np.maximum(tied.sum(1),1),0.0)*(Lq>0)
    return float(credit[q].mean())
def evaluate(sel,struct):
    st=STRUCT[struct]; sel=set(sel)
    P=np.array([np.any([pair_hit(sel,f,r,'held') for f,r in pairs],0) for pairs in st['wm']])
    r=dict(oligos=len(sel),expansions=int(sum(U[i]['nexp'] for i in sel)),accf=acc_full(P,st))
    cov=[]
    for vn,(ref,rows_) in VAL.items():
        phs=ref.phylum[rows_]; M=np.array([pair_hit(sel,f,rr,'val',vn) for f,rr in st['pairs']]).T
        per=[M[phs==p].mean(0) for p in pdz.KEY_PHYLA if (phs==p).sum()>=30]; cov.append(float(np.mean(per))); r[vn+'_all']=float(M.all(1).mean())
    tens=[]
    for db in refs:
        M=np.array([pair_hit(sel,f,rr,'ten',db) for f,rr in st['pairs']]).T; tens.append(float((M.sum(1)>=3).mean()))
    r['cov']=float(np.mean(cov)); r['ten']=float(np.mean(tens)); r['J']=(r['accf']+r['cov']+r['ten'])/3; return r
if __name__=='__main__':
    visited={}
    def ev(sel,st):
        key=(st,frozenset(sel))
        if key not in visited: visited[key]=evaluate(sel,st)
        return visited[key]
    for stn in ('B','C'):
        cur=set(range(len(U)))-STRUCT[stn]['excl']; r=ev(cur,stn); print(stn,'全部',len(cur),{k:round(v,4) for k,v in r.items()},flush=True)
        while len(cur)>8:
            best=None
            for i in sorted(cur):
                loss=ev(cur,stn)['J']-ev(cur-{i},stn)['J']; tie=-U[i]['nexp']
                if best is None or (loss,tie)<(best[0],best[1]): best=(loss,tie,i)
            cur=cur-{best[2]}
            if len(cur) in (24,22,20,18,16,15,14,13,12,11,10): print(stn,len(cur),'删',U[best[2]]['name'],'J %.4f accf %.4f cov %.3f ten %.3f'%tuple(ev(cur,stn)[k] for k in ('J','accf','cov','ten')),flush=True)
    pickle.dump(visited,open(SC+'sc_visited.pkl','wb')); pickle.dump(U,open(SC+'sc_U.pkl','wb'))
