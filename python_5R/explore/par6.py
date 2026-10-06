# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/par_*
import sys; sys.path.insert(0,'/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad')
from par2 import *
import compact_primers as cp
visited,paths=pickle.load(open(SC+'par_visited4.pkl','rb'))
def ev(sel):
    key=frozenset(sel)
    if key not in visited: visited[key]=evaluate(set(sel))
    return visited[key]
base=set(i for i,u in enumerate(U) if 'SNAP' not in u['name'])
# 在 22 条基线（不含借鉴）上做向后消除
cur=set(base); order=[]
while len(cur)>1:
    best=None
    for i in sorted(cur):
        loss=ev(cur)['J']-ev(cur-{i})['J']
        if best is None or (loss,-U[i]['nexp'])<(best[0],best[1]): best=(loss,-U[i]['nexp'],i)
    cur=cur-{best[2]}; order.append((best[2],len(cur),ev(cur)['J'],ev(cur)['expansions']))
rows=[dict(step=k+1,removed=U[i]['name'],slot=U[i]['slot'],seq=U[i]['seq'],nexp=U[i]['nexp'],left=left,J=round(J,4),expansions=e) for k,(i,left,J,e) in enumerate(order)]
O=pd.DataFrame(rows); print(O.to_string()); O.to_csv(SC+'par_base_order.csv',index=False)
# 帕累托前沿（全部状态）
D=pd.DataFrame([dict(sel=sorted(k),**r) for k,r in visited.items()])
def pareto(D,cost):
    D=D.sort_values([cost,'J'],ascending=[True,False]); b=-1; keep=[]
    for i,r in D.iterrows():
        if r['J']>b+1e-12: keep.append(i); b=r['J']
    return D.loc[keep]
P1=pareto(D,'oligos'); P2=pareto(D,'expansions')
names=lambda sel:[U[i]['name'] for i in sel]
for nm,P in (('oligos',P1),('exp',P2)):
    P=P.copy(); P['names']=P['sel'].apply(lambda s:' | '.join(names(s))); P['idx']=P['sel'].apply(lambda s:','.join(map(str,s)))
    P.drop(columns=['sel']).to_csv(SC+f'par_frontier_{nm}.csv',index=False)
# 代表性水平：二聚体、Tm
def stats(sel):
    pool={}
    for i in sel: pool.setdefault(U[i]['slot'],[]).append(U[i]['seq'])
    sev,hp=cp.pool_issues(pool)
    tms=[cp.oligo_tm(U[i]['seq'])[0] for i in sel]
    return len(sev),len(hp),round(min(tms),1),round(max(tms),1)
lv=[]
for k in (24,19,17,16,14,12,11):
    r=P1[P1.oligos==k]
    if len(r)==0: continue
    r=r.iloc[0]; s=stats(r['sel']); lv.append(dict(level=k,**{kk:r[kk] for kk in ('oligos','expansions','J','accf','cov','ten','GG_all','SILVA_all','GG_amp_mean','SILVA_amp_mean','Ten_GG_ge3','Ten_SILVA_ge3')},severe=s[0],hairpins=s[1],Tm_min=s[2],Tm_max=s[3],sel=','.join(map(str,r['sel']))))
L=pd.DataFrame(lv); print(L.drop(columns=['sel']).round(4).to_string()); L.to_csv(SC+'par_levels.csv',index=False)
pickle.dump(([dict(i=i,name=u['name'],slot=u['slot'],seq=u['seq'],nexp=u['nexp']) for i,u in enumerate(U)],),open(SC+'par_U_meta.pkl','wb'))
pickle.dump((visited,paths),open(SC+'par_visited5.pkl','wb'))
