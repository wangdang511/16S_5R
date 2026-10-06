# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/par_*
import sys; sys.path.insert(0,'/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad')
from par2 import *
visited,paths=pickle.load(open(SC+'par_visited3.pkl','rb'))
rows=[]
for key,r in visited.items():
    rows.append(dict(sel=sorted(key),**r))
D=pd.DataFrame(rows); print(len(D),'个状态')
def pareto(D,cost):
    # 最大化 J，最小化 cost
    D=D.sort_values([cost,'J'],ascending=[True,False]); best=-1; keep=[]
    for i,r in D.iterrows():
        if r['J']>best+1e-12: keep.append(i); best=r['J']
    return D.loc[keep]
P1=pareto(D,'oligos'); P2=pareto(D,'expansions')
names=lambda sel:[U[i]['name'] for i in sel]
for nm,P in (('按寡核苷酸数',P1),('按展开序列数',P2)):
    print('==',nm)
    for _,r in P.iterrows(): print(int(r.oligos),int(r.expansions),'J %.4f accf %.4f cov %.3f ten %.3f GG_all %.3f SIL_all %.3f'%(r.J,r.accf,r['cov'],r['ten'],r.GG_all,r.SILVA_all))
# 单条留一（相对完整 24 条与相对 22 条基线）
full=frozenset(range(len(U))); base=frozenset(i for i,u in enumerate(U) if 'SNAP' not in u['name'])
loo=[]
for i in sorted(base):
    r=evaluate(set(base-{i})); b=evaluate(set(base))
    loo.append(dict(name=U[i]['name'],slot=U[i]['slot'],seq=U[i]['seq'],nexp=U[i]['nexp'],dJ=round(b['J']-r['J'],4),daccf=round(b['accf']-r['accf'],4),dcov=round(b['cov']-r['cov'],4),dten=round(b['ten']-r['ten'],4),dall_GG=round(b['GG_all']-r['GG_all'],4),dall_SILVA=round(b['SILVA_all']-r['SILVA_all'],4)))
L=pd.DataFrame(loo).sort_values('dJ'); print(L.to_string()); L.to_csv(SC+'par_loo.csv',index=False)
D['names']=D['sel'].apply(lambda s:' | '.join(names(s)))
D.drop(columns=['sel']).to_csv(SC+'par_all_states.csv',index=False)
for nm,P in (('oligos',P1),('exp',P2)):
    P=P.copy(); P['names']=P['sel'].apply(lambda s:' | '.join(names(s))); P['sel_idx']=P['sel'].apply(lambda s:','.join(map(str,s))); P.drop(columns=['sel']).to_csv(SC+f'par_frontier_{nm}.csv',index=False)
pickle.dump((D,P1,P2,L),open(SC+'par_frontier.pkl','wb'))
