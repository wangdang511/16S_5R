# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/par_*
import sys; sys.path.insert(0,'/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad')
from par2 import *
visited,paths=pickle.load(open(SC+'par_visited3.pkl','rb'))
def ev(sel):
    key=frozenset(sel)
    if key not in visited: visited[key]=evaluate(set(sel))
    return visited[key]
# 局部搜索：同位点互换；在每个规模上从已有最好状态出发
best_by_size={}
for key,r in visited.items():
    k=len(key)
    if k not in best_by_size or r['J']>best_by_size[k][0]: best_by_size[k]=(r['J'],key)
def climb(sel):
    cur=set(sel); improved=True
    while improved:
        improved=False
        for i in sorted(cur):
            for j in byslot[U[i]['slot']]:
                if j in cur: continue
                cand=(cur-{i})|{j}
                if ev(cand)['J']>ev(cur)['J']+1e-9: cur=cand; improved=True; break
            if improved: break
    return frozenset(cur)
for k in range(11,25):
    if k in best_by_size:
        s=climb(best_by_size[k][1]); print(k,'J %.4f'%best_by_size[k][0],'->','%.4f'%ev(s)['J'],flush=True)
# 明确的"借鉴"组合
base=set(i for i,u in enumerate(U) if 'SNAP' not in u['name'])
v1=[i for i,u in enumerate(U) if u['name'].startswith('V1_f')][0]; a1f=[i for i,u in enumerate(U) if u['name']=='A1-F'][0]
v5s=[i for i,u in enumerate(U) if u['name'].startswith('V5_r')][0]; v5a=[i for i,u in enumerate(U) if u['slot']=='V5-R' and 'SNAP' not in u['name']]
tests={'基线 22 条':base,'A1-F 换成 SNAP V1_f':(base-{a1f})|{v1},'V5-R 换成 SNAP 907R（1 条）':(base-set(v5a))|{v5s},'两者都换':(base-{a1f}-set(v5a))|{v1,v5s}}
for nm,sel in tests.items():
    r=ev(sel); print(nm,len(sel),r['expansions'],'J %.4f accf %.4f cov %.3f ten %.3f GG_all %.3f SIL_all %.3f'%(r['J'],r['accf'],r['cov'],r['ten'],r['GG_all'],r['SILVA_all']))
pickle.dump((visited,paths),open(SC+'par_visited4.pkl','wb'))
