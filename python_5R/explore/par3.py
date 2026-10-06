# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/par_*
import sys; sys.path.insert(0,'/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad')
from par2 import *
visited={}
def ev(sel):
    key=frozenset(sel)
    if key not in visited: visited[key]=evaluate(set(sel))
    return visited[key]
allidx=set(range(len(U)))
paths={}
def backward(mode):
    cur=set(allidx); path=[(frozenset(cur),None)]
    while len(cur)>1:
        best=None
        for i in sorted(cur):
            r=ev(cur-{i}); base=ev(cur)
            loss=base['J']-r['J']
            score=loss if mode=='count' else loss/max(U[i]['nexp'],1)
            tie=-U[i]['nexp']
            if best is None or (score,tie)<(best[0],best[1]): best=(score,tie,i)
        cur=cur-{best[2]}; path.append((frozenset(cur),best[2]))
        print(mode,len(cur),'移除',U[best[2]]['name'],'J %.4f'%ev(cur)['J'],flush=True)
    return path
paths['count']=backward('count'); paths['exp']=backward('exp')
# 前向：每个位点先选 1 条最好的，再逐条加
def forward():
    slots=list(byslot)
    cur=set()
    for s in slots: cur.add(byslot[s][0])
    improved=True
    while improved:
        improved=False
        for s in slots:
            best=(ev(cur)['J'],None)
            for i in byslot[s]:
                cand=(cur-set(byslot[s]))|{i}; r=ev(cand)['J']
                if r>best[0]+1e-9: best=(r,cand)
            if best[1] is not None: cur=best[1]; improved=True
    path=[(frozenset(cur),None)]
    while True:
        best=None
        for i in allidx-cur:
            r=ev(cur|{i})['J']
            if best is None or r>best[0]: best=(r,i)
        if best is None: break
        cur=cur|{best[1]}; path.append((frozenset(cur),best[1])); print('forward',len(cur),'加',U[best[1]]['name'],'J %.4f'%best[0],flush=True)
    return path
paths['forward']=forward()
pickle.dump((visited,paths),open(SC+'par_visited3.pkl','wb')); print('visited',len(visited))
