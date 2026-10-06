# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/layout6_pools_20_28_36.csv
import sys; sys.path.insert(0,'/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad')
from sc2 import *
newR=[('16S-A5-R.1','ATGMTGAYTTGACGTCRTC'),('16S-A5-R.2','ATRCTGACYTGACGTCRTC')]
ids5=[add(n,'L6A5-R','R',1188,s) for n,s in newR]
byslot['L6A5-R']=ids5
W=list(WINS); W[4]=(986,1187); st=STRUCT['L6']; st['wins']=W; st['S']=build_S(W); st['lens']=np.array([b-a+1 for a,b in W],np.float32)
base=({L6[n] for n in T.name if not n.startswith('16S-A5-R') and not n.startswith('16S-A1-F')})|{L6['V1_f']}|set(ids5)
cur=set(base); vis={}
def ev(sel):
    k=frozenset(sel)
    if k not in vis: vis[k]=evaluate(sel,'L6')
    return vis[k]
r=ev(cur); print('起点',len(cur),{k:round(v,4) for k,v in r.items() if isinstance(v,float)},flush=True)
pools={len(cur):set(cur)}; hist=[]
while len(cur)>16:
    best=None
    for i in sorted(cur):
        loss=ev(cur)['J']-ev(cur-{i})['J']; tie=-U[i]['nexp']
        if best is None or (loss,tie)<(best[0],best[1]): best=(loss,tie,i)
    cur=cur-{best[2]}; r=ev(cur); hist.append((len(cur),U[best[2]]['name'],r)); pools[len(cur)]=set(cur)
    if len(cur) in (44,40,36,32,28,24,20,16): print(len(cur),'删',U[best[2]]['name'],'exp',r['expansions'],'J %.4f accf %.4f cov %.3f ten %.3f'%(r['J'],r['accf'],r['cov'],r['ten']),flush=True)
pickle.dump((hist,pools,{i:(U[i]['name'],U[i]['seq'],U[i]['s']['start'],U[i]['nexp']) for i in range(base_n) } if False else None),open(SC+'sc4_hist.pkl','wb'))
pickle.dump(dict(hist=hist,pools={k:sorted(v) for k,v in pools.items()},info={i:(U[i]['name'],U[i]['seq'],U[i]['s']['start'],U[i]['s']['orient'],U[i]['nexp']) for i in range(len(U)) if U[i]['slot'].startswith('L6')},res={k:ev(v) for k,v in pools.items()}),open(SC+'sc4_out.pkl','wb'))
print('done')
