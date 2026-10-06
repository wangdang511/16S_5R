# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/layout6_pareto.csv
import sys; sys.path.insert(0,'/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad')
from sc1 import *
import pandas as pd
T=pd.read_csv(SC+'lay3_tm.csv')
base=len(U)
def add(name,slot,orient,start,seq):
    s=dict(orient=orient,start=int(start),L=len(seq),prim=[seq],ext=[0])
    u=dict(name=name,slot=slot,s=s,seq=seq,nexp=len(pdz.expand(seq)))
    u['held']=site_hit(Gh,np.arange(len(held)),s); u['val']={vn:site_hit(ref,rows_,s) for vn,(ref,rows_) in VAL.items()}
    u['ten']={db:site_cov_ten(s,ref,rr) for db,(ref,rr) in refs.items()}
    U.append(u); return len(U)-1
L6={}
for r in T.itertuples():
    L6[r.name]=add(r.name,'L6'+r.site,r.site[-1],r.start,r.seq_new)
s=SNAP['V1_f']; L6['V1_f']=add('V1_f(SNAP 9–27)','L6A1-F','F',s['start'],s['prim'][0])
for i,x in enumerate(U): byslot.setdefault(x['slot'],[]).append(i)
WINS=[(27,245),(333,503),(537,773),(807,953),(986,1173),(1241,1485)]
STRUCT['L6']=dict(excl=set(range(base)),wins=WINS,wm=[[('L6A%d-F'%k,'L6A%d-R'%k)] for k in range(1,7)],pairs=[('L6A%d-F'%k,'L6A%d-R'%k) for k in range(1,7)])
st=STRUCT['L6']; st['S']=build_S(st['wins']); st['lens']=np.array([b-a+1 for a,b in st['wins']],np.float32)
if __name__=='__main__':
    allL6=set(L6.values()); noV1=allL6-{L6['V1_f']}; mainonly={L6[n] for n in T[T.kind=='main'].name}
    withV1=(noV1-{L6[n] for n in T[T.name.str.startswith('16S-A1-F')].name})|{L6['V1_f']}
    for nm,sel in (('51 条（含补充）',noV1),('仅主引物 37 条',mainonly),('51 条 A1-F 换 V1_f',withV1)):
        r=evaluate(sel,'L6'); print(nm,len(sel),{k:round(v,4) if isinstance(v,float) else v for k,v in r.items()},flush=True)
    # 逐条后向删除（从 withV1 起）
    cur=set(withV1); vis={}
    def ev(sel):
        k=frozenset(sel)
        if k not in vis: vis[k]=evaluate(sel,'L6')
        return vis[k]
    hist=[]
    while len(cur)>12:
        best=None
        for i in sorted(cur):
            loss=ev(cur)['J']-ev(cur-{i})['J']; tie=-U[i]['nexp']
            if best is None or (loss,tie)<(best[0],best[1]): best=(loss,tie,i)
        cur=cur-{best[2]}; r=ev(cur); hist.append((len(cur),U[best[2]]['name'],r))
        if len(cur) in (40,36,32,28,24,22,20,18,16,14,12): print(len(cur),'删',U[best[2]]['name'],'exp',r['expansions'],'J %.4f accf %.4f cov %.3f ten %.3f'%(r['J'],r['accf'],r['cov'],r['ten']),flush=True)
    pickle.dump((hist,{k:v for k,v in vis.items()},L6),open(SC+'sc2_hist.pkl','wb'))
    json.dump({str(len(h[2])):1 for h in hist},open(SC+'sc2_ok.json','w'))
