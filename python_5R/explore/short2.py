# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/short_*
import pickle, json, numpy as np, pandas as pd, warnings, itertools
warnings.filterwarnings('ignore'); pd.set_option('display.width',250); pd.set_option('display.max_colwidth',50)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'ext2.py').read().split("rows=[]")[0]
exec(src)
sites=json.load(open(SC+'a5f_final5.json')); ints0=[(sites[f'A{k}-F']['start']+sites[f'A{k}-F']['L'],sites[f'A{k}-R']['start']-1) for k in range(1,6)]
base=ts.accuracy_direct(xall,ints0,gall,vall)
D=pd.read_pickle(SC+'short_scan.pkl')
ok=D[(D.severe==0)&(D.cov_min>=0.93)&(D.tmmin>=54)&(D.tmmax<=67)]
# 每个 (kind,pos,L) 取覆盖率最好、其次展开数最少
best=ok.sort_values(['cov_min'],ascending=False).groupby(['kind','pos','L']).head(1)
print(best.groupby('kind').size())
out=[]
for amp,idx in (('A4',3),('A5',4)):
    F=best[best.kind==amp+'F']; R=best[best.kind==amp+'R']
    for _,f in F.iterrows():
        for _,r in R.iterrows():
            s=f.pos+1; e=r.pos-1; alen=(r.pos+r.L-1)-(f.pos-f.L+1)+1
            if alen>262 or e-s<60: continue
            ints=list(ints0); ints[idx]=(s,e)
            acc=ts.accuracy_direct(xall,ints,gall,vall)
            out.append(dict(amp=amp,F_end=f.pos,F_L=f.L,F_prim=' | '.join(f.prim),F_tm=f.tm,F_cov=round(f.cov_min,3),F_nexp=f.nexp,
                            R_start=r.pos,R_L=r.L,R_prim=' | '.join(r.prim),R_tm=r.tm,R_cov=round(r.cov_min,3),R_nexp=r.nexp,len=alen,loss=round((base-acc)*100,3)))
O=pd.DataFrame(out); O['cov']=O[['F_cov','R_cov']].min(axis=1); O.to_pickle(SC+'short_pairs.pkl'); print(len(O))
for amp in ('A4','A5'):
    g=O[(O.amp==amp)&(O.loss<=0.12)].sort_values(['len','cov'],ascending=[True,False]); print(amp); print(g.head(14)[['F_end','F_L','F_prim','F_cov','R_start','R_L','R_prim','R_cov','len','loss']].to_string())
