# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/v3_*
import pickle, json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore'); pd.set_option('display.width',250); pd.set_option('display.max_colwidth',44)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'ext2.py').read().split("rows=[]")[0]
exec(src)
D1=json.load(open(SC+'short_sites.json'))['D1']
ints0=[(D1[f'A{k}-F']['start']+D1[f'A{k}-F']['L'],D1[f'A{k}-R']['start']-1) for k in range(1,6)]
base=ts.accuracy_direct(xall,ints0,gall,vall)
D=pd.read_pickle(SC+'v3_scan.pkl'); D['kind']=D.kind.str.replace('A4','A2')
ok=D[(D.severe==0)&(D.cov_min>=0.90)&(D.tmmin>=54)&(D.tmmax<=67)]
best=ok.sort_values(['cov_min'],ascending=False).groupby(['kind','pos','L']).head(1)
print(best.groupby('kind').size())
F=best[best.kind=='A2F']; R=best[best.kind=='A2R']
print(F.sort_values('cov_min',ascending=False).head(10)[['pos','L','k','cov_min','nexp','tm','prim']].to_string())
out=[]
for _,f in F.iterrows():
    for _,r in R.iterrows():
        f5=f.pos-f.L+1; r5=r.pos+r.L-1; alen=r5-f5+1
        if f5<286 or r5>535 or alen>270: continue
        ints=list(ints0); ints[1]=(f.pos+1,r.pos-1); acc=ts.accuracy_direct(xall,ints,gall,vall)
        out.append(dict(F_end=f.pos,F_L=f.L,F_prim=' | '.join(f.prim),F_tm=f.tm,F_cov=round(f.cov_min,3),F_nexp=f.nexp,R_start=r.pos,R_L=r.L,R_prim=' | '.join(r.prim),R_tm=r.tm,R_cov=round(r.cov_min,3),R_nexp=r.nexp,f5=f5,r5=r5,len=alen,gain=round((acc-base)*100,3)))
O=pd.DataFrame(out); O['cov']=O[['F_cov','R_cov']].min(axis=1); O.to_pickle(SC+'v3_pairs.pkl'); print(len(O))
for lo,hi in ((187,215),(216,240),(241,270)):
    g=O[(O.len>=lo)&(O.len<=hi)].sort_values(['gain','cov'],ascending=False).head(6); print(lo,hi); print(g[['F_end','F_L','F_prim','F_cov','R_start','R_L','R_prim','R_cov','f5','r5','len','gain']].to_string())
