# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/a5f_*
import pickle, json, itertools, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore'); pd.set_option('display.width',250); pd.set_option('display.max_colwidth',60)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'ext2.py').read().split("rows=[]")[0]
exec(src)
sites=json.load(open(SC+'a2r_final5.json')); s=sites['A5-F']
F1='GCTACACRCRTGCTACAAT'; F2='GCTACACACGTVMTACAAT'
IU={'A':'A','C':'C','G':'G','T':'T','R':'AG','Y':'CT','S':'GC','W':'AT','K':'GT','M':'AC','B':'CGT','D':'AGT','H':'ACT','V':'ACG','N':'ACGT'}
code={frozenset(v):k for k,v in IU.items()}
def broaden(p):
    out={p}
    for i,c in enumerate(p):
        for extra in 'ACGT':
            if extra in IU[c]: continue
            new=frozenset(IU[c]+extra)
            if new in code: out.add(p[:i]+code[new]+p[i+1:])
    return out
rows=[]
for f1 in broaden(F1):
    for f2 in broaden(F2):
        if len(pdz.expand(f1))+len(pdz.expand(f2))>16: continue
        site=dict(s,prim=[f1,f2],ext=[3,3]); sev,hp=cp.pool_issues({'A5-F':[f1,f2]}); kc=keycov(site)
        rows.append(dict(f1=f1,f2=f2,Tm1=round(cp.oligo_tm(f1)[0],1),Tm2=round(cp.oligo_tm(f2)[0],1),nexp=len(pdz.expand(f1))+len(pdz.expand(f2)),SILVA=kc['SILVA'],GG=kc['GG'],severe=len(sev),worst=min([x[4] for x in sev],default=0)))
R=pd.DataFrame(rows); R['cov']=(R.SILVA+R.GG)/2
g=R[(R.severe==0)&(R.Tm1>=55)&(R.Tm2>=55)].sort_values('cov',ascending=False).head(15); print(g.to_string()); R.to_csv(SC+'a5f_fix2.csv',index=False)
print('原',keycov(s))
