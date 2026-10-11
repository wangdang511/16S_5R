# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/a5f_*
import pickle, json, itertools, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore'); pd.set_option('display.width',250); pd.set_option('display.max_colwidth',60)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'ext2.py').read().split("rows=[]")[0]
exec(src)
sites=json.load(open(SC+'a2r_final5.json')); s=sites['A5-F']; print(s)
p0=s['prim'][0]
exps=pdz.expand(p0)
for a in exps:
    h,e=cp.dimer(a,a); print(a,round(h,1),round(e,1))
IU={'R':'AG','Y':'CT','S':'GC','W':'AT','K':'GT','M':'AC','B':'CGT','D':'AGT','H':'ACT','V':'ACG','N':'ACGT'}
sub=lambda c:[c2 for c2 in 'ACGT'+'RYSWKM' if set(IU.get(c2,c2))<set(IU[c])]
vars_={p0:0}
for i,c in enumerate(p0):
    if c in IU:
        for c2 in sub(c):
            q=p0[:i]+c2+p0[i+1:]; vars_[q]=1
for i,j in itertools.combinations([i for i,c in enumerate(p0) if c in IU],2):
    for ci in sub(p0[i]):
        for cj in sub(p0[j]):
            q=list(p0); q[i]=ci; q[j]=cj; vars_[''.join(q)]=2
rows=[]
for q,n in vars_.items():
    prim=list(s['prim']); prim[0]=q; site=dict(s,prim=prim)
    sev,hp=cp.pool_issues({'A5-F':prim}); kc=keycov(site)
    rows.append(dict(seq=q,changes=n,Tm=round(cp.oligo_tm(q)[0],1),nexp=len(pdz.expand(q)),SILVA=kc['SILVA'],GG=kc['GG'],severe=len(sev),worst=min([x[4] for x in sev],default=0)))
R=pd.DataFrame(rows); R['cov']=(R.SILVA+R.GG)/2
print(R.sort_values('cov',ascending=False).head(25).to_string()); R.to_csv(SC+'a5f_fix.csv',index=False)
