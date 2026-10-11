# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/extend_*
import pickle, json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore'); pd.set_option('display.width',250)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'v2_tm.py').read().split("print('pool Tm")[0]
exec(src)
def keycov(s):
    out={}
    for vn,(ref,rows_) in VAL.items():
        phs=ref.phylum[rows_]; h=site_hit(ref,rows_,s)
        out[vn]=round(float(np.mean([h[phs==g].mean() for g in pdz.KEY_PHYLA if (phs==g).sum()>=30])),3)
    return out
cand,pool=pickle.load(open(SC+'ext_cand.pkl','rb'))
rows=[]
for k,v in cand.items():
    s=pool[v['site']]; prim=list(s['prim']); ext=list(s['ext']); prim[v['i']]=v['seq']; ext[v['i']]=v['e']
    site=dict(s,prim=prim,ext=ext); kc=keycov(site)
    rows.append(dict(name=k,site=v['site'],i=v['i'],e=v['e'],seq=v['seq'],nt=len(v['seq']),Tm=round(v['tm'],1),SILVA=kc['SILVA'],GG=kc['GG']))
df=pd.DataFrame(rows); df.to_csv(SC+'ext_cov.csv',index=False); print(df.to_string())
