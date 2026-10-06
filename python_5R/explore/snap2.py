# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/snap_*
import pickle, json, numpy as np, pandas as pd, warnings, collections
warnings.filterwarnings('ignore'); pd.set_option('display.width',250); pd.set_option('display.max_colwidth',60)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'ten1.py').read().split("rows=[]; masks={}")[0]
exec(src)
P=pd.read_csv(SC+'snap_positions.csv')
groups=collections.OrderedDict()
for _,r in P.iterrows():
    g=re.sub(r'_(\d)_r$','_r',r['name']) if False else r['name']
import re
gname=lambda n: re.sub(r'_\d_r$','_r',re.sub(r'_f\d$','_f',n))
for _,r in P.iterrows(): groups.setdefault(gname(r['name']),[]).append(r)
SITES={}
for g,rs in groups.items():
    o='F' if g.endswith('_f') else 'R'
    r0=rs[0]; SITES[g]=dict(orient=o,start=int(r0.start),L=int(r0.nt) if o=='F' else int(r0.nt),prim=[x.seq for x in rs],ext=[0]*len(rs))
# 覆盖：Greengenes/SILVA 验证集（主要门），以及支原体
rows=[]
for g,s in SITES.items():
    kc={}
    for vn,(ref,rows_) in VAL.items():
        phs=ref.phylum[rows_]; h=site_hit(ref,rows_,s)
        per=[h[phs==p].mean() for p in pdz.KEY_PHYLA if (phs==p).sum()>=30]; kc[vn]=(round(float(np.mean(per)),3),round(float(min(per)),3))
    ten={}
    for db,(ref,rr) in refs.items():
        h,valid=site_cov_ten(s,ref,rr); ten[db]=round(float(h[valid].mean()),3) if valid.sum() else np.nan
    tm=[round(cp.oligo_tm(p)[0],1) for p in s['prim']]
    rows.append(dict(site=g,oligos=len(s['prim']),pos=f"{int(P[P.name.isin([x['name'] for x in groups[g]])].start.min())}–{int(P[P.name.isin([x['name'] for x in groups[g]])].end.max())}",expansions=sum(len(pdz.expand(p)) for p in s['prim']),Tm=f"{min(tm)}–{max(tm)}",GG_mean=kc['GG'][0],GG_min=kc['GG'][1],SILVA_mean=kc['SILVA'][0],SILVA_min=kc['SILVA'][1],Ten_GG=ten['GG'],Ten_SILVA=ten['SILVA']))
C=pd.DataFrame(rows); print(C.to_string()); C.to_csv(SC+'snap_site_cov.csv',index=False)
pickle.dump(SITES,open(SC+'snap_sites.pkl','wb'))
