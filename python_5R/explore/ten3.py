# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/ten_*
import pickle, json, numpy as np, pandas as pd, warnings, collections
warnings.filterwarnings('ignore'); pd.set_option('display.width',250); pd.set_option('display.max_colwidth',60)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'ten1.py').read().split("rows=[]; masks={}")[0]
exec(src)
D1j=D1
def addon(site,fold=8,deg=3):
    s=D1j[site]; res={}
    Ms=[];keep=[]
    for db,(ref,rr) in refs.items():
        h,valid=site_cov_ten(s,ref,rr); st,Ln=window(s['orient'],s['start'],s['L'],0)
        M=pdz.site_matrix(ref,st,Ln)[rr]; sel=valid&~h; Ms.append(M[sel])
    M=np.vstack(Ms); 
    side='right' if s['orient']=='F' else 'left'
    top,_=pdz.best_degenerate(M,side,deg,fold)
    prim=top if s['orient']=='F' else pdz.revcomp(top)
    s2=dict(s,prim=s['prim']+[prim],ext=s['ext']+[0]); out=dict(site=site,n_uncovered=len(M),addon=prim,nexp=len(pdz.expand(prim)),Tm=round(cp.oligo_tm(prim)[0],1))
    for db,(ref,rr) in refs.items():
        h,valid=site_cov_ten(s,ref,rr); h2,_=site_cov_ten(s2,ref,rr); out[db+'_before']=round(float(h[valid].mean()),3); out[db+'_after']=round(float(h2[valid].mean()),3)
    return out
rows=[addon(x) for x in ('A2-R','A4-F','A5-F','A1-R','A4-R','A2-F')]
R=pd.DataFrame(rows); print(R.to_string()); R.to_csv(SC+'ten_addon.csv',index=False)
# 加了补充寡核苷酸后的扩增子覆盖（支原体）
add={r['site']:r['addon'] for r in rows}
def mk(site,use_add):
    s=D1j[site]
    return dict(s,prim=s['prim']+([add[site]] if use_add and site in add else []),ext=s['ext']+([0] if use_add and site in add else []))
amps={'V1·V2':('A1-F','A1-R'),'V3':('A2-F','A2-R'),'V4':('A3-F','A3-R'),'V6·V7（缩短）':('A4-F','A4-R'),'V8·V9':('A5-F','A5-R')}
out=[]
for nm,(f,r) in amps.items():
    row=dict(amplicon=nm)
    for db,(ref,rr) in refs.items():
        for tag,ua in (('现有',False),('加补充',True)):
            hf,vf=site_cov_ten(mk(f,ua),ref,rr); hr,vr=site_cov_ten(mk(r,ua),ref,rr); v=vf&vr; row[f'{db}_{tag}']=round(float((hf&hr)[v].mean()),3)
    out.append(row)
O=pd.DataFrame(out); print(O.to_string()); O.to_csv(SC+'ten_addon_amp.csv',index=False)
