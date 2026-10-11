# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/order_list_full.csv 与 final_*
import pickle, json, numpy as np, pandas as pd, warnings, os
warnings.filterwarnings('ignore'); pd.set_option('display.width',250); pd.set_option('display.max_colwidth',60)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'ten1.py').read().split("rows=[]; masks={}")[0]
exec(src)
sites_short,sites_long,sup=pickle.load(open(SC+'final1.pkl','rb'))
SUP={'A1-R':sup['short']['A1-R']['addon'],'A2-F':sup['short']['A2-F']['addon'],'A4-R':sup['short']['A4-R']['addon'],'A5-F':sup['short']['A5-F']['addon']}
def ext_sup(site,p,e):
    s=sites_long[site]; o=s['orient']
    # 补充寡核苷酸占据的窗口
    L=len(p); start=s['start'] if o=='R' else s['start']+s['L']-L   # F：3′ 端对齐；R：从位点起点开始
    if o=='F': start=s['start']+s['L']-L
    if e==0: return p
    cols=range(start-e,start) if o=='F' else range(start+L,start+L+e)
    # 该补充寡核苷酸命中的支原体序列
    cnt=np.zeros((len(cols),4)); 
    s2=dict(s,start=start,L=L,prim=[p],ext=[0])
    for db,(ref,rr) in refs.items():
        h,valid=site_cov_ten(s2,ref,rr)
        for j,c in enumerate(cols):
            col=ref.aln[rr[h],c-1]
            for k,b in enumerate(b'ACGT'): cnt[j,k]+=(col==b).sum()
    add=''.join('ACGT'[int(r.argmax())] for r in cnt)
    if o=='F': return add+p
    return pdz.revcomp(pdz.revcomp(p)+add)
rows=[]; cand={}
for site,p in SUP.items():
    for e in (0,2,3,4):
        q=ext_sup(site,p,e); s=sites_long[site]; o=s['orient']
        L=len(q); sdef=dict(s,start=(s['start']+s['L']-len(p)-e) if o=='F' else s['start'],L=len(q),prim=[q],ext=[0])
        # 位点覆盖（支原体，q 单独）
        r=dict(site=site,e=e,seq=q,nt=len(q),Tm=round(cp.oligo_tm(q)[0],1),nexp=len(pdz.expand(q)))
        for db,(ref,rr) in refs.items():
            comb=dict(s,prim=s['prim']+[q],ext=s['ext']+[0]) if e==0 else None
            # 组合覆盖：核心位点 + 该补充（延长后用其实际窗口）
            h0,v0=site_cov_ten(s,ref,rr)
            wst,wL=(s['start']+s['L']-len(p)-e,len(q)) if o=='F' else (s['start'],len(q))
            M=pdz.site_matrix(ref,wst,wL)[rr]; ok=(M!=0).all(1)&(M!=45).all(1); h1=np.zeros(len(rr),bool); h1[ok]=pdz._match(M[ok],top_of(o,q),'right' if o=='F' else 'left')
            r[db+'_core']=round(float(h0[v0].mean()),3); r[db+'_with']=round(float((h0|h1)[v0].mean()),3)
        rows.append(r); cand[f'{site}.sup.e{e}']=q
R=pd.DataFrame(rows); print(R.to_string()); R.to_csv(SC+'final_sup_ext.csv',index=False)
os.makedirs(SC+'fin2',exist_ok=True)
json.dump({'SUPEXT':[dict(name=k,seq=v) for k,v in cand.items()]},open(SC+'fin2/offtarget_primers.json','w'),ensure_ascii=False)
