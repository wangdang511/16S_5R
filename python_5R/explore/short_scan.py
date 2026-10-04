# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/short_*
import sys; sys.path.insert(0,'/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad')
from v2r_scan import *
def work(args):
    kind,pos,L=args           # F: pos=3′ 端位置；R: pos=位点起点（顶链左端）
    st=pos-L+1 if kind[-1]=='F' else pos
    Md=[it.site_data(r,st,L) for r in (Gd,S12)]
    M=np.vstack([x[0] for x in Md]); phs=np.concatenate([x[1] for x in Md])
    if len(M)<500: return []
    cuts=np.cumsum([0]+[len(x[0]) for x in Md]); blocks=[slice(cuts[i],cuts[i+1]) for i in range(len(Md))]
    side='right' if kind[-1]=='F' else 'left'; out=[]
    for fold,deg in ((4,2),(8,3)):
        tops=it.design_set(M,phs,side,2,fold,deg)
        for k in range(1,len(tops)+1):
            worst,vals=cp.coverage_blocks(M,phs,tops[:k],side,blocks)
            prim=list(tops[:k]) if side=='right' else [pdz.revcomp(x) for x in tops[:k]]; sev,hp=cp.pool_issues({'x':prim})
            out.append(dict(kind=kind,pos=pos,L=L,start=st,fold=fold,k=k,cov_min=worst,cov_gd=vals[0],cov_s12=vals[1],prim=prim,tm=[round(cp.oligo_tm(p)[0],1) for p in prim],nexp=sum(len(pdz.expand(p)) for p in prim),severe=len(sev)))
    return out
if __name__=='__main__':
    tasks=[]
    for L in (19,20):
        tasks+= [('A4F',p,L) for p in range(927,1001)]+[('A4R',p,L) for p in range(1100,1178)]
        tasks+= [('A5F',p,L) for p in range(1241,1291)]+[('A5R',p,L) for p in range(1380,1492)]
    res=[]
    with Pool(4) as p:
        for r in p.imap_unordered(work,tasks): res.extend(r)
    df=pd.DataFrame(res); df['tmmax']=df.tm.apply(max); df['tmmin']=df.tm.apply(min); df.to_pickle(SC+'short_scan.pkl'); print(len(df))
