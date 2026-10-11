# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/layout6_*
import sys; sys.path.insert(0,'/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad')
from v2r_scan import *
def workF(args):
    st,L=args
    Md=[it.site_data(r,st,L) for r in (Gd,S12)]
    M=np.vstack([x[0] for x in Md]); phs=np.concatenate([x[1] for x in Md])
    if len(M)<500: return []
    cuts=np.cumsum([0]+[len(x[0]) for x in Md]); blocks=[slice(cuts[i],cuts[i+1]) for i in range(len(Md))]
    out=[]
    for fold,deg in ((4,2),(8,3)):
        tops=it.design_set(M,phs,'right',3,fold,deg)
        for k in range(1,len(tops)+1):
            worst,vals=cp.coverage_blocks(M,phs,tops[:k],'right',blocks); prim=list(tops[:k]); sev,hp=cp.pool_issues({'x':prim})
            out.append(dict(o='F',start=st,L=L,fold=fold,k=k,cov_min=worst,cov_gd=vals[0],cov_s12=vals[1],prim=prim,tm=[round(cp.oligo_tm(p)[0],1) for p in prim],nexp=sum(len(pdz.expand(p)) for p in prim),severe=len(sev)))
    return out
def workR(args):
    st,L=args
    Md=[it.site_data(r,st,L) for r in (Gd,S12)]
    M=np.vstack([x[0] for x in Md]); phs=np.concatenate([x[1] for x in Md])
    if len(M)<500: return []
    cuts=np.cumsum([0]+[len(x[0]) for x in Md]); blocks=[slice(cuts[i],cuts[i+1]) for i in range(len(Md))]
    out=[]
    for fold,deg in ((4,2),(8,3)):
        tops=it.design_set(M,phs,'left',3,fold,deg)
        for k in range(1,len(tops)+1):
            worst,vals=cp.coverage_blocks(M,phs,tops[:k],'left',blocks); prim=[pdz.revcomp(x) for x in tops[:k]]; sev,hp=cp.pool_issues({'x':prim})
            out.append(dict(o='R',start=st,L=L,fold=fold,k=k,cov_min=worst,cov_gd=vals[0],cov_s12=vals[1],prim=prim,tm=[round(cp.oligo_tm(p)[0],1) for p in prim],nexp=sum(len(pdz.expand(p)) for p in prim),severe=len(sev)))
    return out
def work2(t): return workF(t[1:]) if t[0]=='F' else workR(t[1:])
if __name__=='__main__':
    # F：5' 起点 806–1030（每 1 nt，长度 18/20）; R：5' 端（在 R 链）起点 = 区间起，扫 960–1060
    tasks=[('F',st,L) for st in range(806,1031) for L in (18,20)]+[('R',st,L) for st in range(940,1062) for L in (18,20)]
    res=[]
    with Pool(4) as p:
        for r in p.imap_unordered(work2,tasks): res.extend(r)
    df=pd.DataFrame(res); df.to_pickle(SC+'v5d_scan.pkl'); print(len(df))
