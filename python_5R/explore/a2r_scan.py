# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/a2r_*
import sys; sys.path.insert(0,'/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad')
from v2r_scan import *
def workA(args):
    st,L=args
    Md=[it.site_data(r,st,L) for r in (Gd,S12)]
    M=np.vstack([x[0] for x in Md]); phs=np.concatenate([x[1] for x in Md])
    if len(M)<500: return []
    cuts=np.cumsum([0]+[len(x[0]) for x in Md]); blocks=[slice(cuts[i],cuts[i+1]) for i in range(len(Md))]
    out=[]
    for fold,deg in ((4,2),(8,3)):
        tops=it.design_set(M,phs,'left',2,fold,deg)
        for k in range(1,len(tops)+1):
            worst,vals=cp.coverage_blocks(M,phs,tops[:k],'left',blocks)
            prim=[pdz.revcomp(x) for x in tops[:k]]
            out.append(dict(start=st,L=L,fold=fold,k=k,cov_min=worst,cov_gd=vals[0],cov_s12=vals[1],prim=prim,
                            tm=[round(cp.oligo_tm(p)[0],1) for p in prim],nexp=sum(len(pdz.expand(p)) for p in prim)))
    return out
if __name__=='__main__':
    tasks=[(st,L) for st in range(488,513) for L in range(17,24)]
    res=[]
    with Pool(4) as p:
        for r in p.imap_unordered(workA,tasks): res.extend(r)
    df=pd.DataFrame(res); df['tmmax']=df.tm.apply(max); df['tmmin']=df.tm.apply(min); df.to_pickle(SC+'a2r_scan.pkl')
    g=df[(df.cov_min>=0.93)&(df.tmmax<=65)&(df.tmmin>=57)].sort_values(['k','nexp','cov_min'],ascending=[True,True,False])
    print(len(df),len(g)); print(g.head(25)[['start','L','fold','k','cov_min','cov_gd','cov_s12','tm','nexp','prim']].to_string())
