# 注：此脚本在本次会话中从 scratchpad 运行，依赖 S.pkl / big_labeled.pkl 等中间文件（未入库），路径需自行替换；结果见 docs/primer_design/v1v2_*.csv
import sys; sys.path.insert(0,'/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad')
from v2r_scan import *
def workF(args):
    st,L=args
    Md=[it.site_data(r,st,L) for r in (Gd,S12)]
    M=np.vstack([x[0] for x in Md]); phs=np.concatenate([x[1] for x in Md])
    if len(M)<500: return []
    cuts=np.cumsum([0]+[len(x[0]) for x in Md]); blocks=[slice(cuts[i],cuts[i+1]) for i in range(len(Md))]
    out=[]
    for fold,deg in FOLDS:
        tops=it.design_set(M,phs,'right',6,fold,deg)
        for k in range(1,len(tops)+1):
            worst,vals=cp.coverage_blocks(M,phs,tops[:k],'right',blocks)
            out.append(dict(start=st,L=L,fold=fold,deg=deg,k=k,cov_min=worst,cov_gd=vals[0],cov_s12=vals[1],total_fold=sum(len(pdz.expand(x)) for x in tops[:k]),tops=list(tops[:k])))
    return out
if __name__=='__main__':
    tasks=[(st,L) for st in range(4,15) for L in (17,19,21,23)]
    res=[]
    with Pool(4) as p:
        for r in p.imap_unordered(workF,tasks): res.extend(r)
    df=pd.DataFrame(res); df.to_pickle(SC+'v2f_scan.pkl')
    for k in (1,2,3,4):
        for fold in (4,8,16):
            d=df[(df.k==k)&(df.fold==fold)].sort_values('cov_min',ascending=False)
            if len(d): b=d.iloc[0]; print(k,fold,round(b.cov_min,3),b.start,b.L,b.tops,flush=True)
