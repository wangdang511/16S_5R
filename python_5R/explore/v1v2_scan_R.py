# 注：此脚本在本次会话中从 scratchpad 运行，依赖 S.pkl / big_labeled.pkl 等中间文件（未入库），路径需自行替换；结果见 docs/primer_design/v1v2_*.csv
import pickle, json, time, numpy as np, pandas as pd, warnings
from multiprocessing import Pool
import primer_design as pdz, iterate_5R as it, compact_primers as cp
warnings.filterwarnings('ignore')
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
S=pickle.load(open(SC+'S.pkl','rb')); B=pickle.load(open(SC+'big_labeled.pkl','rb'))
tax=B['tax']; ids=B['ids']; ph=np.array([tax[i][1] if tax[i][1] else 'unclassified' for i in ids])
G=pdz.Reference(B['aln'],ids,np.array(['Bacteria']*len(ids)),ph,S.col_of_pos,S.ref_seq,{})
rng=np.random.default_rng(3); idx12=np.sort(rng.choice(len(S.ids),12000,replace=False))
sub=lambda ref,ix:pdz.Reference(ref.aln[ix],[ref.ids[i] for i in ix],ref.domain[ix],ref.phylum[ix],ref.col_of_pos,ref.ref_seq,{})
S12=sub(S,idx12); gi=np.random.default_rng(1).permutation(len(ids)); Gd=sub(G,gi[:len(ids)//2])
FOLDS=((4,2),(8,3),(16,4))
def work(args):
    st,L=args
    Md=[it.site_data(r,st,L) for r in (Gd,S12)]
    M=np.vstack([x[0] for x in Md]); phs=np.concatenate([x[1] for x in Md])
    if len(M)<500: return []
    cuts=np.cumsum([0]+[len(x[0]) for x in Md]); blocks=[slice(cuts[i],cuts[i+1]) for i in range(len(Md))]
    out=[]
    for fold,deg in FOLDS:
        tops=it.design_set(M,phs,'left',8,fold,deg)
        for k in range(1,len(tops)+1):
            worst,vals=cp.coverage_blocks(M,phs,tops[:k],'left',blocks)
            out.append(dict(start=st,L=L,fold=fold,deg=deg,k=k,cov_min=worst,cov_gd=vals[0],cov_s12=vals[1],total_fold=sum(len(pdz.expand(x)) for x in tops[:k]),
                            tops=[pdz.revcomp(x) for x in tops[:k]]))
    return out
if __name__=='__main__':
    tasks=[(st,L) for st in range(243,281) for L in (17,19,21,23)]
    res=[]; t=time.time()
    with Pool(4) as p:
        for r in p.imap_unordered(work,tasks): res.extend(r)
    df=pd.DataFrame(res); df.to_pickle(SC+'v2r_scan.pkl')
    print('configs',len(df),round(time.time()-t))
    best=df.sort_values('cov_min',ascending=False)
    for k in (1,2,3,4,6,8):
        b=df[df.k==k].sort_values('cov_min',ascending=False).iloc[0]
        print(f'k={k}: best cov_min {b.cov_min:.3f} (GG-design {b.cov_gd:.3f}, SILVA {b.cov_s12:.3f}) at {b.start} L{b.L} fold<= {b.fold} total_fold {b.total_fold}')
