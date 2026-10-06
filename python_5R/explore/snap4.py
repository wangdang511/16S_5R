# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/snap_*
import pickle, json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore'); pd.set_option('display.width',250); pd.set_option('display.max_colwidth',40)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'ten1.py').read().split("rows=[]; masks={}")[0]
exec(src)
SNAP=pickle.load(open(SC+'snap_sites.pkl','rb')); REC=json.load(open(SC+'final_rec_sites.json'))
five=[("R1-F","TGGCGAACGGGTGAGTAA","F",103),("R1-R","CCGTGTCTCAGTCCCARTG","R",314),("R2-F","ACTCCTACGGGAGGCAGC","F",338),("R2-R","GTATTACCGCGGCTGCTG","R",519),("R3-F","GTGTAGCGGTGRAATGCG","F",685),("R3-R","CCCGTCAATTCMTTTGAGTT","R",908),("R4-F","GGAGCATGTGGWTTAATTCGA","F",944),("R4-R","CGTTGCGGGACTTAACCC","R",1087),("R5-F","GGAGGAAGGTGGGGATGAC","F",1175),("R5-R","AAGGCCCGGGAACGTATT","R",1374)]
S5={n:dict(orient=o,start=p,L=len(s),prim=[s],ext=[0]) for n,s,o,p in five}
def cov(s):
    kc={}
    for vn,(ref,rows_) in VAL.items():
        phs=ref.phylum[rows_]; h=site_hit(ref,rows_,s); per=[h[phs==p].mean() for p in pdz.KEY_PHYLA if (phs==p).sum()>=30]; kc[vn]=(float(np.mean(per)),float(min(per)))
    ten={}
    for db,(ref,rr) in refs.items():
        h,valid=site_cov_ten(s,ref,rr); ten[db]=float(h[valid].mean()) if valid.sum() else np.nan
    return kc,ten
rows=[]
for scheme,sites in (('Swift SNAP',SNAP),('5R 原方案',S5),('我们的推荐设计',REC)):
    for n,s in sites.items():
        kc,ten=cov(s); rows.append(dict(scheme=scheme,site=n,oligos=len(s['prim']),expansions=sum(len(pdz.expand(p)) for p in s['prim']),GG_mean=kc['GG'][0],GG_min=kc['GG'][1],SILVA_mean=kc['SILVA'][0],SILVA_min=kc['SILVA'][1],Ten_GG=ten['GG'],Ten_SILVA=ten['SILVA']))
R=pd.DataFrame(rows); R.to_csv(SC+'snap_all_sites.csv',index=False)
S=R.groupby('scheme').agg(sites=('site','size'),oligos=('oligos','sum'),exp=('expansions','sum'),GGm=('GG_mean','mean'),GGmin=('GG_min','min'),SILm=('SILVA_mean','mean'),SILmin=('SILVA_min','min'),Ten_GG=('Ten_GG','mean'),Ten_SIL=('Ten_SILVA','mean'),n_weak=('GG_mean',lambda x:(x<0.8).sum())).round(3)
print(S.to_string()); S.to_csv(SC+'snap_scheme_summary.csv')
print(R.round(3).to_string())
