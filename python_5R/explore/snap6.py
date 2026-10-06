# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/snap_*
import pickle, json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore'); pd.set_option('display.width',250)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
exec(open(SC+'ext2.py').read().split("rows=[]")[0])
SNAP=pickle.load(open(SC+'snap_sites.pkl','rb')); REC=json.load(open(SC+'final_rec_sites.json'))
five=[("R1-F","TGGCGAACGGGTGAGTAA","F",103),("R1-R","CCGTGTCTCAGTCCCARTG","R",314),("R2-F","ACTCCTACGGGAGGCAGC","F",338),("R2-R","GTATTACCGCGGCTGCTG","R",519),("R3-F","GTGTAGCGGTGRAATGCG","F",685),("R3-R","CCCGTCAATTCMTTTGAGTT","R",908),("R4-F","GGAGCATGTGGWTTAATTCGA","F",944),("R4-R","CGTTGCGGGACTTAACCC","R",1087),("R5-F","GGAGGAAGGTGGGGATGAC","F",1175),("R5-R","AAGGCCCGGGAACGTATT","R",1374)]
S5={n:dict(orient=o,start=p,L=len(s),prim=[s],ext=[0]) for n,s,o,p in five}
schemes={'Swift SNAP（假设配对）':(SNAP,[('V1_f','V2_r'),('V3_f','V4_r'),('V4_f','V5_r'),('V6_f','V8_r'),("V6'_f",'V8_r'),('V7_f','V9_r')]),
 '5R 原方案':(S5,[(f'R{i}-F',f'R{i}-R') for i in range(1,6)]),
 '我们的推荐设计':(REC,[(f'A{i}-F',f'A{i}-R') for i in range(1,6)])}
rows=[]
for nm,(sites,pairs) in schemes.items():
    r=dict(scheme=nm,pairs=len(pairs))
    for vn,(ref,rows_) in VAL.items():
        phs=ref.phylum[rows_]; H=[]
        for f,rr in pairs: H.append(site_hit(ref,rows_,sites[f])&site_hit(ref,rows_,sites[rr]))
        H=np.array(H).T; per=[]
        for p in pdz.KEY_PHYLA:
            if (phs==p).sum()>=30: per.append(H[phs==p].mean(0))
        per=np.array(per)                  # phylum × pair
        r[vn+'_pair_mean']=round(float(per.mean()),3); r[vn+'_pair_worst']=round(float(per.min()),3)
        r[vn+'_ge_half']=round(float((H.sum(1)>=np.ceil(len(pairs)/2)).mean()),3); r[vn+'_all']=round(float(H.all(1).mean()),3)
    rows.append(r)
R=pd.DataFrame(rows); print(R.T.to_string()); R.to_csv(SC+'snap_pair_cov.csv',index=False)
