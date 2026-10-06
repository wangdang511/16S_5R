# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/
import pickle, json, numpy as np, pandas as pd, warnings, os
warnings.filterwarnings('ignore'); pd.set_option('display.width',250); pd.set_option('display.max_colwidth',50)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'ten1.py').read().split("rows=[]; masks={}")[0]
exec(src)
REC=json.load(open(SC+'final_rec_sites.json'))
base=REC['A2-F']; sup=base['prim'][1]
cands={'原推荐 ACT+核心（19 nt）':'ACTYCTACGGGWGGCAGCA','ACW+核心（19 nt）':'ACWYCTACGGGWGGCAGCA','CW+核心（18 nt）':'CWYCTACGGGWGGCAGCA','W+核心（17 nt）':'WYCTACGGGWGGCAGCA','核心 16 nt':'YCTACGGGWGGCAGCA','AMW+核心（19 nt）':'AMWYCTACGGGWGGCAGCA'}
rows=[]
for nm,q in cands.items():
    e=len(q)-16; s=dict(base,prim=[q,sup],ext=[e,3])
    r=dict(name=nm,seq=q,nt=len(q),Tm=round(cp.oligo_tm(q)[0],1),nexp=len(pdz.expand(q)))
    for vn,(ref,rows_) in VAL.items():
        phs=ref.phylum[rows_]
        for tag,ss in (('仅主引物',dict(base,prim=[q],ext=[e])),('含补充',s)):
            h=site_hit(ref,rows_,ss); per=[h[phs==p].mean() for p in pdz.KEY_PHYLA if (phs==p).sum()>=30]
            r[f'{vn}_{tag}_mean']=round(float(np.mean(per)),3); r[f'{vn}_{tag}_min']=round(float(min(per)),3)
    rows.append(r)
R=pd.DataFrame(rows); print(R.to_string()); R.to_csv(SC+'a2f_fix2.csv',index=False)
