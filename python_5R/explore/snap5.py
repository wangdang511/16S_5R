# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/snap_*
import pickle, json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore'); pd.set_option('display.width',250)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'ten1.py').read().split("rows=[]; masks={}")[0]
exec(src)
SNAP=pickle.load(open(SC+'snap_sites.pkl','rb')); REC=json.load(open(SC+'final_rec_sites.json')); orig=json.load(open(SC+'a5f_final5.json'))
cands={'SNAP V3_f':SNAP['V3_f'],'我们 A2-F（推荐，含补充）':REC['A2-F'],'我们 A2-F 原版 16 nt':dict(orig['A2-F'],prim=['YCTACGGGWGGCAGCA'],ext=[0],start=341,L=16)}
rows=[]
for nm,s in cands.items():
    for vn,(ref,rows_) in VAL.items():
        phs=ref.phylum[rows_]; h=site_hit(ref,rows_,s)
        for p in pdz.KEY_PHYLA:
            if (phs==p).sum()>=30: rows.append((nm,vn,p,int((phs==p).sum()),round(float(h[phs==p].mean()),3)))
R=pd.DataFrame(rows,columns=['site','db','phylum','n','cov']); P=R.pivot_table(index=['db','phylum','n'],columns='site',values='cov'); print(P.to_string())
