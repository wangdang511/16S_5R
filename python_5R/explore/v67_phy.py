# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/v67_* 与 order_list_5amp_v67short.csv
import pickle, json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore'); pd.set_option('display.width',250)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'ext2.py').read().split("rows=[]")[0]
exec(src)
S=json.load(open(SC+'short_sites.json')); orig=json.load(open(SC+'a5f_final5.json')); D1=S['D1']
rows=[]
for nm,s in (('A4-F 原',orig['A4-F']),('A4-F 新',D1['A4-F']),('A4-R 原',orig['A4-R']),('A4-R 新',D1['A4-R'])):
    for vn,(ref,rows_) in VAL.items():
        phs=ref.phylum[rows_]; h=site_hit(ref,rows_,s)
        for g in pdz.KEY_PHYLA:
            if (phs==g).sum()>=30: rows.append((nm,vn,g,int((phs==g).sum()),round(float(h[phs==g].mean()),3)))
R=pd.DataFrame(rows,columns=['site','db','phylum','n','cov']); P=R.pivot_table(index=['db','phylum','n'],columns='site',values='cov'); print(P.to_string()); P.to_csv(SC+'v67_phylum.csv')
