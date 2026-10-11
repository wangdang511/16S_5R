# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/order_list_full.csv 与 final_*
import pickle, json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore'); pd.set_option('display.width',250)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'ten1.py').read().split("rows=[]; masks={}")[0]
exec(src)
cfgs=pickle.load(open(SC+'final_cfgs.pkl','rb'))
rows=[]
for nm,sites in cfgs.items():
    row=dict(design=nm)
    for db,(ref,rr) in refs.items():
        H=[]; V=[]
        for k in range(1,6):
            hf,vf=site_cov_ten(sites[f'A{k}-F'],ref,rr); hr,vr=site_cov_ten(sites[f'A{k}-R'],ref,rr); H.append(hf&hr); V.append(vf&vr)
            row[f'{db}_A{k}']=round(float((hf&hr)[vf&vr].mean()),3)
        Hm=np.array(H).T; row[db+'_>=3扩增子']=round(float((Hm.sum(1)>=3).mean()),3); row[db+'_>=1']=round(float((Hm.sum(1)>=1).mean()),3)
    rows.append(row)
R=pd.DataFrame(rows); print(R.T.to_string()); R.to_csv(SC+'final_ten_eval.csv',index=False)
