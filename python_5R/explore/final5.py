# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/order_list_full.csv 与 final_*
import pickle, json, numpy as np, pandas as pd, warnings, os
warnings.filterwarnings('ignore'); pd.set_option('display.width',250); pd.set_option('display.max_colwidth',60)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'ten1.py').read().split("rows=[]; masks={}")[0]
exec(src)
src2=open(SC+'ext2.py').read().split("rows=[]")[0].split("src=open")[0]
exec(open(SC+'ext2.py').read().split("rows=[]")[0].replace("exec(src)","exec(src)",1)) if False else None
sites_short,sites_long,sup=pickle.load(open(SC+'final1.pkl','rb'))
E=pd.read_csv(SC+'final_sup_ext_full.csv').set_index(['site','e'])
CH={'A2-F':3,'A5-F':3,'A4-R':2,'A1-R':3}
SUPSEQ={k:E.loc[(k,e)].seq for k,e in CH.items()}
rec=json.loads(json.dumps(sites_long))
# 位点内第 i 条的 5′/3′ 坐标
def coords(s,p,e,sup_len=None):
    o=s['orient']
    if sup_len is None:
        return (s['start']-e,s['start']+s['L']-1) if o=='F' else (s['start'],s['start']+s['L']-1+e)
    return (s['start']+s['L']-sup_len,s['start']+s['L']-1) if o=='F' else (s['start'],s['start']+sup_len-1)
rows=[]
REGION={'A1':'V1·V2','A2':'V3','A3':'V4','A4':'V6·V7','A5':'V8·V9'}
pickle.dump((rec,SUPSEQ,CH),open(SC+'final5.pkl','wb'))
print(SUPSEQ)
