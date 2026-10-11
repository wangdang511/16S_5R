# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/extend_*
import pickle, json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'v2_tm.py').read().split("print('pool Tm")[0]
exec(src)
sj=json.load(open(SC+'v2_tm_sites.json'))
pool={**c4, 'A1-F':sj['A1-F'],'A1-R':sj['A1-R_balanced']}
cand={}; rows=[]
for sn,s in pool.items():
    for i,(p,e0) in enumerate(zip(s['prim'],s['ext'])):
        if len(p)>17: continue
        for e in range(0,7):
            q=p if e==0 else extend(s['orient'],s['start'],s['L'],p,e,[Gd,S12])
            name=f'{sn}|{i}|e{e}'; cand[name]=dict(site=sn,i=i,e=e,seq=q,orient=s['orient'],start=s['start'],L=s['L'],tm=cp.oligo_tm(q)[0])
pickle.dump((cand,pool),open(SC+'ext_cand.pkl','wb'))
json.dump({'延长候选':[dict(name=k,seq=v['seq']) for k,v in cand.items()]},open(SC+'offtarget_primers_ext.json','w'),ensure_ascii=False)
print(len(cand)); 
for k,v in list(cand.items())[:14]: print(k,v['seq'],round(v['tm'],1))
