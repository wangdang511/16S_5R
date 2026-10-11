# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/a2r_*
import pickle, json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore'); pd.set_option('display.width',250)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'ext2.py').read().split("rows=[]")[0]
exec(src)
df=pd.read_pickle(SC+'a2r_scan.pkl')
pick=[(505,17,8,1),(505,18,8,1),(505,19,8,1),(508,18,8,1),(505,18,4,2),(508,17,4,2),(509,17,4,2),(509,18,4,2),(510,17,4,2)]
cands={}
for st,L,fold,k in pick:
    r=df[(df.start==st)&(df.L==L)&(df.fold==fold)&(df.k==k)].iloc[0]
    cands[f'A2R_{st}_L{L}_k{k}_f{fold}']=dict(orient='R',start=st,L=L,prim=list(r.prim),ext=[0]*k)
old=pool['A2-R']; cands['A2R_旧512_L16']=old
ext_old=dict(old,prim=[cand['A2-R|0|e2']['seq']],ext=[2]); cands['A2R_旧512_延长+2']=ext_old
rows=[]
for nm,s in cands.items():
    kc=keycov(s); tms=[round(cp.oligo_tm(p)[0],1) for p in s['prim']]
    sev,hp=cp.pool_issues({'A2-R':s['prim']})
    rows.append(dict(name=nm,start=s['start'],L=s['L'],prim=' | '.join(s['prim']),Tm=tms,nexp=sum(len(pdz.expand(p)) for p in s['prim']),SILVA=kc['SILVA'],GG=kc['GG'],held=round(float(site_hit(Gh,np.arange(len(held)),s).mean()),3),self_severe=len(sev),hp=len(hp)))
R=pd.DataFrame(rows); print(R.to_string()); R.to_csv(SC+'a2r_eval.csv',index=False)
pickle.dump(cands,open(SC+'a2r_cands.pkl','wb'))
json.dump({'A2R候选':[dict(name=f'{n}_{i}',seq=p) for n,s in cands.items() for i,p in enumerate(s['prim'])]},open(SC+'a2r_primers.json','w'),ensure_ascii=False)
