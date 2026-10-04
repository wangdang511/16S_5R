# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/extend_*
import pickle, collections, json, sys
import pandas as pd, numpy as np
sys.path.insert(0,'.')
import primer_design as pdz
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
res,allh=pickle.load(open(SC+'ext/offtarget_hits_延长候选.pkl','rb')); H=allh['延长候选']
cov=pd.read_csv(SC+'ext_cov.csv')
loc={mm:collections.defaultdict(set) for mm in (1,2,3)}
for n,pos,st,lab,mm,L in H:
    if '_' in n: continue
    b=lab.split('#')[0]
    for t in (1,2,3):
        if mm<=t: loc[t][b].add((n,pos//20,st))
mt=collections.defaultdict(list)
for n,pos,st,lab,mm,L in res[('延长候选','rCRS')]: mt[lab.split('#')[0]].append((pos,st,mm))
cov['nexp']=[len(pdz.expand(s)) for s in cov.seq]
for t in (1,2,3): cov[f'loci{t}']=[len(loc[t][n]) for n in cov.name]
cov['per_exp2']=(cov.loci2/cov.nexp).round(1); cov['rCRS_le4']=[sum(1 for h in mt[n] if h[2]<=4) for n in cov.name]
cov.to_csv(SC+'ext_full.csv',index=False)
pd.set_option('display.width',250)
print(cov[['name','nt','Tm','SILVA','GG','nexp','loci1','loci2','loci3','per_exp2','rCRS_le4']].to_string())
