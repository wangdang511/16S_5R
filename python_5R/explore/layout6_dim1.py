# 注：此脚本在本次会话中从 scratchpad 运行（借用 wangdang511/Olivar_primer 的 SADDLE 打分 PrimerSetBadnessFast），依赖未入库的中间文件；结果见 docs/primer_design/layout6_dimer_*.csv
import sys, types, importlib.util, numpy as np, pandas as pd, itertools, collections, warnings
warnings.filterwarnings('ignore')
O='/home/user/wangdang511/olivar_primer/src/olivar'
sys.path.insert(0,O)
try:
    import design as od
    print('olivar.design 导入成功')
except Exception as e:
    print('导入失败',repr(e)); raise
sys.path.insert(0,'/home/user/16S_5R/python_5R'); sys.path.insert(0,'/home/user/16S_5R/python_5R/explore')
import primer_design as pdz
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
# 1) Olivar 示例池基准
base={}
for k in (1,2):
    d=pd.read_csv(f'/home/user/wangdang511/olivar_primer/example_output/olivar-val_pool-{k}.csv')
    fP=[s for n,s in zip(d.name,d.seq) if n.endswith('_fP')]; rP=[s for n,s in zip(d.name,d.seq) if n.endswith('_rP')]
    tot,comp=od.PrimerSetBadnessFast(fP,rP); allc=np.array(comp[0]+comp[1])
    base[k]=(tot,len(fP)+len(rP),allc); print(f'Olivar 示例 pool-{k}: 引物 {len(fP)+len(rP)} 总坏度 {tot:.0f} 每条均值 {allc.mean():.1f} 中位 {np.median(allc):.1f} 最大 {allc.max():.0f}')
# 2) 我们的 32 条池：所有展开序列，按混合浓度加权
T=pd.read_pickle(SC+'pool32_table.pkl')
def vols(T,total=2.5):
    v={}
    for s,g in T.groupby('位点'):
        sup=g[g['名称'].str.endswith('s')]; main=g[~g['名称'].str.endswith('s')]
        for _,r in main.iterrows(): v[r['名称']]=total*(0.9 if len(sup) else 1)*r['展开数']/main['展开数'].sum()
        for _,r in sup.iterrows(): v[r['名称']]=total*0.1*r['展开数']/sup['展开数'].sum()
    return v
V=vols(T)
rows=[]
for _,r in T.iterrows():
    ex=pdz.expand(r['5′→3′序列']); 
    for k,e in enumerate(ex): rows.append(dict(oligo=r['名称'],site=r['位点'],ori=r['位点'][-1],seq=e.lower(),conc=V[r['名称']]/len(ex)))
E=pd.DataFrame(rows); E['conc_rel']=E.conc/E.conc.mean()
print('展开序列',len(E))
F=E[E.ori=='F']; R=E[E.ori=='R']
for lab,conc in (('等浓度（每条展开序列=1）',False),('按混合体积折算的相对浓度',True)):
    tot,comp=od.PrimerSetBadnessFast(list(F.seq),list(R.seq),fP_conc=list(F.conc_rel) if conc else None,rP_conc=list(R.conc_rel) if conc else None)
    E.loc[F.index,'bad']=comp[0]; E.loc[R.index,'bad']=comp[1]
    print(f'我们的池 {lab}: 总坏度 {tot:.0f}，每条展开序列均值 {E.bad.mean():.1f} 中位 {E.bad.median():.1f} 最大 {E.bad.max():.0f}')
    E['bad_'+('w' if conc else 'e')]=E.bad
E.to_pickle(SC+'dim_E.pkl')
pickle=__import__('pickle'); pickle.dump(base,open(SC+'dim_base.pkl','wb'))
# 每个寡核苷酸汇总（按浓度加权）
g=E.groupby('oligo').apply(lambda x: pd.Series(dict(site=x.site.iloc[0],n=len(x),bad_mean=x.bad_e.mean(),bad_max=x.bad_e.max())))
print(g.sort_values('bad_mean',ascending=False).head(12).round(1).to_string())
