# 注：此脚本在本次会话中从 scratchpad 运行（借用 wangdang511/Olivar_primer 的 SADDLE 打分 PrimerSetBadnessFast），依赖未入库的中间文件；结果见 docs/primer_design/layout6_dimer_*.csv
import sys, numpy as np, pandas as pd, collections, primer3
sys.path.insert(0,'/home/user/16S_5R/python_5R')
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
E=pd.read_pickle(SC+'dim_E.pkl').reset_index(drop=True)
kw=dict(mv_conc=50,dv_conc=2,dntp_conc=0.2,dna_conc=250,temp_c=60)
seqs=[s.upper() for s in E.seq]; oli=list(E.oligo)
n=len(seqs); hd=np.zeros((n,n)); es=np.zeros((n,n))
for i in range(n):
    for j in range(i,n):
        hd[i,j]=hd[j,i]=primer3.bindings.calc_heterodimer(seqs[i],seqs[j],**kw).dg/1000.0
        # 3′端稳定性：i 的 3′ 端与 j 配对（不对称，取两个方向最差）
        a=primer3.bindings.calc_end_stability(seqs[i],seqs[j],**kw).dg/1000.0; b=primer3.bindings.calc_end_stability(seqs[j],seqs[i],**kw).dg/1000.0
        es[i,j]=a; es[j,i]=b
np.save(SC+'dim_hd.npy',hd); np.save(SC+'dim_es.npy',es)
print('异二聚体 ΔG（kcal/mol，60 °C）：最差',hd.min().round(1),'; 3′端稳定性最差',es.min().round(1))
# 阈值统计
for th in (-6,-7,-8,-9,-10): print(f'ΔG≤{th}: 全局二聚体 {int((np.triu(hd)<=th).sum())} 对；3′端锚定 {int((es<=th).sum())} 对（有序）')
# 按寡核苷酸对汇总最差
rows=collections.defaultdict(lambda:[0,0])
for i in range(n):
    for j in range(i,n):
        k=tuple(sorted((oli[i],oli[j]))); r=rows[k]
        r[0]=min(r[0],hd[i,j]); r[1]=min(r[1],es[i,j],es[j,i])
df=pd.DataFrame([(a,b,v[0],v[1]) for (a,b),v in rows.items()],columns=['A','B','dG_dimer','dG_3end'])
print('\n全局二聚体最差（kcal/mol）:'); print(df.sort_values('dG_dimer').head(12).round(1).to_string(index=False))
print('\n3′端锚定（可被聚合酶延伸）最差:'); print(df.sort_values('dG_3end').head(12).round(1).to_string(index=False))
df.to_csv(SC+'dim_pairs.csv',index=False)
