# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/v5_* 与 order_list_full_v5.csv
import pandas as pd, numpy as np
import primer_design as pdz, compact_primers as cp
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
D=pd.read_csv('/home/user/16S_5R/docs/primer_design/order_list_full.csv')
D=D[D.order!='SMURF5-A4-F'].copy()          # 去掉 906 正向
def stats(p):
    ts=[cp.oligo_tm(x)[0] for x in pdz.expand(p)]; gc=sum(c in 'GCS' for c in p)/len(p); return round(cp.oligo_tm(p)[0],1),round(min(ts),1),round(max(ts),1),round(gc,2),len(pdz.expand(p))
new=[('SMURF5-A4-F.1','A4-F','V6·V7','F','核心','必订（V5 版）','MWACGCGARGAACCTTACC','967–985',1.0,0.54),
     ('SMURF5-A4-F.2','A4-F','V6·V7','F','核心','必订（V5 版）','MWACMCGAAGAACCTTACC','967–985',2.8,0.36),
     ('SMURF5-A4-F.3（补充）','A4-F','V6·V7','F','支原体补充','建议','RTACMCGAARAACCTTACC','967–985',1.8,0.10),
     ('SMURF5-V5-R.1','V5-R','V5','R','核心','必订（V5 版）','CCCGTCAATTYCTTTRAGTTT','907–927',2.8,0.5),
     ('SMURF5-V5-R.2','V5-R','V5','R','核心','必订（V5 版）','CCCGTCAATTCCTTTGAGYYY','907–927',0.4,0.5)]
rows=[]
for o,site,reg,ori,grp,rec,seq,pos,pe,sh in new:
    t,lo,hi,gc,ne=stats(seq); rows.append(dict(order=o,site=site,region=reg,orient=ori,group=grp,recommend=rec,seq=seq,nt=len(seq),pos=pos,Tm=t,Tm_lo=lo,Tm_hi=hi,GC=gc,expansions=ne,human_loci_le2=int(round(pe*ne)),human_per_exp=pe,is_sup=grp!='核心',share_with_sup=sh,share_no_sup=(sh if grp=='核心' else 0),uL_with_sup=round(sh*2.5,2),uL_no_sup=round((sh if grp=='核心' else 0)*2.5,2)))
D=pd.concat([D,pd.DataFrame(rows)],ignore_index=True)
# 不含补充时 A4-F 核心 .6/.4
D.loc[D.order=='SMURF5-A4-F.1','share_no_sup']=0.6; D.loc[D.order=='SMURF5-A4-F.1','uL_no_sup']=1.5
D.loc[D.order=='SMURF5-A4-F.2','share_no_sup']=0.4; D.loc[D.order=='SMURF5-A4-F.2','uL_no_sup']=1.0
# 备选里的 967 缩短版不再是备选，改成备选 906
D=D[~D.order.str.contains('alt967')]
D=pd.concat([D,pd.DataFrame([dict(order='SMURF5-A4-F.alt906',site='A4-F',region='V6·V7',orient='F',group='备选',recommend='不加 V5 时用：V6·V7 的 906 正向引物（扩增子 906–1195 = 290 bp）；与 V5-R 的 907 位点互补，二者不能同时订',seq='RAAACTCAAAGGAATTGACGGR',nt=22,pos='906–927',Tm=60.5,Tm_lo=59.9,Tm_hi=61.2,GC=0.36,expansions=4,human_loci_le2=2,human_per_exp=0.5,is_sup=False)])],ignore_index=True)
D.to_csv(SC+'order_list_full_v5.csv',index=False)
print(D[['order','recommend','seq','pos','expansions','uL_with_sup','uL_no_sup']].to_string()); print('总 µL',D.uL_with_sup.sum(),D.uL_no_sup.sum())
sev,hp=cp.pool_issues({r.order:[r.seq] for r in D[D.group!='备选'].itertuples()}); print('推荐池严重二聚体',sev,hp)
