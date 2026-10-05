# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/ten_*
import pickle, numpy as np, pandas as pd, collections, warnings
warnings.filterwarnings('ignore'); pd.set_option('display.width',250)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
masks,gen0,fam0,ggrows=pickle.load(open(SC+'ten_masks.pkl','rb'))
def pair(a,b,db):
    ha,va=masks[(a,db)]; hb,vb=masks[(b,db)]; v=va&vb; return (ha&hb)[v].mean() if v.sum() else np.nan, int(v.sum())
amps={'V1·V2':('A1-F','A1-R'),'V3':('A2-F','A2-R'),'V4':('A3-F','A3-R'),'V6·V7（缩短，967）':('A4-F','A4-R'),'V6·V7（原，906）':('A4-F(原 906–927)','A4-R'),'V8·V9':('A5-F','A5-R')}
rows=[]
for nm,(f,r) in amps.items():
    row=dict(amplicon=nm)
    for db in ('GG','SILVA'): c,n=pair(f,r,db); row[db]=round(c,3); row[db+'_n']=n
    rows.append(row)
print(pd.DataFrame(rows).to_string())
# 5R 各扩增子
for i in range(1,6):
    row={}
    for db in ('GG','SILVA'): c,n=pair(f'5R R{i}-F',f'5R R{i}-R',db); row[db]=round(c,3)
    print('5R R%d'%i,row)
# 全部可用扩增子（各自都能扩出的比例，GG，用各位点 valid 的交集）
for db in ('GG','SILVA'):
    sets=[('V1·V2',('A1-F','A1-R')),('V3',('A2-F','A2-R')),('V4',('A3-F','A3-R')),('V6V7',('A4-F','A4-R')),('V8V9',('A5-F','A5-R'))]
    n=len(masks[('A2-F',db)][0]); H=np.zeros((n,len(sets)),bool)
    for j,(nm,(f,r)) in enumerate(sets): H[:,j]=masks[(f,db)][0]&masks[(r,db)][0]
    print(db,'扩出的扩增子数分布（0..5）:',np.bincount(H.sum(1),minlength=6).tolist(),'至少 3 个：%.2f 至少 1 个：%.2f'%((H.sum(1)>=3).mean(),(H.sum(1)>=1).mean()))
# GG 属/科分解：每个位点未被覆盖的 Tenericutes 属
db='GG'
for site in ('A2-R','A4-F','A5-F','A1-R','A4-R','A2-F'):
    h,v=masks[(site,db)]; miss=~h&v; c=collections.Counter(gen0[ggrows][miss]); tot=collections.Counter(gen0[ggrows][v])
    print(site,'未覆盖',int(miss.sum()),'/',int(v.sum()),'前几位属:',[(g if g else '(无属名)',m,tot[g]) for g,m in c.most_common(6)])
print('GG Tenericutes 属组成:',collections.Counter(gen0[ggrows]).most_common(8))
rows=[]
for nm,(f,r) in amps.items():
    row=dict(amplicon=nm)
    for db in ('GG','SILVA'): c,n=pair(f,r,db); row[db]=round(c,3); row[db+'_n']=n
    rows.append(row)
for i in range(1,6):
    row=dict(amplicon=f'5R R{i}')
    for db in ('GG','SILVA'): c,n=pair(f'5R R{i}-F',f'5R R{i}-R',db); row[db]=round(c,3); row[db+'_n']=n
    rows.append(row)
pd.DataFrame(rows).to_csv(SC+'ten_amp_cov.csv',index=False)
gen=[]
for site in ('A2-R','A4-F','A5-F','A1-R','A4-R','A2-F','A4-F(原 906–927)'):
    h,v=masks[(site,'GG')]; miss=~h&v; c=collections.Counter(gen0[ggrows][miss]); tot=collections.Counter(gen0[ggrows][v])
    for g in ('Mycoplasma','Candidatus Phytoplasma','Asteroleplasma','Acholeplasma','Ureaplasma'):
        gen.append(dict(site=site,genus=g,n=tot[g],covered=tot[g]-c[g]))
pd.DataFrame(gen).to_csv(SC+'ten_genus_cov.csv',index=False)
