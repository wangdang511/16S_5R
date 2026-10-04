# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/short_*
import pickle, json, numpy as np, pandas as pd, warnings, itertools
warnings.filterwarnings('ignore'); pd.set_option('display.width',250)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'ext2.py').read().split("rows=[]")[0]
exec(src)
sites=json.load(open(SC+'a5f_final5.json'))
def interior(k):
    f,r=sites[f'A{k}-F'],sites[f'A{k}-R']
    return (f['start']+f['L'],r['start']-1)   # 与 evaluate 里一致：F 引物 3′ 端之后 到 R 引物起点之前
ints0=[interior(k) for k in range(1,6)]; print(ints0)
base=ts.accuracy_direct(xall,ints0,gall,vall); print('base ideal',round(base,4), len(gall))
res=[]
for a in range(0,3):
    pass
A4=ints0[3]; A5=ints0[4]
rows=[]
for s4 in range(A4[0],A4[0]+90,10):
    for e4 in range(A4[1],A4[1]-110,-10):
        ints=list(ints0); ints[3]=(s4,e4)
        rows.append(('A4',s4,e4,e4-s4+1,ts.accuracy_direct(xall,ints,gall,vall)))
for s5 in range(A5[0],A5[0]+80,10):
    for e5 in range(A5[1],A5[1]-110,-10):
        ints=list(ints0); ints[4]=(s5,e5)
        rows.append(('A5',s5,e5,e5-s5+1,ts.accuracy_direct(xall,ints,gall,vall)))
R=pd.DataFrame(rows,columns=['amp','s','e','len','acc']); R['loss']=(base-R.acc)*100
R.to_csv(SC+'short_grid.csv',index=False)
for a in ('A4','A5'):
    g=R[R.amp==a]; print(a,'原',ints0[3 if a=='A4' else 4]); print(g[g.loss<=0.15].sort_values('len').head(12).to_string())
