# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/lvl_ten_valid.csv
import sys; sys.path.insert(0,'/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad')
from par2 import *
SETS=pickle.load(open(SC+'lvl_sets.pkl','rb'))
names=['V1·V2','V3','V4','V5 (A3-F×V5-R)','V6·V7','V8·V9']
rows=[]
for lv,sel in SETS.items():
    for db in ('GG','SILVA'):
        for j,(f,r) in enumerate(PAIRS):
            a=[i for i in byslot[f] if i in sel]; b=[i for i in byslot[r] if i in sel]
            h=np.any([U[i]['ten'][db][0] for i in a],0)&np.any([U[i]['ten'][db][0] for i in b],0)
            v=U[a[0]]['ten'][db][1]&U[b[0]]['ten'][db][1]
            rows.append(dict(level=lv,db=db,amplicon=names[j],cov_valid=float(h[v].mean()),n_valid=int(v.sum())))
R=pd.DataFrame(rows); print(R.pivot_table(index=['db','amplicon'],columns='level',values='cov_valid').round(3).to_string()); R.to_csv(SC+'lvl_ten_valid.csv',index=False)
