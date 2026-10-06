# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/v5_* 与 order_list_full_v5.csv
import json, pandas as pd, itertools
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
P=pd.read_csv(SC+'snap_positions.csv')
F={'V1_f':(9,27),'V3_f':(341,357),'V4_f':(517,533),"V6'_f":(1055,1070),'V6_f':(967,985),'V7_f':(1099,1114)}
R={'V2_r':(372,391),'V4_r':(785,803),'V5_r':(907,926),'V8_r':(1390,1407),'V9_r':(1492,1507)}
def combos(F,R):
    rows=[]
    for (fn,(fa,fb)),(rn,(ra,rb)) in itertools.product(F.items(),R.items()):
        if ra>fb: rows.append((fn,rn,rb-fa+1))
    return sorted(rows,key=lambda x:x[2])
for nm,(F_,R_) in {'SNAP':(F,R),'5R':({f'R{i}-F':(a,b) for i,(a,b) in zip(range(1,6),[(103,120),(338,355),(685,702),(944,964),(1175,1193)])},{f'R{i}-R':(a,b) for i,(a,b) in zip(range(1,6),[(314,332),(519,536),(908,927),(1087,1104),(1374,1391)])}),
  '我们':({f'A{i}-F':(a,b) for i,(a,b) in zip(range(1,6),[(8,24),(338,356),(556,576),(906,927),(1223,1241)])},{f'A{i}-R':(a,b) for i,(a,b) in zip(range(1,6),[(246,265),(516,535),(785,804),(1177,1195),(1492,1510)])})}.items():
    c=combos(F_,R_); print(nm,len(c),'种组合；≤300 bp:',[x for x in c if x[2]<=300]); print('   全部:',[f'{a}×{b}:{l}' for a,b,l in c])
