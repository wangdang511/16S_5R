# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/layout6_pools_20_28_36.csv
import sys; sys.path.insert(0,'/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad')
from sc2 import *
import json
cands=json.load(open(SC+'a5r_off2/offtarget_primers.json'))
base_sel=({L6[n] for n in T.name if not n.startswith('16S-A5-R') and not n.startswith('16S-A1-F')})|{L6['V1_f']}
r0=evaluate(base_sel|{L6[n] for n in T.name if n.startswith('16S-A5-R')},'L6'); print('当前 A5-R(1174)',{k:round(v,4) for k,v in r0.items() if isinstance(v,float)},flush=True)
out={}
for nm,ol in cands.items():
    st=int(nm[1:5]); ids_=[]
    for o in ol: ids_.append(add(o['name'],'L6A5-R','R',st,o['seq']))
    for i in ids_: byslot.setdefault('L6A5-R',[]).append(i) if i not in byslot.get('L6A5-R',[]) else None
    # 窗口：A5 的 R 起点变化 -> 重建窗口
    W=list(WINS); W[4]=(986,st-1); STRUCT['L6']['wins']=W; s_=STRUCT['L6']; s_['S']=build_S(W); s_['lens']=np.array([b-a+1 for a,b in W],np.float32)
    # 仅保留当前候选的 A5-R
    for j in list(byslot['L6A5-R']): pass
    sel=set(base_sel)|set(ids_)
    r=evaluate(sel,'L6'); out[nm]=r; print(nm,'amp A5 长度',st+18-967+1,{k:round(v,4) for k,v in r.items() if isinstance(v,float)},flush=True)
    byslot['L6A5-R']=[i for i in byslot['L6A5-R'] if i not in ids_]
