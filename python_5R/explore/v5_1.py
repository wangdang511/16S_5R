# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/v5_* 与 order_list_full_v5.csv
import pickle, json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore'); pd.set_option('display.width',250)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
exec(open(SC+'ext2.py').read().split("rows=[]")[0])
REC=json.load(open(SC+'final_rec_sites.json'))
def merge(iv):
    iv=sorted((a,b) for a,b in iv if b>=a); out=[]
    for a,b in iv:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
ints=[(REC[f'A{k}-F']['start']+REC[f'A{k}-F']['L'],REC[f'A{k}-R']['start']-1) for k in range(1,6)]
base=ts.accuracy_direct(xall,merge(ints),gall,vall); print('当前（完整内部）',round(base,4),ints)
# 在当前设计上加不同的 V5 区间
tests={'V5 核心 822–879':[(822,879)],'V5 区 805–905':[(805,905)],'785–906（A3-R 到 A4-F 之间全部）':[(785,906)],'V5 窗口（R 在 907 往前读 130）777–906':[(777,906)],'V5 窗口（R 在 905 往前读 120）785–885':[(785,885)],'V5 窗口 822–906':[(822,906)]}
for nm,iv in tests.items():
    a=ts.accuracy_direct(xall,merge(ints+iv),gall,vall); print(f'{nm}: {a*100:.2f}%  (+{(a-base)*100:.2f})')
# 去掉 A3 的 R 引物后：A3 变成 556–905/926 的双端读段，读段窗口 577–696（F 往后 120）与 R 往前读
for nm,(F,Rp) in {'A3\'=556–926，读段 F 往后 120 + R 往前 120（含 UDP）':((577,696),(787,906)),'A3\'=556–905，R 5′ 905':((577,696),(786,885))}.items():
    ints2=[i for k,i in enumerate(ints) if k!=2]+[F,Rp]
    a=ts.accuracy_direct(xall,merge(ints2),gall,vall); print(nm,f'{a*100:.2f}%  ({(a-base)*100:+.2f})','总碱基',sum(b-a_+1 for a_,b in merge(ints2)))
