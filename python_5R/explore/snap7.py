# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/snap_pe_compare.csv
import pickle, json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore'); pd.set_option('display.width',250)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
exec(open(SC+'ext2.py').read().split("rows=[]")[0])
def merge(iv):
    iv=sorted((a,b) for a,b in iv if b>=a); out=[]
    for a,b in iv:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
P=pd.read_csv(SC+'snap_positions.csv')
def snap_windows(RL):
    iv=[]
    for _,r in P.drop_duplicates('name').iterrows():
        n=r['name']
        iv.append((r.end+1,r.end+RL) if '_f' in n else (r.start-RL,r.start-1))
    return merge(iv)
VR={'V1':(69,99),'V2':(137,242),'V3':(433,497),'V4':(576,682),'V5':(822,879),'V6':(986,1043),'V7':(1117,1173),'V8':(1243,1294),'V9':(1435,1465)}
def vcov(iv):
    iv=merge(iv); return {k:round(sum(max(0,min(b,y)-max(a,x)+1) for x,y in iv)/(b-a+1),2) for k,(a,b) in VR.items()}
REC=json.load(open(SC+'final_rec_sites.json')); SB,V5R=pickle.load(open(SC+'v5_sb.pkl','rb'))
RL150=140
def ours(sites,extra=()):
    iv=[]
    for k in range(1,6):
        F,R=sites[f'A{k}-F'],sites[f'A{k}-R']; f3=F['start']+F['L']-1; f5=min(F['start']-e for e in F['ext']); r5=max(R['start']+R['L']-1+e for e in R['ext']); rs=R['start']
        iv+=[(f3+1,min(rs-1,f5+RL150-1)),(max(f3+1,r5-RL150+1),rs-1)]
    return merge(iv+list(extra))
rows=[]
items=[('Swift SNAP，PE250（读段 230 nt，估计）',snap_windows(230)),('Swift SNAP，PE250（读段 200 nt，保守估计）',snap_windows(200)),('Swift SNAP，PE150（读段 130 nt）',snap_windows(130)),
       ('我们的推荐设计，PE150（内联 UDP）',ours(REC)),('我们的方案 B（加 V5），PE150（内联 UDP）',ours(SB,[(788,906)]))]
for nm,iv in items:
    a=ts.accuracy_direct(xall,merge(iv),gall,vall); r=dict(scheme=nm,bases=sum(b-x+1 for x,b in merge(iv)),ideal_acc=round(a,4)); r.update(vcov(iv)); rows.append(r)
R=pd.DataFrame(rows); print(R.to_string()); R.to_csv(SC+'snap_pe_compare.csv',index=False)
for RL in (230,130): print(RL,snap_windows(RL))
