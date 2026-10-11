# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/snap_*
import pickle, json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore'); pd.set_option('display.width',250)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
exec(open(SC+'ext2.py').read().split("rows=[]")[0])
def merge(iv):
    iv=sorted((a,b) for a,b in iv if b>=a); out=[]
    for a,b in iv:
        if out and a<=out[-1][1]+0: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
P=pd.read_csv(SC+'snap_positions.csv'); gname=lambda n:n
RL=130
snap=[]
for _,r in P.drop_duplicates('name').iterrows():
    n=r['name']
    if n.endswith('_f') or '_f' in n: snap.append((r.end+1,r.end+RL))
    else: snap.append((r.start-RL,r.start-1))
snapm=merge(snap); print('SNAP 读段窗口（合并）',snapm)
# 5R
five=[("R1-F","TGGCGAACGGGTGAGTAA","F",103),("R1-R","CCGTGTCTCAGTCCCARTG","R",314),("R2-F","ACTCCTACGGGAGGCAGC","F",338),("R2-R","GTATTACCGCGGCTGCTG","R",519),("R3-F","GTGTAGCGGTGRAATGCG","F",685),("R3-R","CCCGTCAATTCMTTTGAGTT","R",908),("R4-F","GGAGCATGTGGWTTAATTCGA","F",944),("R4-R","CGTTGCGGGACTTAACCC","R",1087),("R5-F","GGAGGAAGGTGGGGATGAC","F",1175),("R5-R","AAGGCCCGGGAACGTATT","R",1374)]
w5=[]; 
for k in range(5):
    f=five[2*k]; r=five[2*k+1]; fe=f[3]+len(f[1])-1; rs=r[3]
    a=(fe+1,min(rs-1,fe+RL)); b=(max(fe+1,rs-RL),rs-1); w5+= [a,b]
w5m=merge(w5); print('5R 窗口',w5m)
# 我们的推荐设计（含 UDP 内联 10 bp 的 2×150）：窗口 140−引物长度
REC=json.load(open(SC+'final_rec_sites.json'))
def ends(s):
    o=s['orient']; return [ (s['start']-e) if o=='F' else (s['start']+s['L']-1+e) for e in s['ext']]
ours=[]; ours_full=[]
for k in range(1,6):
    F,R=REC[f'A{k}-F'],REC[f'A{k}-R']; f3=F['start']+F['L']-1; rs=R['start']; f5=min(ends(F)); r5=max(ends(R))
    ours_full.append((f3+1,rs-1)); ours.append((f3+1,min(rs-1,f5+139))); ours.append((max(f3+1,r5-139),rs-1))
oursm=merge(ours); print('我们的窗口',oursm)
def accs(iv):
    iv=merge(iv); return ts.accuracy_direct(xall,iv,gall,vall), sum(b-a+1 for a,b in iv)
VR={'V1':(69,99),'V2':(137,242),'V3':(433,497),'V4':(576,682),'V5':(822,879),'V6':(986,1043),'V7':(1117,1173),'V8':(1243,1294),'V9':(1435,1465)}
def vcov(iv):
    iv=merge(iv); out={}
    for k,(a,b) in VR.items():
        cov=sum(max(0,min(b,y)-max(a,x)+1) for x,y in iv); out[k]=round(cov/(b-a+1),2)
    return out
rows=[]
for nm,iv in (('Swift SNAP（读段窗口，每读 130 nt）',snap),('5R 原方案（读段窗口，每读 130 nt）',w5),('我们的推荐设计（含 10 bp UDP 内联，2×150）',ours),('我们的推荐设计（完整内部）',ours_full)):
    a,n=accs(iv); r=dict(scheme=nm,bases=n,ideal_acc=round(a,4)); r.update(vcov(iv)); rows.append(r)
R=pd.DataFrame(rows); print(R.to_string()); R.to_csv(SC+'snap_compare.csv',index=False)
json.dump({'snap':snapm,'five':w5m,'ours':oursm,'ours_full':merge(ours_full)},open(SC+'snap_windows.json','w'))
