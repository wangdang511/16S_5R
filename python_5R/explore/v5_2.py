# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/v5_* 与 order_list_full_v5.csv
import pickle, json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore'); pd.set_option('display.width',250); pd.set_option('display.max_colwidth',60)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
exec(open(SC+'ext2.py').read().split("rows=[]")[0])
REC=json.load(open(SC+'final_rec_sites.json')); D967=json.load(open(SC+'short_sites.json'))['D1']
SNAP=pickle.load(open(SC+'snap_sites.pkl','rb'))
def merge(iv):
    iv=sorted((a,b) for a,b in iv if b>=a); out=[]
    for a,b in iv:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
acc=lambda iv: ts.accuracy_direct(xall,merge(iv),gall,vall)
def keycov(s):
    out={}
    for vn,(ref,rows_) in VAL.items():
        phs=ref.phylum[rows_]; h=site_hit(ref,rows_,s); per=[h[phs==p].mean() for p in pdz.KEY_PHYLA if (phs==p).sum()>=30]; out[vn]=(round(float(np.mean(per)),3),round(float(min(per)),3))
    return out
# 成熟引物位点 vs 我们的
mk=lambda o,start,L,prims: dict(orient=o,start=start,L=L,prim=prims,ext=[0]*len(prims))
V5R=mk('R',907,21,['CCCGTCAATTYCTTTRAGTTT','CCCGTCAATTCCTTTGAGYYY'])
cmp={'515F（SNAP V4_f，517–533，1 条）':SNAP['V4_f'],'我们的 A3-F（556–576，2 条）':REC['A3-F'],'907R 型 V5 反向（907–927，2 条）':V5R,'SNAP V5_r（907–926，1 条）':SNAP['V5_r'],'我们的 A4-F（906–927，1 条）':REC['A4-F'],'967F 型 SNAP V6_f（6 条）':SNAP['V6_f'],'我们的 967 版 A4-F（2 条）':D967['A4-F'],'我们的 A3-R（785–804，1 条）':REC['A3-R'],'SNAP V4_r（785–803，1 条）':SNAP['V4_r']}
rows=[]
for nm,s in cmp.items():
    kc=keycov(s); rows.append(dict(site=nm,oligos=len(s['prim']),expansions=sum(len(pdz.expand(p)) for p in s['prim']),GG=kc['GG'],SILVA=kc['SILVA'],Tm=[round(cp.oligo_tm(p)[0],1) for p in s['prim']]))
print(pd.DataFrame(rows).to_string()); pd.DataFrame(rows).to_csv(SC+'v5_site_compare.csv',index=False)
# 窗口法理想准确率
UDP=10; RL=150-UDP
def win(F,R):   # F: (5′,3′end)  R: (5′end, start)
    f5,f3=F; r5,rs=R; return [(f3+1,min(rs-1,f5+RL-1)),(max(f3+1,r5-RL+1),rs-1)]
ints_cur=[]; 
for k in range(1,6):
    F,R=REC[f'A{k}-F'],REC[f'A{k}-R']; f3=F['start']+F['L']-1; f5=min(F['start']-e for e in F['ext']); r5=max(R['start']+R['L']-1+e for e in R['ext']); ints_cur+=win((f5,f3),(r5,R['start']))
print('当前 REC 窗口',merge(ints_cur),round(acc(ints_cur)*100,2))
# S_B：加 V5R，A4-F→967（A4 变成 229 bp），V5 由 A3-F×V5R 提供
F3=REC['A3-F']; f3_3=F3['start']+F3['L']-1; f5_3=min(F3['start']-e for e in F3['ext'])
v5win=win((f5_3,f3_3),(927,907)); print('V5 产物 A3-F×V5R 窗口',v5win)
ints_B=[]
for k in range(1,6):
    F,R=(D967['A4-F'] if k==4 else REC[f'A{k}-F']),REC[f'A{k}-R']; f3=F['start']+F['L']-1; f5=min(F['start']-e for e in F['ext']); r5=max(R['start']+R['L']-1+e for e in R['ext']); ints_B+=win((f5,f3),(r5,R['start']))
ints_B+=v5win
print('S_B',round(acc(ints_B)*100,2),merge(ints_B))
# S_C：成熟引物位点重组（去掉 V3 反向 A2-R；515F、806R 型、907R 型）
ints_C=[]
for k in (1,4,5):
    F,R=(D967['A4-F'] if k==4 else REC[f'A{k}-F']),REC[f'A{k}-R']; f3=F['start']+F['L']-1; f5=min(F['start']-e for e in F['ext']); r5=max(R['start']+R['L']-1+e for e in R['ext']); ints_C+=win((f5,f3),(r5,R['start']))
A2F=REC['A2-F']; a2f3=A2F['start']+A2F['L']-1; a2f5=min(A2F['start']-e for e in A2F['ext'])
r806=(804,785); ints_C+=win((a2f5,a2f3),r806)       # A2-F × 806R 型：338–804
ints_C+=win((517,533),r806)                          # 515F × 806R：517–804
ints_C+=win((517,533),(927,907))                     # 515F × 907R：517–927
print('S_C',round(acc(ints_C)*100,2),merge(ints_C))
pickle.dump((ints_cur,ints_B,ints_C),open(SC+'v5_ints.pkl','wb'))
