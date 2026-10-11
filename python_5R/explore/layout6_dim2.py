# 注：此脚本在本次会话中从 scratchpad 运行（借用 wangdang511/Olivar_primer 的 SADDLE 打分 PrimerSetBadnessFast），依赖未入库的中间文件；结果见 docs/primer_design/layout6_dimer_*.csv
import sys, numpy as np, pandas as pd, collections
sys.path.insert(0,'/home/user/16S_5R/python_5R')
import primer_design as pdz
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
E=pd.read_pickle(SC+'dim_E.pkl').reset_index(drop=True)
rc=lambda s: pdz.revcomp(s.upper()).lower()
GC=set('cg'); END={4:1,5:4,6:20}; MID={7:100,8:500}
# 复刻 Olivar 的打分，但按“伙伴引物（q）”分解：p 的反向互补窗口 c[j:j+k] 与 q 的末端/中间 k-mer 相同
seqs=list(E.seq); oli=list(E.oligo); conc=list(E.conc_rel)
n=len(seqs)
endidx={k:collections.defaultdict(list) for k in END}; midx={k:collections.defaultdict(list) for k in MID}
for qi,q in enumerate(seqs):
    l=len(q)
    for k in END: endidx[k][q[-k:]].append(qi)
    for k in MID:
        for j in range(l-k+1): midx[k][q[j:j+k]].append((qi,j,l))
pair=collections.defaultdict(float); detail=collections.defaultdict(lambda: collections.defaultdict(float))
for pi,p in enumerate(seqs):
    c=rc(p); l=len(c)
    for k,w in END.items():
        for j in range(l-k+1):
            km=c[j:j+k]
            for qi in endidx[k].get(km,[]):
                val=conc[qi]*w/(j+1)*2**sum(b in GC for b in km)*conc[pi]
                pair[(pi,qi)]+=val; detail[(pi,qi)][f'end{k}:{km}@{j}']+=val
    for k,w in MID.items():
        for j in range(l-k+1):
            km=c[j:j+k]
            for qi,jj,ql in midx[k].get(km,[]):
                val=conc[qi]*w*(1/(ql-jj-k+1))/(j+1)*2**sum(b in GC for b in km)*conc[pi]
                pair[(pi,qi)]+=val; detail[(pi,qi)][f'mid{k}:{km}@{j}/{jj}']+=val
# 校验：每条 p 的总和应与 Olivar 逐条坏度同量级（Olivar 本身对 p 的 conc 乘一次、对 hash 里的 q 乘 conc）
tot=sum(pair.values()); print('分解后总坏度',round(tot),'(Olivar 加权总坏度 171002)')
# 汇总到寡核苷酸对（无序）
og=collections.defaultdict(float)
for (pi,qi),v in pair.items(): og[tuple(sorted((oli[pi],oli[qi])))]+=v
top=sorted(og.items(),key=lambda x:-x[1])[:25]
print('坏度最高的寡核苷酸对（无序，含自身）:')
for (a,b),v in top: print(f'  {a:14s} × {b:14s} {v:10.0f}  占总 {100*v/tot:4.1f}%')
# 每个高坏度寡核苷酸的主要来源 k-mer
bad=E.groupby('oligo').bad_w.mean().sort_values(ascending=False)
print('\n坏度最高的寡核苷酸及主要原因 k-mer:')
for o in bad.index[:8]:
    ids=[i for i in range(n) if oli[i]==o]; agg=collections.defaultdict(float); partners=collections.defaultdict(float)
    for pi in ids:
        for (a,b),d in detail.items():
            if a!=pi: continue
            for kk,v in d.items(): agg[(kk,oli[b])]+=v
    t=sorted(agg.items(),key=lambda x:-x[1])[:3]
    print(f'  {o} 均值 {bad[o]:.0f}:',[(kk,part,round(v)) for (kk,part),v in t])
import pickle; pickle.dump((og,),open(SC+'dim_og.pkl','wb'))
