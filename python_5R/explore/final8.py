# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/order_list_full.csv 与 final_*
import pickle, collections, json, sys, numpy as np, pandas as pd
sys.path.insert(0,'.')
import explore.offtarget as ot, primer_design as pdz, compact_primers as cp
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
REC=json.load(open(SC+'final_rec_sites.json')); sites_short,sites_long,sup=pickle.load(open(SC+'final1.pkl','rb')); rec,SUPSEQ,CH=pickle.load(open(SC+'final5.pkl','rb'))
# --- 整池脱靶（推荐池）
ol=json.load(open(SC+'fin3/offtarget_primers.json'))['REC']; res,allh=pickle.load(open(SC+'fin3/offtarget_hits_REC.pkl','rb')); H=[h for h in allh['REC'] if '_' not in h[0]]
summ=[]
for mm in (1,2,3):
    P=[p for p in ot.pairs(H,1500,mm) if p[3]>=80]; uniq=set((p[0],p[1]//50,p[2]//50) for p in P); loc=set((h[0],h[1]//20,h[2]) for h in H if h[4]<=mm); summ.append(dict(mm=mm,loci=len(loc),products=len(uniq)))
pd.DataFrame(summ).to_csv(SC+'final_offtarget_summary.csv',index=False); print(pd.DataFrame(summ))
print('rCRS≤4',sorted({(h[1],h[2],h[3].split('#')[0],h[4]) for h in res[('REC','rCRS')] if h[4]<=4}))
def per(hitsd,names):
    out={}
    byo=collections.defaultdict(set)
    for h in hitsd:
        if h[4]<=2: byo[h[3].split('#')[0]].add((h[0],h[1]//20,h[2]))
    for n in names: out[n]=len(byo[n])
    return out
recper=per(H,[o['name'] for o in ol])
res2,allh2=pickle.load(open(SC+'fin/offtarget_hits_ALL.pkl','rb')); H2=[h for h in allh2['ALL'] if '_' not in h[0]]
allo=json.load(open(SC+'final_all_oligos.json')); altper=per(H2,list(allo))
# --- 订购清单行
REGION={'A1':'V1·V2','A2':'V3','A3':'V4','A4':'V6·V7','A5':'V8·V9'}
def pos(s,e_or_len,p,sup=False):
    o=s['orient']; e=e_or_len
    return (s['start']-e,s['start']+s['L']-1) if o=='F' else (s['start'],s['start']+s['L']-1+e)
def stats(p):
    ts=[cp.oligo_tm(x)[0] for x in pdz.expand(p)]; gc=sum(c in 'GCS' for c in p)/len(p)
    return round(cp.oligo_tm(p)[0],1),round(min(ts),1),round(max(ts),1),round(gc,2),len(pdz.expand(p))
rows=[]
SUPTAG={'A1-R':'可选','A2-F':'建议','A4-R':'建议','A5-F':'建议'}
for sn in ['A1-F','A1-R','A2-F','A2-R','A3-F','A3-R','A4-F','A4-R','A5-F','A5-R']:
    s=REC[sn]; core=sites_long[sn]
    n_core=len(core['prim'])
    for i,(p,e) in enumerate(zip(s['prim'],s['ext'])):
        issup=i>=n_core; a,b=pos(s,e,p)
        t,lo,hi,gc,ne=stats(p); nm=f"SMURF5-{sn}"+('' if len(s['prim'])==1 else f".{i+1}")+("（补充）" if issup else "")
        loci=recper[f'{sn}_{i}']
        rows.append(dict(order=nm,site=sn,region=REGION[sn[:2]],orient=sn[3],group='支原体补充' if issup else '核心',recommend=(SUPTAG[sn] if issup else '必订'),seq=p,nt=len(p),pos=f'{a}–{b}',Tm=t,Tm_lo=lo,Tm_hi=hi,GC=gc,expansions=ne,human_loci_le2=loci,human_per_exp=round(loci/ne,1),is_sup=issup))
# 备选
alt=[('A2-R.alt短版','A2-R','V3 短版反向引物（505–523，扩增子 338–523 = 186 bp）',sites_short['A2-R']['prim'][0],505,523),
     ('A4-F.alt967.1','A4-F','V6·V7 缩短版正向引物（扩增子 967–1195 = 229 bp，需 .1 + .2 两条）',allo['A4-F.alt967.1'],967,985),
     ('A4-F.alt967.2','A4-F','V6·V7 缩短版正向引物（同上）',allo['A4-F.alt967.2'],967,985)]
for k,sn,desc,p,a,b in alt:
    t,lo,hi,gc,ne=stats(p); loci=altper[k]
    rows.append(dict(order='SMURF5-'+k,site=sn,region=REGION[sn[:2]],orient=sn[3],group='备选',recommend=desc,seq=p,nt=len(p),pos=f'{a}–{b}',Tm=t,Tm_lo=lo,Tm_hi=hi,GC=gc,expansions=ne,human_loci_le2=loci,human_per_exp=round(loci/ne,1),is_sup=False))
D=pd.DataFrame(rows)
# 混合比例：每个位点总量 2.5 µL（100 µM 母液 / 100 µL 10× 混合液）；补充占 10%，核心按位点内覆盖比例分
core_share={'A1-R':[0.45,0.35,0.20],'A3-F':[0.9,0.1],'A5-F':[0.5,0.5]}
def shares(sn,with_sup):
    g=D[(D.site==sn)&(D.group!='备选')]; core=g[g.group=='核心']; sup=g[g.group=='支原体补充']
    base=core_share.get(sn,[1.0]*1); 
    if len(core)==1: base=[1.0]
    out={}
    sf=0.10 if (with_sup and len(sup)) else 0.0
    for (idx,r),b in zip(core.iterrows(),base): out[idx]=round(b*(1-sf),3)
    for idx,r in sup.iterrows(): out[idx]=round(0.10 if with_sup else 0.0,3)
    return out
D['share_with_sup']=np.nan; D['share_no_sup']=np.nan
for sn in D.site.unique():
    for idx,v in shares(sn,True).items(): D.loc[idx,'share_with_sup']=v
    for idx,v in shares(sn,False).items(): D.loc[idx,'share_no_sup']=v
D['uL_with_sup']=(D.share_with_sup*2.5).round(2); D['uL_no_sup']=(D.share_no_sup*2.5).round(2)
D.to_csv(SC+'order_list_full.csv',index=False)
pd.set_option('display.width',250); pd.set_option('display.max_colwidth',40)
print(D[['order','recommend','seq','nt','pos','Tm','expansions','human_per_exp','share_with_sup','uL_with_sup','uL_no_sup']].to_string())
print('总母液 µL（含补充）',D.uL_with_sup.sum(),'不含',D.uL_no_sup.sum())
