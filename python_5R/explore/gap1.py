# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/gap_*
import pickle, json, numpy as np, pandas as pd, warnings, itertools
warnings.filterwarnings('ignore'); pd.set_option('display.width',250)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'ext2.py').read().split("rows=[]")[0]
exec(src)
S=json.load(open(SC+'short_sites.json')); orig=json.load(open(SC+'a5f_final5.json'))
def mk(o,start,L,prims): return dict(orient=o,start=start,L=L,prim=prims,ext=[0]*len(prims))
E1=json.loads(json.dumps(S['D1'])); E1['A5-R']=mk('R',1483,19,['GTTACGACTTMRCCCCART','GTTACGACTTCACCYTCHT'])   # V6V7 缩短 + V8V9 反向 1483
E2=json.loads(json.dumps(S['D1'])); E2['A5-R']=mk('R',1485,19,['TTGTTACGACTTMRCCCCM','TTGTTACGACTTCRCCYTM'])  # 1485
designs={'原（A4 290 / A5 288）':orig,'V6V7 缩短（229 / 288）':S['D1'],'V6V7 + V8 缩短（229 / 185）':S['D2'],'V6V7 缩短 + V8V9 反向 1483（229 / ≈279）':E1,'V6V7 缩短 + V8V9 反向 1485（229 / ≈281）':E2}
READ=150; UDP=10; EFF=READ-UDP
def ends(s):
    o=s['orient']; return [ (s['start']-e) if o=='F' else (s['start']+s['L']-1+e) for e in s['ext']]
rows=[]; info={}
for nm,sites in designs.items():
    sp=[]
    for k in range(1,6):
        F,R=sites[f'A{k}-F'],sites[f'A{k}-R']; fe,re_=ends(F),ends(R)
        Ls=[r-f+1 for f in fe for r in re_]
        f3=F['start']+F['L']-1; r3=R['start']          # 3′ 端（F）/ 引物起点（R）= 内部区边界
        interior=(f3+1,r3-1)
        gaps=[max(0,L-2*EFF) for L in Ls]
        sp.append(dict(amp=k,Lmin=min(Ls),Lmax=max(Ls),gap_min=min(gaps),gap_max=max(gaps),interior=interior,fe=fe,re=re_))
    info[nm]=sp
    # 精度：缺口按最坏组合（最长扩增子）切掉
    ints_full=[a['interior'] for a in sp]
    ints_gap=[]
    for a in sp:
        f5=max(a['fe']) if False else min(a['fe']); r5=max(a['re'])      # 最长扩增子（最坏）
        lo,hi=a['interior']; end1=min(hi,f5+EFF-1); start2=max(lo,r5-EFF+1)
        ints_gap+= [(lo,hi)] if end1>=start2-1 else [(lo,end1),(start2,hi)]
    d=dict(design=nm,ideal_full=round(ts.accuracy_direct(xall,ints_full,gall,vall),4),ideal_gap=round(ts.accuracy_direct(xall,ints_gap,gall,vall),4))
    S2={k:v for k,v in sites.items()}; ev,sev,hp=evaluate((nm,0),S2)
    # 现实准确率（缺口掉位点）
    dd=[]; 
    for k in range(1,6):
        F,R=sites[f'A{k}-F'],sites[f'A{k}-R']; hit=site_hit(Gh,np.arange(len(held)),F)&site_hit(Gh,np.arange(len(held)),R); a=info[nm][k-1]
        f5=min(a['fe']); r5=max(a['re']); lo,hi=a['interior']; end1=min(hi,f5+EFF-1); start2=max(lo,r5-EFF+1)
        if end1>=start2-1: dd.append(((lo,hi),hit))
        else: dd+= [((lo,end1),hit),((start2,hi),hit)]
    for lab,a_,f_ in SET: d['realistic_gap_'+lab]=round(realistic2(dd,a_,f_),4)
    d.update(oligos=ev['oligos'],SILVA_amp=ev['SILVA_amp_cov'],GG_amp=ev['GG_amp_cov'],SILVA_all=round(ev['SILVA_p_all'],3),GG_all=round(ev['GG_p_all'],3),Tm=f"{ev['Tm_min']:.1f}-{ev['Tm_max']:.1f}",severe=ev['severe'])
    rows.append(d); print(d,flush=True)
    for a in sp: print('   A%d 扩增子 %d–%d bp，缺口 %d–%d bp，内部 %s'%(a['amp'],a['Lmin'],a['Lmax'],a['gap_min'],a['gap_max'],a['interior']))
R=pd.DataFrame(rows); R.to_csv(SC+'gap_design_eval.csv',index=False)
pickle.dump(info,open(SC+'gap_info.pkl','wb'))
json.dump({'E1':E1,'E2':E2},open(SC+'gap_sites.json','w'),ensure_ascii=False,indent=1)
