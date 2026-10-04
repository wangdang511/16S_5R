# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/v1v2_tm_*.csv/json
import pickle, json, itertools, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore'); pd.set_option('display.width',260)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'robust_realistic.py').read().split("D={'S1'")[0]
exec(src)
rf=open(SC+'refine.py').read()
exec(rf[rf.index('def extend'):rf.index('def site_cov')])
S=pickle.load(open(SC+'S.pkl','rb')); B=pickle.load(open(SC+'big_labeled.pkl','rb'))
tax=B['tax']; ids=B['ids']; ph0=np.array([tax[i][1] if tax[i][1] else 'unclassified' for i in ids])
rng=np.random.default_rng(3); idx12=np.sort(rng.choice(len(S.ids),12000,replace=False))
sub=lambda ref,ix:pdz.Reference(ref.aln[ix],[ref.ids[i] for i in ix],ref.domain[ix],ref.phylum[ix],ref.col_of_pos,ref.ref_seq,{})
S12=sub(S,idx12); G=pdz.Reference(B['aln'],ids,np.array(['Bacteria']*len(ids)),ph0,S.col_of_pos,S.ref_seq,{})
gi=np.random.default_rng(1).permutation(len(ids)); Gd=sub(G,gi[:len(ids)//2])
SET=[('abs60',60,0.0),('abs300',300,0.0),('frac0.8',0,0.8)]
c4=store[('4 扩增子',0.95)]['sites']
print('pool Tm (4-amp 95%):',{n:[round(cp.oligo_tm(p)[0],1) for p in s['prim']] for n,s in c4.items()})
R0=dict(orient='R',start=246,L=17,prim=['TACCYCACCAACWARCT','TACCCCRCCAACTABCT','TACCYTACCAACTARYT'],ext=[0,0,0])
F0=dict(orient='F',start=8,L=17,prim=['AGRGTTTGATYMTGGCT'],ext=[0])
def hit_site(s,ref=None):
    return site_hit(Gh,np.arange(len(held)),s)
def keycov(s):
    out={}
    for vn,(ref,rows_) in VAL.items():
        phs=ref.phylum[rows_]; h=site_hit(ref,rows_,s)
        out[vn]=round(float(np.mean([h[phs==g].mean() for g in pdz.KEY_PHYLA if (phs==g).sum()>=30])),3)
    return out
print('R0 cov',keycov(R0),'F0 cov',keycov(F0),'Tm R0',[round(cp.oligo_tm(p)[0],1) for p in R0['prim']],'F0',[round(cp.oligo_tm(p)[0],1) for p in F0['prim']])
def tm_extend(s,TMIN):
    new=list(s['prim']); es=list(s['ext'])
    for i,p in enumerate(new):
        if cp.oligo_tm(p)[0]<TMIN:
            for e in range(1,7):
                q=extend(s['orient'],s['start'],s['L'],p,e,[Gd,S12])
                new[i],es[i]=q,e
                if cp.oligo_tm(q)[0]>=TMIN: break
    return dict(s,prim=new,ext=es)
res={}
for TMIN in (52,54,56):
    R1=tm_extend(R0,TMIN)
    print(TMIN,'R ext',R1['ext'],R1['prim'],[round(cp.oligo_tm(p)[0],1) for p in R1['prim']],'cov',keycov(R1))
    res[TMIN]=R1
pickle.dump(res,open(SC+'v2_tm.pkl','wb'))
mid=open(SC+'final_eval.py').read()
R54=res[54]; rows=[]
for nm,A in (('V1V2 未均衡 (3条/22)',{'A1-F':F0,'A1-R':R0}),('V1V2 Tm 均衡 (≥54)',{'A1-F':F0,'A1-R':R54})):
    sites={**A,**{f'A{int(k[1])+1}-{k[3]}':v for k,v in c4.items()}}
    ev,sev,hp=evaluate((nm,0),sites)
    d=build(sites); r=dict(design=nm,oligos=ev['oligos'],expansions=ev['expansions'],Tm_min=round(ev['Tm_min'],1),Tm_max=round(ev['Tm_max'],1),severe=ev['severe'],hairpins=ev['hairpins'],
        SILVA_amp=ev['SILVA_amp_cov'],GG_amp=ev['GG_amp_cov'],SILVA_all=round(ev['SILVA_p_all'],3),GG_all=round(ev['GG_p_all'],3),ideal_acc=round(ev['ideal_acc'],4))
    for lab,a,f in SET: r[lab]=round(realistic2(d,a,f),4)
    rows.append(r); print(r,flush=True); print(' severe:',sev,' hairpins:',hp,flush=True)
ev,sev,hp=evaluate(('base',0),c4); print('base4 oligos',ev['oligos'],ev['SILVA_amp_cov'],ev['GG_amp_cov'],round(ev['SILVA_p_all'],3),round(ev['GG_p_all'],3),round(ev['ideal_acc'],4))
pd.DataFrame(rows).to_csv(SC+'v2_tm_design.csv',index=False)
json.dump({'A1-F':F0,'A1-R_unbalanced':R0,'A1-R_balanced':R54},open(SC+'v2_tm_sites.json','w'),ensure_ascii=False,indent=1)
