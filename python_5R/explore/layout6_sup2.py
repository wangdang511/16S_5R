# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/layout6_*
import pickle, numpy as np, pandas as pd, warnings, collections, sys, json
warnings.filterwarnings('ignore'); pd.set_option('display.width',250)
sys.path.insert(0,'/home/user/16S_5R/python_5R/explore')
import primer_design as pdz, iterate_5R as it, compact_primers as cp
exec(open('/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/loc1.py').read().split("hit={}; dat={}; MAT={}")[0])
rng=np.random.default_rng(3); idx12=np.sort(rng.choice(len(S.ids),12000,replace=False))
Sb=pdz.Reference(S.aln[idx12],[S.ids[i] for i in idx12],S.domain[idx12],S.phylum[idx12],S.col_of_pos,S.ref_seq,{})
gen_G=genus; gen_S=np.array(['']*len(idx12),dtype=object); 
TARG={
('R',246):['Bifidobacterium','Prevotella','Roseburia','Dialister','Blautia','Megasphaera','Veillonella'],
('F',1222):['Neisseria','Oscillospira','Ruminococcus','Faecalibacterium'],
('R',954):['Mycoplasma','Ureaplasma','Veillonella','Streptococcus','Staphylococcus','Lactobacillus','Gemella'],
('F',788):['Streptococcus','Staphylococcus','Lactobacillus','Gemella'],
('F',967):['Streptococcus','Veillonella','Enterococcus','Staphylococcus','Mycoplasma','Gemella'],
('F',314):['Veillonella','Mycoplasma','Peptoniphilus','Haemophilus','Streptococcus','Enterococcus'],
('R',504):['Haemophilus','Streptococcus'],
('R',1174):['Mycoplasma','Veillonella'],
}

res={}
CASES=[('R',246,'Bacteroidetes',['Prevotella'],8,3),('R',954,'Tenericutes',['Mycoplasma','Ureaplasma'],8,3),('R',246,'Actinobacteria',['Bifidobacterium'],8,3)]
for o,st,phy,tg,fold,deg in CASES:
    ps=SITES[(o,st)]; side='right' if o=='F' else 'left'
    Ms=[];phs=[];gs=[]
    for ref,gg in ((G,gen_G),(Sb,gen_S)):
        M=pdz.site_matrix(ref,st,22); ok=(M!=0).all(1)&(M!=ord('-'))[:,:].all(1)&(ref.domain=='Bacteria')
        Ms.append(M[ok]); phs.append(ref.phylum[ok]); gs.append(gg[ok])
    M=np.vstack(Ms); ph_=np.concatenate(phs); g_=np.concatenate(gs)
    tops=[p if o=='F' else pdz.revcomp(p) for p in ps]
    h0=np.zeros(len(M),bool)
    for t in tops: h0|=pdz._match(M[:,:len(t)],t,side)
    tm=(ph_==phy)&((g_=='')|np.isin(g_,tg) if phy=='Tenericutes' else np.isin(g_,tg)|(g_=='') if phy=='Bacteroidetes' else np.isin(g_,tg))
    tmask=(ph_==phy) if phy!='Actinobacteria' else np.isin(g_,tg)
    h=h0.copy(); new=[]
    for _ in range(2):
        wi=(tmask&~h).astype(float)
        if wi.sum()<=0: break
        top,_=pdz.best_degenerate(M[:,:L],side,deg,fold,weights=wi)
        if top is None: break
        hh=pdz._match(M[:,:L],top,side)
        if (hh&~h).sum()<3: break
        new.append(top); h|=hh
    prim=[t if o=='F' else pdz.revcomp(t) for t in new]
    keyc=lambda hh:np.mean([hh[ph_==k].mean() for k in pdz.KEY_PHYLA if (ph_==k).sum()>=30])
    out=dict(case=f'{o}{st}-{phy}',new=prim,nexp=[len(pdz.expand(p)) for p in prim],tm=[round(cp.oligo_tm(p)[0],1) for p in prim],
       target_before=round(h0[tmask].mean(),3),target_after=round(h[tmask].mean(),3),key_before=round(keyc(h0),4),key_after=round(keyc(h),4))
    gen={}
    for g in tg:
        m=(g_==g)
        if m.sum()>=3: gen[g]=(int(m.sum()),round(100*h0[m].mean()),round(100*h[m].mean()))
    out['genus']=gen; print(json.dumps(out,ensure_ascii=False),flush=True)
