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
for (o,st),tg in TARG.items():
    ps=SITES[(o,st)]; side='right' if o=='F' else 'left'
    Ms=[];phs=[];gs=[]
    for ref,gg in ((G,gen_G),(Sb,gen_S)):
        M=pdz.site_matrix(ref,st,22); ok=(M!=0).all(1)&(M!=ord('-')).all(1)&(ref.domain=='Bacteria')
        Ms.append(M[ok]); phs.append(ref.phylum[ok]); gs.append(gg[ok])
    M=np.vstack(Ms); ph_=np.concatenate(phs); g_=np.concatenate(gs)
    tops=[p if o=='F' else pdz.revcomp(p) for p in ps]
    # 末端 L 与引物长度不同：用引物长度截取
    def hitset(tops_):
        h=np.zeros(len(M),bool)
        for t in tops_:
            h|=pdz._match(M[:,:len(t)],t,side)
        return h
    h0=hitset(tops)
    lab=np.where(np.isin(ph_,pdz.KEY_PHYLA),ph_,'__other__'); w=np.zeros(len(M))
    for gname in np.unique(lab): m=lab==gname; w[m]=1.0/m.sum()
    tmask=np.isin(g_,tg)|((ph_=='Tenericutes')&(('Mycoplasma' in tg)))
    new=[]; h=h0.copy()
    for _ in range(2):
        wi=w*(~h)+ (3.0/max(1,(tmask&~h).sum()))*(tmask&~h)
        if (wi.sum()<=0): break
        top,_=pdz.best_degenerate(M[:,:L],side,2,4,weights=wi)
        if top is None: break
        hh=pdz._match(M[:,:L],top,side)
        if (hh&~h).sum()<5: break
        new.append(top); h|=hh
    prim=[t if o=='F' else pdz.revcomp(t) for t in new]
    keyc=lambda hh:np.mean([hh[ph_==k].mean() for k in pdz.KEY_PHYLA if (ph_==k).sum()>=30])
    out=dict(site=f'{o}{st}',new=prim,nexp=[len(pdz.expand(p)) for p in prim],tm=[round(cp.oligo_tm(p)[0],1) for p in prim],
        key_before=round(keyc(h0),4),key_after=round(keyc(h),4),all_before=round(h0.mean(),4),all_after=round(h.mean(),4),
        ten_before=round(h0[ph_=='Tenericutes'].mean(),3),ten_after=round(h[ph_=='Tenericutes'].mean(),3),
        main_tm=[round(cp.oligo_tm(p)[0],1) for p in ps])
    gen={}
    for g in tg:
        m=(g_==g)
        if m.sum()>=3: gen[g]=(int(m.sum()),round(100*h0[m].mean()),round(100*h[m].mean()))
    out['genus']=gen; res[f'{o}{st}']=out
    print(json.dumps(out,ensure_ascii=False),flush=True)
json.dump(res,open(SC+'sup1.json','w'),ensure_ascii=False)
