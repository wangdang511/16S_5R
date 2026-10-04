import pickle, json, itertools, numpy as np, pandas as pd, warnings
import primer_design as pdz, iterate_5R as it, compact_primers as cp, tiling_search as ts
warnings.filterwarnings('ignore'); pd.set_option('display.width',260)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
res=pickle.load(open(SC+'pool_opt2.pkl','rb'))
S=pickle.load(open(SC+'S.pkl','rb')); B=pickle.load(open(SC+'big_labeled.pkl','rb'))
tax=B['tax']; ids=B['ids']; A0=B['aln']; gen0=np.array([tax[i][5] for i in ids]); ph0=np.array([tax[i][1] if tax[i][1] else 'unclassified' for i in ids])
rng=np.random.default_rng(3); idx12=np.sort(rng.choice(len(S.ids),12000,replace=False)); rest=np.setdiff1d(np.arange(len(S.ids)),idx12)
sub=lambda ref,ix:pdz.Reference(ref.aln[ix],[ref.ids[i] for i in ix],ref.domain[ix],ref.phylum[ix],ref.col_of_pos,ref.ref_seq,{})
Srest=sub(S,rest)
G=pdz.Reference(A0,ids,np.array(['Bacteria']*len(ids)),ph0,S.col_of_pos,S.ref_seq,{})
gi=np.random.default_rng(1).permutation(len(ids)); Gv=sub(G,gi[len(ids)//2:])
full=lambda ref:np.flatnonzero(((ref.aln[:,19:1510]!=0).all(1))&(ref.domain=='Bacteria'))
VAL={'SILVA':(Srest,full(Srest)),'GG':(Gv,full(Gv))}
# 理想准确率: 全部 5899 条带标签序列
D=pickle.load(open(SC+'search_data.pkl','rb')); xall=D['x']; gall=D['gen']; vall=(pd.Series(gall).map(pd.Series(gall).value_counts()).values>=2)
# 留出集（不在 Gd 里）用于考虑扩增失败的准确率
full0=(A0[:,19:1510]!=0).all(1); idx=np.flatnonzero(full0)
cnt=pd.Series(gen0[idx]).value_counts(); idx=idx[np.isin(gen0[idx],cnt[cnt>=2].index)]
design_ids=set(gi[:len(ids)//2]); held=np.array([i for i in idx if i not in design_ids])
cnt=pd.Series(gen0[held]).value_counts(); held=held[np.isin(gen0[held],cnt[cnt>=2].index)]
Ah=A0[held]; genh=gen0[held]
Gh=pdz.Reference(Ah,[ids[i] for i in held],np.array(['Bacteria']*len(held)),ph0[held],S.col_of_pos,S.ref_seq,{})
xh=np.full(Ah.shape,4,np.int64)
for k,v in {ord('A'):0,ord('C'):1,ord('G'):2,ord('T'):3}.items(): xh[Ah==k]=v
def hit(ref,rows,c):
    M=pdz.site_matrix(ref,c['start'],c['L'])[rows]; ok=(M!=0).all(1)&(M!=45).all(1); h=np.zeros(len(M),bool)
    h[ok]=it.hit_any(M[ok],[p if c['orient']=='F' else pdz.revcomp(p) for p in c['prim']],'right' if c['orient']=='F' else 'left'); return h
def realistic(design,minshared=60):
    n=len(xh); Pres=np.zeros((n,xh.shape[1]),bool)
    for (a,b),h in design: Pres[h,a-1:b]=True
    ident=Pres.any(1); Sm=np.zeros((n,n),np.float32); Lm=np.zeros((n,n),np.float32); cols=np.flatnonzero(Pres.any(0))
    for s0 in range(0,len(cols),300):
        c=cols[s0:s0+300]; Pm=Pres[:,c].astype(np.float32); Lm+=Pm@Pm.T
        X=np.zeros((n,len(c)*4),np.float32)
        for k,p in enumerate(c):
            m=Pres[:,p]&(xh[:,p]<4); X[np.where(m)[0],k*4+xh[m,p]]=1
        Sm+=X@X.T
    rate=np.where(Lm>=minshared,(Lm-Sm)/np.maximum(Lm,1),np.inf); np.fill_diagonal(rate,np.inf)
    mn=rate.min(1,keepdims=True); tied=(rate==mn)&np.isfinite(rate); same=(genh[None,:]==genh[:,None])
    credit=np.where(tied.any(1),(tied&same).sum(1)/np.maximum(tied.sum(1),1),0.0)*ident
    q=pd.Series(genh).map(pd.Series(genh).value_counts()).values>=2
    return float(credit[q].mean()),float(ident[q].mean())
rows=[]; detail={}
for (nm,T),r in res.items():
    if 'error' in r: rows.append(dict(scheme=nm,target=T,error=r['error'])); continue
    ch=r['chosen']; amps=sorted(ch)
    row=dict(scheme=nm,target=T,oligos=sum(c['n'] for f,rr in ch.values() for c in (f,rr)),sites=2*len(amps),
             expansions=sum(c['fold'] for f,rr in ch.values() for c in (f,rr)),severe_dimers=len(r['severe']),hairpins=len(r['hairpins']))
    tms=[cp.oligo_tm(p)[0] for f,rr in ch.values() for c in (f,rr) for p in c['prim']]
    row['Tm_min']=min(tms); row['Tm_max']=max(tms); row['Tm_mean']=float(np.mean(tms))
    spans=[(f['start'],rr['start']+rr['L']-1) for f,rr in (ch[k] for k in amps)]
    row['spans']=' '.join(f'{a}-{b}({b-a+1})' for a,b in spans)
    for vn,(ref,rows_) in VAL.items():
        phs=ref.phylum[rows_]; H=[]; per_site=[]
        for k in amps:
            hf=hit(ref,rows_,ch[k][0]); hr=hit(ref,rows_,ch[k][1]); H.append(hf&hr)
            for h in (hf,hr):
                v=[h[phs==g].mean() for g in pdz.KEY_PHYLA if (phs==g).sum()>=30]; per_site.append((np.mean(v),np.min(v)))
        H=np.array(H).T; cnt_=H.sum(1); allok=H.all(1)
        pa=[min(H[phs==g,i].mean() for g in pdz.KEY_PHYLA if (phs==g).sum()>=30) for i in range(len(amps))]
        pall=[allok[phs==g].mean() for g in pdz.KEY_PHYLA if (phs==g).sum()>=30]
        row[f'{vn}_site_mean_min']=min(x[0] for x in per_site); row[f'{vn}_site_worstphylum_min']=min(x[1] for x in per_site)
        row[f'{vn}_amp_cov']='/'.join(f'{x*100:.0f}' for x in H.mean(0)); row[f'{vn}_amp_worst']='/'.join(f'{x*100:.0f}' for x in pa)
        row[f'{vn}_p_all']=float(allok.mean()); row[f'{vn}_p_ge_half']=float((cnt_>=np.ceil(len(amps)/2)).mean()); row[f'{vn}_p_all_worstphylum']=float(min(pall)); row[f'{vn}_mean_amps']=float(cnt_.mean())
        row[f'{vn}_site_cov']='/'.join(f'{x[0]*100:.0f}' for x in per_site)
    ints=[(f['start']+f['L'],rr['start']-1) for f,rr in (ch[k] for k in amps)]
    row['ideal_acc']=ts.accuracy_direct(xall,ints,gall,vall); row['bases']=sum(b-a+1 for a,b in ints)
    d=[((f['start']+f['L'],rr['start']-1),hit(Gh,np.arange(len(held)),f)&hit(Gh,np.arange(len(held)),rr)) for f,rr in (ch[k] for k in amps)]
    row['realistic_acc'],row['realistic_identified']=realistic(d)
    # 跨区产物：F_i + R_j (j>i)
    cross=[]
    for i,j in itertools.combinations(range(len(amps)),2):
        f=ch[amps[i]][0]; rr=ch[amps[j]][1]; L=rr['start']+rr['L']-1-f['start']+1
        both=hit(Srest,VAL['SILVA'][1],f)&hit(Srest,VAL['SILVA'][1],rr)
        cross.append((L,float(both.mean()),f'A{amps[i]}F+A{amps[j]}R'))
    cross.sort(); row['shortest_cross']=f'{cross[0][2]} {cross[0][0]} bp'; row['shortest_cross_len']=cross[0][0]
    row['cross_pairs_lt600']=' '.join(f'{c[2]}:{c[0]}bp' for c in cross if c[0]<600)
    rows.append(row)
    detail[(nm,T)]=dict(amps=amps,sites=[dict(name=f'A{k}-{o}',orient=o,start=c['start'],L=c['L'],shift=c['shift'],n=c['n'],prim=c['prim'],fold=c['fold'],
        tm=[list(map(lambda x:round(x,1),cp.oligo_tm(p))) for p in c['prim']],gc=[round(pdz.gc_content(p),2) for p in c['prim']],cov_design=[round(x,3) for x in c['cov_design']]) for k in amps for o,c in (('F',ch[k][0]),('R',ch[k][1]))],severe=r['severe'],hairpins=r['hairpins'])
    print(nm,T,{k:(round(v,3) if isinstance(v,float) else v) for k,v in row.items() if k in ('oligos','expansions','ideal_acc','realistic_acc','SILVA_p_all','GG_p_all','Tm_min','Tm_max','severe_dimers')},flush=True)
df=pd.DataFrame(rows); df.to_csv(SC+'compact_eval2.csv',index=False); pickle.dump(detail,open(SC+'compact_detail2.pkl','wb'))
