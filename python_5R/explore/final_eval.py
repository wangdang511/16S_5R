import pickle, json, itertools, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore'); pd.set_option('display.width',280); pd.set_option('display.max_colwidth',60)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'eval_final.py').read().split('rows=[]; detail={}')[0]
exec(src)
RF=pickle.load(open(SC+'refine.pkl','rb')); POOL=pickle.load(open(SC+'pool_opt2.pkl','rb'))
def final_sites(key):
    ch=POOL[key]['chosen']; sites={}
    for k,(f,rr) in ch.items():
        for o,c in (('F',f),('R',rr)):
            r=RF[key][f'A{k}-{o}']; prims=list(r['before']); es=[0]*len(prims)
            for i,(p,e) in enumerate(zip(r['after'],r['ext'])):
                # 只接受 A3-F 的延长（覆盖率基本不变）；A5-F 的延长使覆盖率下降 5 个百分点，不采用
                if e>0 and f'A{k}-{o}'=='A3-F': prims[i]=p; es[i]=e
            sites[f'A{k}-{o}']=dict(orient=o,start=c['start'],L=c['L'],prim=prims,ext=es)
    return sites
def window(o,start,L,e): return (start-e,L+e) if o=='F' else (start,L+e)
def ohit(ref,rows,o,start,L,top):
    M=pdz.site_matrix(ref,start,L)[rows]; ok=(M!=0).all(1)&(M!=45).all(1); h=np.zeros(len(M),bool)
    h[ok]=pdz._match(M[ok],top,'right' if o=='F' else 'left'); return h
top_of=lambda o,p: p if o=='F' else pdz.revcomp(p)
def site_hit(ref,rows,s):
    h=np.zeros(len(rows),bool)
    for p,e in zip(s['prim'],s['ext']):
        st,Ln=window(s['orient'],s['start'],s['L'],e); h|=ohit(ref,rows,s['orient'],st,Ln,top_of(s['orient'],p))
    return h
def evaluate(key,sites):
    amps=sorted(set(int(n[1]) for n in sites)); row=dict(scheme=key[0],target=key[1])
    pool={n:s['prim'] for n,s in sites.items()}
    sev,hp=cp.pool_issues(pool)
    row['oligos']=sum(len(s['prim']) for s in sites.values()); row['expansions']=sum(len(pdz.expand(p)) for s in sites.values() for p in s['prim'])
    tms=[cp.oligo_tm(p)[0] for s in sites.values() for p in s['prim']]; row['Tm_min']=min(tms); row['Tm_max']=max(tms)
    row['severe']=len(sev); row['hairpins']=len(hp)
    spans=[(min(window(sites[f'A{k}-F']['orient'],sites[f'A{k}-F']['start'],sites[f'A{k}-F']['L'],max(sites[f'A{k}-F']['ext']))[0],sites[f'A{k}-F']['start']),
            sites[f'A{k}-R']['start']+sites[f'A{k}-R']['L']+max(sites[f'A{k}-R']['ext'])-1) for k in amps]
    row['spans']=' '.join(f'{a}-{b}({b-a+1})' for a,b in spans)
    for vn,(ref,rows_) in VAL.items():
        phs=ref.phylum[rows_]; H=[]; site_cov={}
        for k in amps:
            hf=site_hit(ref,rows_,sites[f'A{k}-F']); hr=site_hit(ref,rows_,sites[f'A{k}-R']); H.append(hf&hr)
            for nm,h in ((f'A{k}-F',hf),(f'A{k}-R',hr)):
                v=[h[phs==g].mean() for g in pdz.KEY_PHYLA if (phs==g).sum()>=30]; site_cov[nm]=(float(np.mean(v)),float(np.min(v)))
        H=np.array(H).T; allok=H.all(1); cnt=H.sum(1)
        pa=[min(H[phs==g,i].mean() for g in pdz.KEY_PHYLA if (phs==g).sum()>=30) for i in range(len(amps))]
        pall=[allok[phs==g].mean() for g in pdz.KEY_PHYLA if (phs==g).sum()>=30]
        row[f'{vn}_amp_cov']='/'.join(f'{x*100:.0f}' for x in H.mean(0)); row[f'{vn}_amp_worst']='/'.join(f'{x*100:.0f}' for x in pa)
        row[f'{vn}_p_all']=float(allok.mean()); row[f'{vn}_p_all_worst']=float(min(pall)); row[f'{vn}_p_ge_half']=float((cnt>=np.ceil(len(amps)/2)).mean()); row[f'{vn}_mean_amps']=float(cnt.mean())
        row[f'{vn}_site_cov']=site_cov
    ints=[(sites[f'A{k}-F']['start']+sites[f'A{k}-F']['L'],sites[f'A{k}-R']['start']-1) for k in amps]
    row['ideal_acc']=ts.accuracy_direct(xall,ints,gall,vall); row['bases']=sum(b-a+1 for a,b in ints)
    d=[((a,b),site_hit(Gh,np.arange(len(held)),sites[f'A{k}-F'])&site_hit(Gh,np.arange(len(held)),sites[f'A{k}-R'])) for (a,b),k in zip(ints,amps)]
    row['realistic_acc'],row['realistic_identified']=realistic(d)
    cross=[]
    rS=VAL['SILVA']
    for i,j in itertools.combinations(range(len(amps)),2):
        f=sites[f'A{amps[i]}-F']; r=sites[f'A{amps[j]}-R']; L=r['start']+r['L']-1-f['start']+1
        both=site_hit(rS[0],rS[1],f)&site_hit(rS[0],rS[1],r); cross.append((L,float(both.mean()),f'A{amps[i]}F+A{amps[j]}R'))
    cross.sort(); row['cross']=' '.join(f'{c[2]}:{c[0]}bp({c[1]*100:.0f}%)' for c in cross[:3])
    return row,sev,hp
rows=[]; store={}
for key in [('4 扩增子',0.95),('3 扩增子',0.95)]:
    sites=final_sites(key); row,sev,hp=evaluate(key,sites); rows.append(row); store[key]=dict(sites=sites,severe=sev,hairpins=hp)
    print(key,{k:(round(v,3) if isinstance(v,float) else v) for k,v in row.items() if not k.endswith('site_cov')},'severe',sev,'hairpins',hp,flush=True)
# 90% 版本与 5 扩增子（用原优化结果，未做延长）
for key in [('4 扩增子',0.9),('3 扩增子',0.9),('5 扩增子',0.9)]:
    ch=POOL[key]['chosen']; sites={}
    for k,(f,rr) in ch.items():
        for o,c in (('F',f),('R',rr)): sites[f'A{k}-{o}']=dict(orient=o,start=c['start'],L=c['L'],prim=list(c['prim']),ext=[0]*len(c['prim']))
    row,sev,hp=evaluate(key,sites); rows.append(row); store[key]=dict(sites=sites,severe=sev,hairpins=hp)
    print(key,{k:(round(v,3) if isinstance(v,float) else v) for k,v in row.items() if not k.endswith('site_cov')},flush=True)
pd.DataFrame(rows).to_pickle(SC+'final_eval_df.pkl'); pickle.dump(store,open(SC+'final_store.pkl','wb'))
