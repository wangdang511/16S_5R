import pickle, json, itertools, numpy as np, pandas as pd, warnings
import primer_design as pdz, iterate_5R as it, compact_primers as cp, tiling_search as ts
warnings.filterwarnings('ignore'); pd.set_option('display.width',260); pd.set_option('display.max_colwidth',50)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
exec(open(SC+'eval_final.py').read().split('rows=[]; detail={}')[0].split("res=pickle.load")[0])  # imports/consts
res=pickle.load(open(SC+'pool_opt2.pkl','rb'))
S=pickle.load(open(SC+'S.pkl','rb')); B=pickle.load(open(SC+'big_labeled.pkl','rb'))
tax=B['tax']; ids=B['ids']; A0=B['aln']; gen0=np.array([tax[i][5] for i in ids]); ph0=np.array([tax[i][1] if tax[i][1] else 'unclassified' for i in ids])
rng=np.random.default_rng(3); idx12=np.sort(rng.choice(len(S.ids),12000,replace=False)); rest=np.setdiff1d(np.arange(len(S.ids)),idx12)
sub=lambda ref,ix:pdz.Reference(ref.aln[ix],[ref.ids[i] for i in ix],ref.domain[ix],ref.phylum[ix],ref.col_of_pos,ref.ref_seq,{})
S12,Srest=sub(S,idx12),sub(S,rest)
G=pdz.Reference(A0,ids,np.array(['Bacteria']*len(ids)),ph0,S.col_of_pos,S.ref_seq,{})
gi=np.random.default_rng(1).permutation(len(ids)); Gd,Gv=sub(G,gi[:len(ids)//2]),sub(G,gi[len(ids)//2:])
full=lambda ref:np.flatnonzero(((ref.aln[:,19:1510]!=0).all(1))&(ref.domain=='Bacteria'))
VAL={'SILVA':(Srest,full(Srest)),'GG':(Gv,full(Gv))}
def window(o,start,L,e): return (start-e,L+e) if o=='F' else (start,L+e)
def ohit(ref,rows,o,start,L,top):
    M=pdz.site_matrix(ref,start,L)[rows]; ok=(M!=0).all(1)&(M!=45).all(1); h=np.zeros(len(M),bool)
    h[ok]=pdz._match(M[ok],top,'right' if o=='F' else 'left'); return h
def top_of(o,prim): return prim if o=='F' else pdz.revcomp(prim)
def extend(o,start,L,prim,e,refs):
    """沿模板把 5′ 端延长 e 个碱基，新增位置取该寡核苷酸所覆盖序列的多数碱基"""
    s0,L0=window(o,start,L,0); rows=[(r,np.arange(len(r.aln))) for r in refs]
    base=[];
    for r,_ in rows:
        M=pdz.site_matrix(r,s0,L0); ok=(M!=0).all(1)&(M!=45).all(1)&(r.domain=='Bacteria'); h=np.zeros(len(M),bool); h[ok]=pdz._match(M[ok],top_of(o,prim),'right' if o=='F' else 'left')
        base.append((r,h))
    add=''
    cols=range(start-e,start) if o=='F' else range(start+L,start+L+e)
    for p in cols:
        cnt=np.zeros(4)
        for r,h in base:
            col=r.aln[h,p-1]
            for k,b in enumerate(b'ACGT'): cnt[k]+=(col==b).sum()
        add+= 'ACGT'[int(cnt.argmax())]
    t=top_of(o,prim)
    newtop=(add+t) if o=='F' else (t+add)
    return newtop if o=='F' else pdz.revcomp(newtop)   # 返回 5'->3' 引物
def site_cov(o,start,L,prims,e_list,vname):
    ref,rows=VAL[vname]; phs=ref.phylum[rows]; h=np.zeros(len(rows),bool)
    for p,e in zip(prims,e_list):
        s,Ln=window(o,start,L,e); h|=ohit(ref,rows,o,s,Ln,top_of(o,p) if e==0 else top_of(o,p))
    per=[h[phs==g].mean() for g in pdz.KEY_PHYLA if (phs==g).sum()>=30]; return float(np.mean(per)),float(np.min(per)),h
TMIN,TMAX=56.0,64.0
out={}
for key in [('4 扩增子',0.95),('3 扩增子',0.95)]:
    ch=res[key]['chosen']; sites={}
    for k,(f,rr) in ch.items():
        for o,c in (('F',f),('R',rr)):
            prims=list(c['prim']); es=[0]*len(prims); new=list(prims)
            for i,p in enumerate(prims):
                if cp.oligo_tm(p)[0]<TMIN:
                    best=(p,0)
                    for e in range(1,7):
                        q=extend(o,c['start'],c['L'],p,e,[Gd,S12])   # 只用设计集的序列确定新增碱基
                        if cp.oligo_tm(q)[0]>=TMIN: best=(q,e); break
                    new[i],es[i]=best
            # 覆盖率是否掉
            before={v:site_cov(o,c['start'],c['L'],prims,[0]*len(prims),v)[:2] for v in VAL}
            def cov_ext(newp,es_):
                res_={}
                for v in VAL:
                    ref,rows=VAL[v]; phs=ref.phylum[rows]; h=np.zeros(len(rows),bool)
                    for p,e in zip(newp,es_):
                        s,Ln=window(o,c['start'],c['L'],e); h|=ohit(ref,rows,o,s,Ln,top_of(o,p))
                    per=[h[phs==g].mean() for g in pdz.KEY_PHYLA if (phs==g).sum()>=30]; res_[v]=(float(np.mean(per)),float(np.min(per)))
                return res_
            after=cov_ext(new,es)
            sites[f'A{k}-{o}']=dict(orient=o,start=c['start'],L=c['L'],before=prims,after=new,ext=es,cov_before=before,cov_after=after,tm_before=[round(cp.oligo_tm(p)[0],1) for p in prims],tm_after=[round(cp.oligo_tm(p)[0],1) for p in new])
    out[key]=sites
    print('==',key)
    for nm,s in sites.items():
        if any(s['ext']): print(nm,s['ext'],s['tm_before'],'->',s['tm_after'],'cov',{v:(round(s['cov_before'][v][0],3),round(s['cov_after'][v][0],3)) for v in VAL},'min-phylum',{v:(round(s['cov_before'][v][1],2),round(s['cov_after'][v][1],2)) for v in VAL},s['after'])
pickle.dump(out,open(SC+'refine.pkl','wb'))
