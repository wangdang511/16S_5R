# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/layout6_*
import pickle, sys, json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0,'/home/user/16S_5R/python_5R/explore')
import primer_design as pdz, compact_primers as cp
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
S=pickle.load(open(SC+'S.pkl','rb')); B=pickle.load(open(SC+'big_labeled.pkl','rb'))
ids=B['ids']; tax=B['tax']; ph=np.array([tax[i][1] if tax[i][1] else 'unclassified' for i in ids])
G=pdz.Reference(B['aln'],ids,np.array(['Bacteria']*len(ids)),ph,S.col_of_pos,S.ref_seq,{})
rng=np.random.default_rng(3); idx12=np.sort(rng.choice(len(S.ids),12000,replace=False))
Sb=pdz.Reference(S.aln[idx12],[S.ids[i] for i in idx12],S.domain[idx12],S.phylum[idx12],S.col_of_pos,S.ref_seq,{})
REF=S.ref_seq.upper()
df=pd.read_csv(SC+'lay2_oligos.csv')
def pop(st,n,side_start):
    out=[]
    for ref in (G,Sb):
        M=pdz.site_matrix(ref,st,n)
        if M is None: return None
        ok=(M!=0).all(1)&(M!=ord('-')).all(1)&(ref.domain=='Bacteria'); out.append((M[ok],ref.phylum[ok]))
    return np.vstack([o[0] for o in out]),np.concatenate([o[1] for o in out])
def tmr(seq):
    t=[cp.oligo_tm(e)[0] for e in pdz.expand(seq)]; return min(t),max(t)
def mk(row,d):
    """d>0 向 5′ 外延伸 d nt（取 E. coli 参照碱基）；d<0 从 5′ 端截短。返回 (seq, start, nt)"""
    seq=row.seq; st=int(row.start); o=row.site[-1]; n=len(seq)
    if o=='F':
        if d>=0: ext=REF[st-1-d:st-1]; return ext+seq,st-d,n+d
        return seq[-d:],st-d,n+d
    else:
        if d>=0: ext=pdz.revcomp(REF[st-1+n:st-1+n+d]); return ext+seq,st,n+d
        return seq[-d:],st,n+d
def keycov(M,phs,seqs,o):
    h=np.zeros(len(M),bool); side='right' if o=='F' else 'left'
    for s in seqs:
        top=s if o=='F' else pdz.revcomp(s); h|=pdz._match(M[:,:len(top)] if o=='F' else M[:,-len(top):],top,side) if False else pdz._match(M[:,M.shape[1]-len(top):] if o=='F' else M[:,:len(top)],top,side)
    ks=[h[phs==k].mean() for k in pdz.KEY_PHYLA if (phs==k).sum()>=30]; return float(np.mean(ks)),h

TARGET=60.0
import collections
rows=[]
def maj(col):
    c=collections.Counter(bytes(col).decode()); c.pop('-',None); c.pop('N',None)
    return c.most_common(1)[0][0] if c else 'N'
for site,g in df.groupby('site',sort=False):
    st=int(g.start.iloc[0]); o=site[-1]
    W0=(st-8 if st>8 else st) if o=='F' else st; Wn=((st-W0)+22) if o=='F' else 34
    P=pop(W0,Wn,0); M,phs=P
    keys=[k for k in pdz.KEY_PHYLA if (phs==k).sum()>=30]
    def cov(h): return float(np.mean([h[phs==k].mean() for k in keys]))
    def mwin(n,stt):
        if o=='F': a=stt-W0; return M[:,a:a+n]
        return M[:,0:n]
    def hit(seq,stt):
        top=seq if o=='F' else pdz.revcomp(seq); return pdz._match(mwin(len(seq),stt),top,'right' if o=='F' else 'left')
    def mk2(r,d):
        seq=r.seq; n=len(seq); stt=int(r.start)
        if d<=0:
            if o=='F': return seq[-d:] if d<0 else seq, stt-d, n+d
            return (seq[-d:] if d<0 else seq), stt, n+d
        if o=='F':
            cols=[maj(M[:,(stt-W0)-j]) for j in range(1,d+1)]       # 向左
            return ''.join(reversed(cols))+seq, stt-d, n+d
        cols=[maj(M[:,n+j-1]) for j in range(1,d+1)]               # 向右（正链）
        return pdz.revcomp(''.join(cols))+seq, stt, n+d
    cur={r.name:(r.seq,int(r.start)) for r in g.itertuples()}
    def setcov():
        h=np.zeros(len(M),bool)
        for sq,stt in cur.values(): h|=hit(sq,stt)
        return cov(h)
    base=setcov(); k0=base; chosen={}
    for r in g.itertuples():
        cands=[]
        for d in range(-4,7):
            seq,stt,n=mk2(r,d)
            if n<19 or n>26: continue
            if o=='F' and (stt<W0 or stt+n>W0+Wn): continue
            if o=='R' and n>Wn: continue
            lo,hi=tmr(seq); mid=(lo+hi)/2
            cands.append((abs(mid-TARGET)+0.3*max(0,hi-lo-6),d,seq,stt,n,lo,hi))
        cands.sort(key=lambda x:x[0]); pick=None
        old=cur[r.name]
        for c in cands:
            cur[r.name]=(c[2],c[3])
            if setcov()>=base-0.002: pick=c; break
            cur[r.name]=old
        if pick is None: pick=min(cands,key=lambda x:abs(x[1])); cur[r.name]=(pick[2],pick[3])
        chosen[r.name]=pick
    k1=setcov()
    print(site,'关键门平均覆盖 %.4f → %.4f'%(k0,k1),flush=True)
    for r in g.itertuples():
        b=chosen[r.name]
        rows.append(dict(name=r.name,site=site,kind=r.kind,seq_old=r.seq,seq_new=b[2],d=b[1],nt=b[4],start=b[3],tm_min=round(b[5],1),tm_max=round(b[6],1),exp=len(pdz.expand(b[2])),site_cov_old=round(k0,4),site_cov_new=round(k1,4)))
T=pd.DataFrame(rows); T.to_csv(SC+'lay3_tm.csv',index=False)
print(T[['name','seq_old','seq_new','d','tm_min','tm_max','exp']].to_string(index=False))
print('Tm 范围',T.tm_min.min(),T.tm_max.max(),'展开',T.exp.sum())
