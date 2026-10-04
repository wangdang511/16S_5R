import pickle, json, itertools, numpy as np, pandas as pd, warnings
import primer_design as pdz, iterate_5R as it, tiling_search as ts
warnings.filterwarnings('ignore'); pd.set_option('display.width',250)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
S=pickle.load(open(SC+'S.pkl','rb')); B=pickle.load(open(SC+'big_labeled.pkl','rb')); sets=json.load(open(SC+'final5b_sets.json'))
tax=B['tax']; ids=B['ids']; ph=np.array([tax[i][1] if tax[i][1] else 'unclassified' for i in ids])
G=pdz.Reference(B['aln'],ids,np.array(['Bacteria']*len(ids)),ph,S.col_of_pos,S.ref_seq,{})
rng=np.random.default_rng(3); idx12=np.sort(rng.choice(len(S.ids),12000,replace=False)); rest=np.setdiff1d(np.arange(len(S.ids)),idx12)
sub=lambda ref,ix:pdz.Reference(ref.aln[ix],[ref.ids[i] for i in ix],ref.domain[ix],ref.phylum[ix],ref.col_of_pos,ref.ref_seq,{})
Srest=sub(S,rest); gi=np.random.default_rng(1).permutation(len(ids)); Gv=sub(G,gi[len(ids)//2:])
full=lambda ref:np.flatnonzero(((ref.aln[:,19:1510]!=0).all(1))&(ref.domain=='Bacteria'))
names={1:'V1V2',2:'V3',3:'V4',4:'V6V7',5:'V8V9'}
def hits(ref,rows,n,which):
    d=sets[n]; M=pdz.site_matrix(ref,d['st'],d['L'])[rows]; ok=(M!=0).all(1)&(M!=45).all(1)
    tops=[p if d['o']=='F' else pdz.revcomp(p) for p in d[which]]; h=np.zeros(len(M),bool); h[ok]=it.hit_any(M[ok],tops,'right' if d['o']=='F' else 'left'); return h
rows=[]
for vn,ref in [('SILVA',Srest),('GG',Gv)]:
    rr=full(ref); phs=ref.phylum[rr]
    for which in ['single','multi']:
        H={k:hits(ref,rr,f'A{k}-F',which)&hits(ref,rr,f'A{k}-R',which) for k in range(1,6)}
        for r in range(1,6):
            for comb in itertools.combinations(range(1,6),r):
                ok=np.all([H[k] for k in comb],axis=0); cnt=np.sum([H[k] for k in comb],axis=0)
                per=[ok[phs==g].mean() for g in pdz.KEY_PHYLA if (phs==g).sum()>=30]
                rows.append(dict(db=vn,primers=which,combo='+'.join(names[k] for k in comb),n=r,p_all=ok.mean(),p_ge_half=(cnt>=np.ceil(r/2)).mean(),p_ge1=(cnt>=1).mean(),worst_phylum_all=min(per)))
df=pd.DataFrame(rows); df.to_csv(SC+'subset_amp.csv',index=False)
a=pd.read_csv(SC+'subset_acc.csv')
m=df[(df.primers=='multi')].pivot_table(index='combo',columns='db',values=['p_all','p_ge1','worst_phylum_all']).round(3)
m.columns=[f'{x}_{y}' for x,y in m.columns]
j=a.set_index('combo').join(m)
print(j.sort_values(['n','acc_all'],ascending=[True,False]).round(3).to_string())
j.to_csv(SC+'subset_joined.csv')
# 5R baseline accuracy on all labelled sequences
D=pickle.load(open(SC+'search_data.pkl','rb')); x=D['x']; gen=D['gen']; va=(pd.Series(gen).map(pd.Series(gen).value_counts()).values>=2)
FIVE=[(103,18,314,19),(338,18,519,18),(685,18,908,20),(944,21,1087,18),(1175,19,1374,18)]
iv5=[(fs+fl,rs-1) for fs,fl,rs,rl in FIVE]
print('5R current on all labelled:',ts.accuracy_direct(x,iv5,gen,va),' full 20-1510:',ts.accuracy_direct(x,[(20,1510)],gen,va))
