import pickle, json, sys, time, numpy as np, pandas as pd, warnings
from multiprocessing import Pool
import primer_design as pdz, iterate_5R as it, compact_primers as cp
warnings.filterwarnings('ignore')
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
S=pickle.load(open(SC+'S.pkl','rb')); B=pickle.load(open(SC+'big_labeled.pkl','rb')); sets=json.load(open(SC+'final5b_sets.json'))
tax=B['tax']; ids=B['ids']; ph=np.array([tax[i][1] if tax[i][1] else 'unclassified' for i in ids])
G=pdz.Reference(B['aln'],ids,np.array(['Bacteria']*len(ids)),ph,S.col_of_pos,S.ref_seq,{})
rng=np.random.default_rng(3); idx12=np.sort(rng.choice(len(S.ids),12000,replace=False))
sub=lambda ref,ix:pdz.Reference(ref.aln[ix],[ref.ids[i] for i in ix],ref.domain[ix],ref.phylum[ix],ref.col_of_pos,ref.ref_seq,{})
S12=sub(S,idx12); gi=np.random.default_rng(1).permutation(len(ids)); Gd=sub(G,gi[:len(ids)//2])
SITES=[(f'A{k}-{o}',o,sets[f'A{k}-{o}']['st'],sets[f'A{k}-{o}']['L']) for k in (1,2,3,4,5) for o in 'FR']
def work(args):
    site,T=args; t=time.time()
    c=cp.enumerate_configs(site,[Gd,S12],T,shifts=range(-2,3),lengths=(16,17,18,19,20,22,24))
    return site[0],T,c,time.time()-t
if __name__=='__main__':
    tasks=[(s,T) for T in (0.90,0.95) for s in SITES]
    out={}
    with Pool(4) as p:
        for name,T,c,dt in p.imap_unordered(work,tasks):
            ns=pd.Series([x['n'] for x in c]).value_counts().sort_index().to_dict() if c else {}
            print(name,T,'configs',len(c),ns,round(dt),flush=True); out[(name,T)]=c
    pickle.dump(out,open(SC+'compact_cfgs2.pkl','wb')); print('done')
