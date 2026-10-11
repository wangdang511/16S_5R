import pickle, json, sys, time, itertools, numpy as np, pandas as pd, warnings
import primer_design as pdz, iterate_5R as it, compact_primers as cp, tiling_search as ts
warnings.filterwarnings('ignore'); pd.set_option('display.width',250)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
TARGET_TM=60.5
cfgs=pickle.load(open(SC+'compact_cfgs2.pkl','rb'))
SCHEMES={'5 扩增子':[1,2,3,4,5],'4 扩增子':[2,3,4,5],'3 扩增子':[2,3,5]}
def tm_cost(c):
    return sum(max(0,abs(cp.oligo_tm(p)[0]-TARGET_TM)-2)**2 for p in c['prim'])
def span_ok(chosen):
    # chosen: {amp k: (F cfg, R cfg)}
    sp=sorted((f['start'],r['start']+r['L']-1) for f,r in chosen.values())
    if any(b-a+1>290 or b-a+1<150 for a,b in sp): return False
    return all(sp[i+1][0]-sp[i][1]-1>=20 for i in range(len(sp)-1))
def pool(chosen):
    return {f'A{k}-{o}':c['prim'] for k,(f,r) in chosen.items() for o,c in (('F',f),('R',r))}
def cost(chosen):
    if not span_ok(chosen): return 1e9,None,None
    sev,hp=cp.pool_issues(pool(chosen))
    tmc=sum(tm_cost(c) for f,r in chosen.values() for c in (f,r))
    tot=sum(c['fold'] for f,r in chosen.values() for c in (f,r)); n=sum(c['n'] for f,r in chosen.values() for c in (f,r))
    return 1000*len(sev)+10*len(hp)+20*tmc+3*n+0.2*tot,sev,hp
def optimize(amps,T):
    cand={}
    for k in amps:
        for o in 'FR':
            c=cfgs.get((f'A{k}-{o}',T),[])
            if not c: return None,f'A{k}-{o} 在 {T:.0%} 下找不到满足目标的引物集合'
            nmin=min(x['n'] for x in c); cand[(k,o)]=[x for x in c if x['n']==nmin]
    chosen={k:(min(cand[(k,'F')],key=tm_cost),min(cand[(k,'R')],key=tm_cost)) for k in amps}
    cur,sev,hp=cost(chosen)
    for sweep in range(3):
        changed=False
        for k in amps:
            for oi,o in enumerate('FR'):
                best=(cur,chosen[k][oi])
                for c in cand[(k,o)]:
                    trial=dict(chosen); pair=list(chosen[k]); pair[oi]=c; trial[k]=tuple(pair)
                    v,_,_=cost(trial)
                    if v<best[0]-1e-9: best=(v,c)
                if best[1] is not chosen[k][oi]:
                    pair=list(chosen[k]); pair[oi]=best[1]; chosen[k]=tuple(pair); cur=best[0]; changed=True
        print('   sweep',sweep,round(cur,1),flush=True)
        if not changed: break
    cur,sev,hp=cost(chosen)
    return (chosen,cur,sev,hp),None
if __name__=='__main__':
    res={}
    for T in (0.90,0.95):
        for nm,amps in SCHEMES.items():
            t=time.time(); print(nm,T,flush=True)
            r,err=optimize(amps,T)
            if r is None: print('  FAIL',err,flush=True); res[(nm,T)]=dict(error=err); continue
            chosen,c,sev,hp=r
            res[(nm,T)]=dict(chosen=chosen,cost=c,severe=sev,hairpins=hp,amps=amps)
            print('  oligos',sum(x['n'] for f,r_ in chosen.values() for x in (f,r_)),'severe dimers',len(sev),'hairpins',len(hp),round(time.time()-t),flush=True)
    pickle.dump(res,open(SC+'pool_opt2.pkl','wb')); print('done')
