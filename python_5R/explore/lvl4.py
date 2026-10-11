# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/lvl4_*.csv
import sys; sys.path.insert(0,'/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad')
from par2 import *
import compact_primers as cp, os, json
L=pd.read_csv(SC+'par_levels.csv').set_index('level')
SETS={lv:set(int(x) for x in str(L.loc[lv,'sel']).split(',')) for lv in (24,16,12)}
a1f=[i for i,u in enumerate(U) if u['name']=='A1-F'][0]; v1=[i for i,u in enumerate(U) if u['name'].startswith('V1_f')][0]
SETS['12v']=(SETS[12]-{a1f})|{v1}
def acc_credit(P):
    Pm=P.T.astype(np.float32); Lq=(Pm*lens).sum(1); Sq=np.zeros((n,n),np.float32)
    for w in range(K): Sq+=Pm[:,w][:,None]*Sw[w]
    rate=np.where(Lq[:,None]>0,(Lq[:,None]-Sq)/np.maximum(Lq[:,None],1),np.inf); np.fill_diagonal(rate,np.inf)
    mn=rate.min(1,keepdims=True); tied=(rate==mn)&np.isfinite(rate)
    credit=np.where(tied.any(1),(tied&same).sum(1)/np.maximum(tied.sum(1),1),0.0)*(Lq>0)
    return credit
def amp_hits(sel,vn):
    M=[]
    for f,r in PAIRS:
        a=[i for i in byslot.get(f,[]) if i in sel]; b=[i for i in byslot.get(r,[]) if i in sel]
        M.append(np.zeros(len(U[0]['val'][vn]),bool) if not a or not b else (np.any([U[i]['val'][vn] for i in a],0)&np.any([U[i]['val'][vn] for i in b],0)))
    return np.array(M).T
def ten_hits(sel,db):
    M=[]
    for f,r in PAIRS:
        a=[i for i in byslot.get(f,[]) if i in sel]; b=[i for i in byslot.get(r,[]) if i in sel]
        ln=len(U[0]['ten'][db][0])
        M.append(np.zeros(ln,bool) if not a or not b else (np.any([U[i]['ten'][db][0] for i in a],0)&np.any([U[i]['ten'][db][0] for i in b],0)))
    return np.array(M).T
data={}
for lv,sel in SETS.items():
    P=amp_masks(sel); data[lv]=dict(credit=acc_credit(P),GG=amp_hits(sel,'GG'),SILVA=amp_hits(sel,'SILVA'),tGG=ten_hits(sel,'GG'),tSILVA=ten_hits(sel,'SILVA'))
phyl={vn:VAL[vn][0].phylum[VAL[vn][1]] for vn in VAL}
def metrics(d,idx):
    accf=d['credit'][q&idx['acc']].mean() if False else d['credit'][idx['acc']][q[idx['acc']]].mean()
    cov=[]
    for vn in ('GG','SILVA'):
        M=d[vn][idx[vn]]; ph=phyl[vn][idx[vn]]
        per=[M[ph==p].mean(0) for p in pdz.KEY_PHYLA if (ph==p).sum()>=30]; cov.append(np.mean(per))
    ten=[ (d['t'+db][idx['t'+db]].sum(1)>=3).mean() for db in ('GG','SILVA')]
    return accf,float(np.mean(cov)),float(np.mean(ten))
rng=np.random.default_rng(0); B=400
full={'acc':np.arange(n),'GG':np.arange(len(phyl['GG'])),'SILVA':np.arange(len(phyl['SILVA'])),'tGG':np.arange(len(data[24]['tGG'])),'tSILVA':np.arange(len(data[24]['tSILVA']))}
pt={lv:metrics(data[lv],full) for lv in SETS}
boot={lv:[] for lv in SETS}
for b in range(B):
    idx={k:rng.integers(0,len(v),len(v)) for k,v in full.items()}
    for lv in SETS: boot[lv].append(metrics(data[lv],idx))
boot={lv:np.array(v) for lv,v in boot.items()}
Jb={lv:boot[lv].mean(1) for lv in SETS}
res=[]
for lv in SETS:
    a,c,t=pt[lv]; res.append(dict(level=lv,accf=a,cov=c,ten=t,J=(a+c+t)/3,J_lo=np.percentile(Jb[lv],2.5),J_hi=np.percentile(Jb[lv],97.5)))
print(pd.DataFrame(res).round(4).to_string())
for a,b in ((24,16),(16,12),(24,12),('12v',12),(16,'12v'),(24,'12v')):
    d=Jb[a]-Jb[b]; print(f'J({a})-J({b}): {d.mean():.4f}  95% [{np.percentile(d,2.5):.4f}, {np.percentile(d,97.5):.4f}]  P(差>0)={(d>0).mean():.3f}')
    for j,nm in enumerate(('accf','cov','ten')):
        dd=boot[a][:,j]-boot[b][:,j]; print('   ',nm,f'{dd.mean():.4f} [{np.percentile(dd,2.5):.4f},{np.percentile(dd,97.5):.4f}]')
pd.DataFrame(res).to_csv(SC+'lvl4_J.csv',index=False)
pickle.dump((pt,{lv:boot[lv] for lv in SETS}),open(SC+'lvl4_boot.pkl','wb'))
# 每个扩增子的覆盖
names=['V1·V2','V3','V4','V5 (A3-F×V5-R)','V6·V7','V8·V9']
rows=[]
for lv in SETS:
    d=data[lv]
    for vn in ('GG','SILVA'):
        ph=phyl[vn]; per=np.array([d[vn][ph==p].mean(0) for p in pdz.KEY_PHYLA if (ph==p).sum()>=30]); worst=per.min(0); mean=per.mean(0)
        for j,nm in enumerate(names): rows.append(dict(level=lv,db=vn,amplicon=nm,mean=mean[j],worst=worst[j]))
    for db in ('GG','SILVA'):
        M=d['t'+db]
        for j,nm in enumerate(names): rows.append(dict(level=lv,db='Ten_'+db,amplicon=nm,mean=M[:,j].mean(),worst=np.nan))
        rows.append(dict(level=lv,db='Ten_'+db,amplicon='>=3 个扩增子',mean=(M.sum(1)>=3).mean(),worst=np.nan))
    for vn in ('GG','SILVA'): rows.append(dict(level=lv,db=vn,amplicon='全部 6 个扩增子',mean=d[vn].all(1).mean(),worst=np.nan))
pd.DataFrame(rows).to_csv(SC+'lvl4_amp.csv',index=False)
# 池内二聚体 / Tm / 展开
st=[]
for lv,sel in SETS.items():
    pool={}
    for i in sel: pool.setdefault(U[i]['slot'],[]).append(U[i]['seq'])
    sev,hp=cp.pool_issues(pool); tms=[cp.oligo_tm(U[i]['seq'])[0] for i in sel]
    st.append(dict(level=lv,oligos=len(sel),expansions=sum(U[i]['nexp'] for i in sel),severe=len(sev),hairpins=len(hp),Tm_min=round(min(tms),1),Tm_max=round(max(tms),1)))
print(pd.DataFrame(st)); pd.DataFrame(st).to_csv(SC+'lvl4_pool.csv',index=False)
os.makedirs(SC+'lvl',exist_ok=True)
for lv,sel in SETS.items():
    json.dump({f'L{lv}':[dict(name=U[i]['name'],seq=U[i]['seq']) for i in sorted(sel)]},open(SC+f'lvl/offtarget_primers_{lv}.json','w'),ensure_ascii=False)
pickle.dump(SETS,open(SC+'lvl4_sets.pkl','wb'))
