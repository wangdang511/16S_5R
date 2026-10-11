# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换
import sys; sys.path.insert(0,'/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad')
from sc2 import *
newR=[('16S-A5-R.1','ATGMTGAYTTGACGTCRTC'),('16S-A5-R.2','ATRCTGACYTGACGTCRTC')]
ids5=[add(n,'L6A5-R','R',1188,s) for n,s in newR]; byslot['L6A5-R']=ids5
W=list(WINS); W[4]=(986,1187); st=STRUCT['L6']; st['wins']=W; st['S']=build_S(W); st['lens']=np.array([b-a+1 for a,b in W],np.float32)
Tp,res=pickle.load(open(SC+'pools_final.pkl','rb'))
name2id={U[i]['name']:i for i in range(len(U)) if U[i]['slot'].startswith('L6')}
name2id['16S-A1-F.SNAP']=name2id['V1_f(SNAP 9–27)']; name2id['16S-A5-R.1']=ids5[0]; name2id['16S-A5-R.2']=ids5[1]
def pool(n): return [r['名称'] for _,r in Tp.iterrows() if r[f'池{n}']=='✓']
ADD=['16S-A1-R.7s','16S-A1-R.9s','16S-A6-F.5s']; ADD2=ADD+['16S-A1-R.8s','16S-A6-F.4s']
P28=pool(28); P31=P28+ADD; P33=P28+ADD2; P36=pool(36)
pdz_ids=lambda names:{name2id[n] for n in names}
out={}
for nm,names in (('28',P28),('28+3',P31),('33',P33),('36',P36)):
    r=evaluate(pdz_ids(names),'L6'); out[nm]=r
    print(nm,len(names),'展开',r['expansions'],{k:round(v,4) for k,v in r.items() if isinstance(v,float)},flush=True)
# 属级匹配率：每个位点 hit(GG)；扩增子 = F 且 R
AMP=[('A1-F','A1-R'),('A2-F','A2-R'),('A3-F','A3-R'),('A4-F','A4-R'),('A5-F','A5-R'),('A6-F','A6-R')]
def sitehit(names,site):
    h=np.zeros(len(G.aln),bool); dat=np.zeros(len(G.aln),bool)
    for n in names:
        if n.replace('16S-','').split('.')[0]!=site: continue
        i=name2id[n]; s=U[i]['s']; M=pdz.site_matrix(G,s['start'],s['L']); ok=(M!=0).all(1)&(M!=45).all(1)
        top=U[i]['seq'] if s['orient']=='F' else pdz.revcomp(U[i]['seq'])
        hh=np.zeros(len(M),bool); hh[ok]=pdz._match(M[ok],top,'right' if s['orient']=='F' else 'left'); h|=hh; dat|=ok
    return h,dat
gen=gen0
GEN=['Prevotella','Bifidobacterium','Neisseria','Roseburia','Dialister','Blautia','Megasphaera','Veillonella','Oscillospira','Ruminococcus','Faecalibacterium','Streptococcus','Staphylococcus','Haemophilus','Mycoplasma','Ureaplasma','Lactobacillus','Bacteroides','Enterococcus','Gemella']
rows=[]
for g in GEN:
    m=gen==g
    r={'属':g,'n':int(m.sum())}
    for nm,names in (('28',P28),('28+3',P31),('33',P33)):
        hs={s:sitehit(names,s) for s in {x for p in AMP for x in p}}
        vals=[]
        for f,rr in AMP:
            ok=m&hs[f][1]&hs[rr][1]
            vals.append(round(100*(hs[f][0]&hs[rr][0])[ok].mean()) if ok.sum()>=3 else None)
        r[nm]=vals
    rows.append(r)
for r in rows:
    if r['28']!=r['33']: print('%-16s n=%3d  28:%s | 28+3:%s | 33:%s'%(r['属'],r['n'],r['28'],r['28+3'],r['33']))
pickle.dump((out,rows),open(SC+'sc6_out.pkl','wb')); print('done'); pickle.dump(dict(P33=P33,P31=P31,info={i:(U[i]['name'],U[i]['seq'],U[i]['s']['start'],U[i]['s']['orient'],U[i]['nexp']) for i in range(len(U)) if U[i]['slot'].startswith('L6')}),open(SC+'sc6_pools.pkl','wb'))
