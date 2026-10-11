# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/lvl4_*.csv
import sys; sys.path.insert(0,'/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad')
from par2 import *
import compact_primers as cp
SETS=pickle.load(open(SC+'lvl4_sets.pkl','rb'))
O5=pd.read_csv(SC+'order_list_full_v5.csv'); hum={r.seq:r.human_per_exp for r in O5.itertuples()}
hum['GAGTTTGATCMTGGCTCAG']=25.0; hum['CCGTCAATTCMTTTGAGTTT']=11.0
def pos(u):
    s=u['s']; e=s['ext'][0]; return (s['start']-e,s['start']+s['L']-1) if s['orient']=='F' else (s['start'],s['start']+s['L']-1+e)
rows=[]
for lv,sel in SETS.items():
    for slot in ['A1-F','A1-R','A2-F','A2-R','A3-F','A3-R','V5-R','A4-F','A4-R','A5-F','A5-R']:
        idx=[i for i in byslot[slot] if i in sel]
        if not idx: continue
        sup=[i for i in idx if '(补)' in U[i]['name']]; core=[i for i in idx if i not in sup]
        cov={i:U[i]['held'].mean() for i in core}; tot=sum(cov.values()); sf=0.10*len(sup) if core else 1.0/len(sup)
        for i in idx:
            share=(0.10 if i in sup else cov[i]/tot*(1-sf)) if core else 1.0/len(idx)
            a,b=pos(U[i]); p=U[i]['seq']; ts=[cp.oligo_tm(x)[0] for x in pdz.expand(p)]
            rows.append(dict(level=lv,i=i,name=U[i]['name'],slot=slot,orient=slot[-1],seq=p,nt=len(p),pos=f'{a}–{b}',Tm=round(cp.oligo_tm(p)[0],1),Tm_lo=round(min(ts),1),Tm_hi=round(max(ts),1),GC=round(sum(c in 'GCS' for c in p)/len(p),2),expansions=U[i]['nexp'],human_per_exp=hum.get(p,np.nan),is_sup=i in sup,share=round(share,3),uL=round(share*2.5,2)))
D=pd.DataFrame(rows); D.to_csv(SC+'order_levels_long4.csv',index=False)
print(D[D.human_per_exp.isna()][['name','seq']])
for lv in SETS: print(lv,D[D.level==lv].uL.sum(), len(D[D.level==lv]))
print(D[D.level==16][['name','seq','pos','Tm','expansions','human_per_exp','share','uL']].to_string())
