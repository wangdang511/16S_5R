# 注：此脚本在本次会话中从 scratchpad 运行，依赖 S.pkl / big_labeled.pkl 等中间文件（未入库），路径需自行替换；结果见 docs/primer_design/v1v2_*.csv
import pickle, json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore'); pd.set_option('display.width',260)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'robust_realistic.py').read().split("D={'S1'")[0]
exec(src)
df=pd.read_pickle(SC+'v2r_scan.pkl')
def pick(k,fold):
    b=df[(df.k==k)&(df.fold==fold)].sort_values('cov_min',ascending=False).iloc[0]
    return dict(orient='R',start=int(b.start),L=int(b.L),prim=list(b.tops),ext=[0]*k)
CAND={'旧 A1-R (3 条/展开 16)':A1['A1-R'],'2 条/32':pick(2,16),'3 条/22':pick(3,8),'3 条/48':pick(3,16),'4 条/30':pick(4,8),'4 条/64':pick(4,16),'6 条/96':pick(6,16)}
FC={'旧F(19nt,1条)':A1['A1-F'],'F17(1条)':dict(orient='F',start=8,L=17,prim=['AGRGTTTGATYMTGGCT'],ext=[0]),
    'F17(2条)':dict(orient='F',start=8,L=17,prim=['AGRGTTTGATYMTGGCT','YGWGTTTGATCCTGGCK'],ext=[0,0])}
rows=[]; res={}
import itertools
for (fn,F),(nm,r) in itertools.product(FC.items(),CAND.items()):
    sites={'A1-F':F,'A1-R':r}; nm0=nm; nm=fn+' + '+nm
    ev,sev,hp=evaluate((nm,0),sites)
    d=build(sites)
    row=dict(A1R=nm,start=r['start'],L=r['L'],oligos=ev['oligos'],exp=ev['expansions'],Tm_min=ev['Tm_min'],Tm_max=ev['Tm_max'],severe=ev['severe'],
             SILVA_amp=ev.get('SILVA_amp_cov'),GG_amp=ev.get('GG_amp_cov'),held_amp=float(d[0][1].mean()))
    row['F']=fn; rows.append(row); res[nm]=sites; print(row,flush=True)
pd.DataFrame(rows).to_csv(SC+'v2_a1r.csv',index=False)
# 整体设计
OUT=[]
c4=store[('4 扩增子',0.95)]['sites']
SET=[('abs60',60,0.0),('abs150',150,0.0),('abs300',300,0.0),('frac0.5',0,0.5),('frac0.8',0,0.8)]
for nm in ('旧F(19nt,1条) + 旧 A1-R (3 条/展开 16)','F17(1条) + 2 条/32','F17(1条) + 3 条/22','F17(1条) + 4 条/64','F17(2条) + 4 条/64'):
    A=res[nm]
    for base,bs in (('S1',s1),('紧凑 4 扩增子 95%',c4)):
        sites={**A,**{f'A{int(k[1])+1}-{k[3]}':v for k,v in bs.items()}}
        d=build(sites); r=dict(design=f'{base} + V1V2[{nm}]')
        for lab,a,f in SET: r[lab]=realistic2(d,a,f)
        allh=[((a,b),np.ones(len(held),bool)) for (a,b),_ in d]; r['ideal_on_set']=realistic2(allh,60,0.0)
        ev,_,_=evaluate((nm,0),{k:v for k,v in sites.items()}) if False else (None,None,None)
        OUT.append(r); print({k:(round(v,3) if isinstance(v,float) else v) for k,v in r.items()},flush=True)
pd.DataFrame(OUT).to_csv(SC+'v2_designs.csv',index=False)
