# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/v3_*
import pickle, json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore'); pd.set_option('display.width',250); pd.set_option('display.max_colwidth',60)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'ext2.py').read().split("rows=[]")[0]
exec(src)
D1=json.load(open(SC+'short_sites.json'))['D1']; O=pd.read_pickle(SC+'v3_pairs.pkl'); S=pd.read_pickle(SC+'v3_scan.pkl')
def mk(o,start,L,prims): return dict(orient=o,start=start,L=L,prim=list(prims),ext=[0]*len(prims))
Rrow=S[(S.kind=='A4R')&(S.pos==516)&(S.L==20)&(S.severe==0)].sort_values('cov_min',ascending=False)
print(Rrow[['k','fold','cov_min','tm','nexp','prim']].head(6).to_string())
rr=Rrow[Rrow.k==1].iloc[0]; Rn=mk('R',516,20,rr.prim)
# 变体 b：只延长 R（F 不变）；变体 c：F 换到 345、R 516
f345=S[(S.kind=='A4F')&(S.pos==345)&(S.L==19)&(S.severe==0)].sort_values('cov_min',ascending=False).iloc[0]; Fn=mk('F',345-19+1,19,f345.prim)
B=json.loads(json.dumps(D1)); Cc=json.loads(json.dumps(D1)); Cc['A2-R']=Rn; Cc2=json.loads(json.dumps(D1)); Cc2['A2-R']=Rn; Cc2['A2-F']=Fn
rows=[]
for nm,sites in (('原（V3 338–523，186 bp）',B),('V3 反向延长（338–535，198 bp）',Cc),('V3 两端延长（327–535，209 bp）',Cc2)):
    ev,sev,hp=evaluate((nm,0),sites); d=build(sites)
    r=dict(design=nm,oligos=ev['oligos'],expansions=ev['expansions'],Tm_min=round(ev['Tm_min'],1),Tm_max=round(ev['Tm_max'],1),severe=ev['severe'],hairpins=ev['hairpins'],spans=ev['spans'],SILVA_amp=ev['SILVA_amp_cov'],GG_amp=ev['GG_amp_cov'],SILVA_all=round(ev['SILVA_p_all'],3),GG_all=round(ev['GG_p_all'],3),ideal=round(ev['ideal_acc'],4))
    for lab,a,f in SET: r[lab]=round(realistic2(d,a,f),4)
    rows.append(r); print({k:v for k,v in r.items() if k!='spans'},'\n   ',r['spans'],[ (x[1],x[3],x[4]) for x in sev],flush=True)
pd.DataFrame(rows).to_csv(SC+'v3_design_eval_k1.csv',index=False)
json.dump({'B':Cc,'C':Cc2},open(SC+'v3_sites_k1.json','w'),ensure_ascii=False,indent=1)
