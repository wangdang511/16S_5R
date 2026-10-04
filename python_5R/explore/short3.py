# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/short_*
import pickle, json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore'); pd.set_option('display.width',250); pd.set_option('display.max_colwidth',60)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'ext2.py').read().split("rows=[]")[0]
exec(src)
old=json.load(open(SC+'a5f_final5.json'))
def mk(o,start,L,prims): return dict(orient=o,start=start,L=L,prim=prims,ext=[0]*len(prims))
A4F=mk('F',967,19,['MWACGCGARGAACCTTACC','MWACMCGAAGAACCTTACC']); A4R=mk('R',1177,19,['ACGTCRTCCYCACCTTCCY'])
A5R=mk('R',1389,19,['GACGGGCGGTGTGTACAAD','GACGGGCGGTGTGTRCAAS'])
D1=json.loads(json.dumps(old)); D1['A4-F']=A4F; D1['A4-R']=A4R
D2=json.loads(json.dumps(D1)); D2['A5-R']=A5R
D1b=json.loads(json.dumps(old)); D1b['A4-F']=A4F   # 只换 F，R 保持原来的 A4-R
rows=[]
for nm,sites in (('原（A4 290 bp，A5 288 bp）',old),('A4 缩短（229 bp）',D1),('A4 缩短，保留原 A4-R',D1b),('A4 + A5 都缩短（229 / 185 bp）',D2)):
    for pn,keep in (('5',lambda k:True),('4',lambda k:not k.startswith('A1'))):
        S={k:v for k,v in sites.items() if keep(k)}
        ev,sev,hp=evaluate((nm,0),S); d=build(S)
        r=dict(design=f'{pn} 扩增子 · {nm}',oligos=ev['oligos'],expansions=ev['expansions'],Tm_min=round(ev['Tm_min'],1),Tm_max=round(ev['Tm_max'],1),severe=ev['severe'],hairpins=ev['hairpins'],spans=ev['spans'],SILVA_amp=ev['SILVA_amp_cov'],GG_amp=ev['GG_amp_cov'],SILVA_all=round(ev['SILVA_p_all'],3),GG_all=round(ev['GG_p_all'],3),ideal=round(ev['ideal_acc'],4))
        for lab,a,f in SET: r[lab]=round(realistic2(d,a,f),4)
        rows.append(r); print({k:v for k,v in r.items() if k!='spans'},'\n   ',r['spans'],' sev',[(x[1],x[3],x[4]) for x in sev],flush=True)
pd.DataFrame(rows).to_csv(SC+'short_design_eval.csv',index=False)
json.dump({'D1':D1,'D1b':D1b,'D2':D2},open(SC+'short_sites.json','w'),ensure_ascii=False,indent=1)
