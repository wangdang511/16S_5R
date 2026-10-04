# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/a5f_*
import pickle, json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore'); pd.set_option('display.width',250)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'ext2.py').read().split("rows=[]")[0]
exec(src)
old=json.load(open(SC+'a2r_final5.json')); new=json.loads(json.dumps(old)); new['A5-F']['prim'][0]='GCKACACRCRTGCTACAAT'
rows=[]
for nm,sites in (('原（A5-F.1 = GCTRCACRCRTGCTACAAT）',old),('新（A5-F.1 = GCKACACRCRTGCTACAAT）',new)):
    for pn,keep in (('5',lambda k:True),('4',lambda k:not k.startswith('A1'))):
        S={k:v for k,v in sites.items() if keep(k)}
        ev,sev,hp=evaluate((nm,0),S); d=build(S)
        r=dict(design=f'{pn} 扩增子 · {nm}',oligos=ev['oligos'],expansions=ev['expansions'],Tm_min=round(ev['Tm_min'],1),Tm_max=round(ev['Tm_max'],1),severe=ev['severe'],hairpins=ev['hairpins'],SILVA_amp=ev['SILVA_amp_cov'],GG_amp=ev['GG_amp_cov'],SILVA_all=round(ev['SILVA_p_all'],3),GG_all=round(ev['GG_p_all'],3),ideal=round(ev['ideal_acc'],4))
        for lab,a,f in SET: r[lab]=round(realistic2(d,a,f),4)
        rows.append(r); print(r,flush=True); print('  sev',sev,flush=True)
pd.DataFrame(rows).to_csv(SC+'a5f_design_eval.csv',index=False)
json.dump(new,open(SC+'a5f_final5.json','w'),ensure_ascii=False,indent=1)
mk=lambda sites:[dict(name=f'{sn}_{i}',seq=p) for sn,s in sites.items() for i,p in enumerate(s['prim'])]
import os; os.makedirs(SC+'a5f',exist_ok=True)
json.dump({'5最终':mk(new)},open(SC+'a5f/offtarget_primers.json','w'),ensure_ascii=False,indent=1)
