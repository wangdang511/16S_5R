# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/
import pickle, json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore'); pd.set_option('display.width',250)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
exec(open(SC+'ext2.py').read().split("rows=[]")[0])
REC=json.load(open(SC+'final_rec_sites.json')); old=json.load(open(SC+'final_rec_sites_v1.json'))
rows=[]
for nm,sites in (('推荐（改前 A2-F ACT…）',old),('推荐（改后 A2-F ACW…）',REC)):
    ev,sev,hp=evaluate((nm,0),sites); d=build(sites)
    r=dict(design=nm,oligos=ev['oligos'],expansions=ev['expansions'],Tm_min=round(ev['Tm_min'],1),Tm_max=round(ev['Tm_max'],1),severe=ev['severe'],hairpins=ev['hairpins'],SILVA_amp=ev['SILVA_amp_cov'],GG_amp=ev['GG_amp_cov'],SILVA_all=round(ev['SILVA_p_all'],3),GG_all=round(ev['GG_p_all'],3),ideal=round(ev['ideal_acc'],4))
    for lab,a,f in SET: r[lab]=round(realistic2(d,a,f),4)
    rows.append(r); print(r,[ (x[1],x[3],x[4],x[5]) for x in sev],hp)
pd.DataFrame(rows).to_csv(SC+'final_design_eval3.csv',index=False)
