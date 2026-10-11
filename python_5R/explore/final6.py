# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/order_list_full.csv 与 final_*
import pickle, json, numpy as np, pandas as pd, warnings, os
warnings.filterwarnings('ignore'); pd.set_option('display.width',250); pd.set_option('display.max_colwidth',60)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
exec(open(SC+'ext2.py').read().split("rows=[]")[0])
rec,SUPSEQ,CH=pickle.load(open(SC+'final5.pkl','rb'))
sites_short,sites_long,sup=pickle.load(open(SC+'final1.pkl','rb'))
def add_sup(sites,names):
    S=json.loads(json.dumps(sites))
    for n in names:
        q=SUPSEQ[n]; S[n]['prim']=S[n]['prim']+[q]; S[n]['ext']=S[n]['ext']+[len(q)-S[n]['L']]
    return S
cfg={'核心（不含补充）':sites_long,'推荐：核心 + 4 条支原体补充（延长后）':add_sup(sites_long,['A1-R','A2-F','A4-R','A5-F']),'核心 + 3 条补充（不含 A1-R）':add_sup(sites_long,['A2-F','A4-R','A5-F'])}
rows=[]
for nm,sites in cfg.items():
    ev,sev,hp=evaluate((nm,0),sites); d=build(sites)
    r=dict(design=nm,oligos=ev['oligos'],expansions=ev['expansions'],Tm_min=round(ev['Tm_min'],1),Tm_max=round(ev['Tm_max'],1),severe=ev['severe'],hairpins=ev['hairpins'],SILVA_amp=ev['SILVA_amp_cov'],GG_amp=ev['GG_amp_cov'],SILVA_all=round(ev['SILVA_p_all'],3),GG_all=round(ev['GG_p_all'],3),ideal=round(ev['ideal_acc'],4))
    for lab,a,f in SET: r[lab]=round(realistic2(d,a,f),4)
    rows.append(r); print({k:v for k,v in r.items()},[ (x[1],x[3],x[4],x[5]) for x in sev],hp,flush=True)
pd.DataFrame(rows).to_csv(SC+'final_design_eval2.csv',index=False)
REC=cfg['推荐：核心 + 4 条支原体补充（延长后）']; json.dump(REC,open(SC+'final_rec_sites.json','w'),ensure_ascii=False,indent=1); pickle.dump(cfg,open(SC+'final_cfgs2.pkl','wb'))
os.makedirs(SC+'fin3',exist_ok=True)
json.dump({'REC':[dict(name=f'{sn}_{i}',seq=p) for sn,s in REC.items() for i,p in enumerate(s['prim'])]},open(SC+'fin3/offtarget_primers.json','w'),ensure_ascii=False)
