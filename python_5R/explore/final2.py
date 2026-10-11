# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/order_list_full.csv 与 final_*
import pickle, json, numpy as np, pandas as pd, warnings, collections, os
warnings.filterwarnings('ignore'); pd.set_option('display.width',250); pd.set_option('display.max_colwidth',70)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'ext2.py').read().split("rows=[]")[0]
exec(src)
sites_short,sites_long,sup=pickle.load(open(SC+'final1.pkl','rb'))
S967=json.load(open(SC+'short_sites.json'))['D1']['A4-F']
SUP={'A1-R':sup['short']['A1-R']['addon'],'A2-F':sup['short']['A2-F']['addon'],'A4-R':sup['short']['A4-R']['addon'],'A5-F':sup['short']['A5-F']['addon']}
SUP_ALT={'A2-R(短版)':sup['short']['A2-R']['addon'],'A4-F(967)':'RTACMCGAARAACCTTACC'}
core=json.loads(json.dumps(sites_long))        # 推荐：V3 长版、V6V7 906
def with_sup(sites,names):
    S=json.loads(json.dumps(sites))
    for n in names: S[n]['prim']=S[n]['prim']+[SUP[n]]; S[n]['ext']=S[n]['ext']+[0]
    return S
cfgs={'核心（不含支原体补充）':core,'核心 + 4 条支原体补充':with_sup(core,['A1-R','A2-F','A4-R','A5-F']),'核心 + A4-R、A5-F 补充':with_sup(core,['A4-R','A5-F']),'V3 短版 + 4 条补充（含短版补充）':None}
S_short=json.loads(json.dumps(sites_short)); 
S_short=with_sup(S_short,['A1-R','A2-F','A4-R','A5-F']); S_short['A2-R']['prim']=S_short['A2-R']['prim']+[SUP_ALT['A2-R(短版)']]; S_short['A2-R']['ext']=S_short['A2-R']['ext']+[0]
cfgs['V3 短版 + 补充']=S_short; del cfgs['V3 短版 + 4 条补充（含短版补充）']
rows=[]
for nm,sites in cfgs.items():
    ev,sev,hp=evaluate((nm,0),sites); d=build(sites)
    r=dict(design=nm,oligos=ev['oligos'],expansions=ev['expansions'],Tm_min=round(ev['Tm_min'],1),Tm_max=round(ev['Tm_max'],1),severe=ev['severe'],hairpins=ev['hairpins'],spans=ev['spans'],SILVA_amp=ev['SILVA_amp_cov'],GG_amp=ev['GG_amp_cov'],SILVA_all=round(ev['SILVA_p_all'],3),GG_all=round(ev['GG_p_all'],3),ideal=round(ev['ideal_acc'],4))
    for lab,a,f in SET: r[lab]=round(realistic2(d,a,f),4)
    rows.append(r); print({k:v for k,v in r.items() if k!='spans'},'\n   ',r['spans'],[ (x[1],x[3],x[4],x[5]) for x in sev],hp,flush=True)
pd.DataFrame(rows).to_csv(SC+'final_design_eval.csv',index=False)
# 全体寡核苷酸（含备选）两两检查
allo={}
for s,v in sites_long.items():
    for i,p in enumerate(v['prim']): allo[f'{s}.{i+1}']=p
allo['A2-R.alt短版']=sites_short['A2-R']['prim'][0]
for i,p in enumerate(S967['prim']): allo[f'A4-F.alt967.{i+1}']=p
for k,p in SUP.items(): allo[f'{k}.补充']=p
for k,p in SUP_ALT.items(): allo[f'{k}.补充']=p
sev,hp=cp.pool_issues({k:[p] for k,p in allo.items()}); print('全部寡核苷酸（含备选）严重二聚体:',[(x[0],x[2],x[4],x[5]) for x in sev],'发夹:',hp)
json.dump(allo,open(SC+'final_all_oligos.json','w'),ensure_ascii=False,indent=1)
os.makedirs(SC+'fin',exist_ok=True)
json.dump({'ALL':[dict(name=k,seq=v) for k,v in allo.items()]},open(SC+'fin/offtarget_primers.json','w'),ensure_ascii=False)
pickle.dump(cfgs,open(SC+'final_cfgs.pkl','wb')); print(len(allo))
