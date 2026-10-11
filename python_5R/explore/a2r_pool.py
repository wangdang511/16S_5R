# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/a2r_*
import pickle, json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore'); pd.set_option('display.width',250)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'ext2.py').read().split("rows=[]")[0]
exec(src)
fin=json.load(open(SC+'ext_final_sites.json')); cands=pickle.load(open(SC+'a2r_cands.pkl','rb'))
new=cands['A2R_505_L19_k1_f8']; new2=cands['A2R_509_L18_k2_f4']
rows=[]
variants={'A2-R 原（512，16 nt）':None,'A2-R 延长+2（69.6 °C）':'ext','A2-R 均衡（505，19 nt）':new,'A2-R 均衡 2 条（509，18 nt）':new2}
sites_out={}
for vn,v in variants.items():
    for pn,keep in (('4',lambda k:not k.startswith('A1')),('5',lambda k:True)):
        sites={k:s for k,s in fin.items() if keep(k)}
        if v=='ext': pass
        elif v is None: sites['A2-R']=pool['A2-R']
        else: sites['A2-R']=v
        nm=f'{pn} 扩增子 · {vn}'
        ev,sev,hp=evaluate((nm,0),sites); d=build(sites)
        r=dict(design=nm,oligos=ev['oligos'],expansions=ev['expansions'],Tm_min=round(ev['Tm_min'],1),Tm_max=round(ev['Tm_max'],1),severe=ev['severe'],SILVA_amp=ev['SILVA_amp_cov'],GG_amp=ev['GG_amp_cov'],SILVA_all=round(ev['SILVA_p_all'],3),GG_all=round(ev['GG_p_all'],3),ideal=round(ev['ideal_acc'],4))
        for lab,a,f in SET: r[lab]=round(realistic2(d,a,f),4)
        rows.append(r); print(r,flush=True)
        if v is not None and v!='ext' or v=='ext': sites_out[nm]=sites
pd.DataFrame(rows).to_csv(SC+'a2r_design_eval.csv',index=False)
json.dump(sites_out['4 扩增子 · A2-R 均衡（505，19 nt）'],open(SC+'a2r_final4.json','w'),ensure_ascii=False,indent=1)
json.dump(sites_out['5 扩增子 · A2-R 均衡（505，19 nt）'],open(SC+'a2r_final5.json','w'),ensure_ascii=False,indent=1)
mk=lambda sites:[dict(name=f'{sn}_{i}',seq=p) for sn,s in sites.items() for i,p in enumerate(s['prim'])]
json.dump({'4均衡':mk(sites_out['4 扩增子 · A2-R 均衡（505，19 nt）']),'V1V2均衡':mk({k:v for k,v in sites_out['5 扩增子 · A2-R 均衡（505，19 nt）'].items() if k.startswith('A1')})},open(SC+'a2r/offtarget_primers.json','w'),ensure_ascii=False,indent=1)
