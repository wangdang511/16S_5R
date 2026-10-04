# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/extend_*
import pickle, json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore'); pd.set_option('display.width',250)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'ext2.py').read().split("rows=[]")[0]
exec(src)
cand,pool=pickle.load(open(SC+'ext_cand.pkl','rb'))
CHOICE={('A2-F',0):3,('A2-R',0):2,('A4-R',0):3,('A5-F',0):3,('A5-F',1):3,('A1-F',0):0,('A1-R',0):1,('A1-R',1):1}
def apply(sites,choice):
    out={}
    for sn,s in sites.items():
        prim=list(s['prim']); ext=list(s['ext'])
        for (cn,i),e in choice.items():
            if cn==sn: prim[i]=cand[f'{sn}|{i}|e{e}']['seq']; ext[i]=e
        out[sn]=dict(s,prim=prim,ext=ext)
    return out
base5=pool                           # 含 A1（V1·V2 均衡版）
ext5=apply(pool,CHOICE)
base4={k:v for k,v in pool.items() if not k.startswith('A1')}; ext4={k:v for k,v in ext5.items() if not k.startswith('A1')}
rows=[]
for nm,sites in (('4 扩增子（原）',base4),('4 扩增子（延长后）',ext4),('5 扩增子（原）',base5),('5 扩增子（延长后）',ext5)):
    ev,sev,hp=evaluate((nm,0),sites); d=build(sites)
    r=dict(design=nm,oligos=ev['oligos'],expansions=ev['expansions'],Tm_min=round(ev['Tm_min'],1),Tm_max=round(ev['Tm_max'],1),severe=ev['severe'],hairpins=ev['hairpins'],SILVA_amp=ev['SILVA_amp_cov'],GG_amp=ev['GG_amp_cov'],SILVA_all=round(ev['SILVA_p_all'],3),GG_all=round(ev['GG_p_all'],3),ideal=round(ev['ideal_acc'],4))
    for lab,a,f in SET: r[lab]=round(realistic2(d,a,f),4)
    rows.append(r); print(r,flush=True); print('  severe',sev,'hp',hp,flush=True)
pd.DataFrame(rows).to_csv(SC+'ext_design_eval.csv',index=False)
oligos=[]
for nm,sites in (('4 扩增子（延长后）',ext4),('V1V2 追加（延长后）',{k:v for k,v in ext5.items() if k.startswith('A1')})):
    for sn,s in sites.items():
        for p in s['prim']: oligos.append(dict(set=nm,site=sn,seq=p))
for nm,sites in (('4 扩增子（延长后）',ext4),('V1V2 追加（延长后）',{k:v for k,v in ext5.items() if k.startswith('A1')})): pass
json.dump({'4扩增子延长后':[dict(name=f"{o['site']}_{i}",seq=o['seq']) for i,o in enumerate(o for o in oligos if o['set'].startswith('4'))],
           'V1V2追加延长后':[dict(name=f"{o['site']}_{i}",seq=o['seq']) for i,o in enumerate(o for o in oligos if o['set'].startswith('V1'))]},open(SC+'ext2/offtarget_primers.json','w'),ensure_ascii=False,indent=1)
json.dump({k:v for k,v in ext5.items()},open(SC+'ext_final_sites.json','w'),ensure_ascii=False,indent=1)
for o in oligos: print(o)
