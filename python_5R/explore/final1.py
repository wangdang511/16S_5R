# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/order_list_full.csv 与 final_*
import pickle, json, numpy as np, pandas as pd, warnings, collections
warnings.filterwarnings('ignore'); pd.set_option('display.width',250); pd.set_option('display.max_colwidth',60)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'ten3.py').read().split("rows=[addon(x)")[0]
exec(src)       # 含 ten1 的 refs/site_cov_ten、D1、addon 函数
orig=json.load(open(SC+'a5f_final5.json')); V3L=json.load(open(SC+'v3_sites_k1.json'))['B']
base=json.loads(json.dumps(D1))
base['A4-F']=orig['A4-F']                      # 保留 906–927
sites_short=json.loads(json.dumps(base))       # V3 短版 505–523
sites_long=json.loads(json.dumps(base)); sites_long['A2-R']=V3L['A2-R']   # V3 长版 516–535
# 补充寡核苷酸（支原体）：对每个弱位点，重新用当前位点定义设计
def addon_for(sites,site,fold=8,deg=3):
    global D1j; D1j=sites; return addon(site,fold,deg)
sup={}
for nm,sites in (('short',sites_short),('long',sites_long)):
    for site in ('A1-R','A2-F','A2-R','A4-R','A5-F'):
        if site in sup.get(nm,{}): continue
        r=addon_for(sites,site); sup.setdefault(nm,{})[site]=r
        print(nm,site,r['addon'],r['nexp'],r['Tm'],r['GG_before'],'->',r['GG_after'],r['SILVA_before'],'->',r['SILVA_after'])
pickle.dump((sites_short,sites_long,sup),open(SC+'final1.pkl','wb'))
json.dump({'short':sites_short,'long':sites_long},open(SC+'final1_sites.json','w'),ensure_ascii=False,indent=1)
