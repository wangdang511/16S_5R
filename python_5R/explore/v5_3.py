# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/v5_* 与 order_list_full_v5.csv
import pickle, json, numpy as np, pandas as pd, warnings, itertools
warnings.filterwarnings('ignore'); pd.set_option('display.width',250); pd.set_option('display.max_colwidth',60)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'ten1.py').read().split("rows=[]; masks={}")[0]
exec(src)
exec(open(SC+'ext2.py').read().split("rows=[]")[0].split("sites=")[0]) if False else None
REC=json.load(open(SC+'final_rec_sites.json')); D967=json.load(open(SC+'short_sites.json'))['D1']
mk=lambda o,start,L,prims: dict(orient=o,start=start,L=L,prim=list(prims),ext=[0]*len(prims))
V5R=mk('R',907,21,['CCCGTCAATTYCTTTRAGTTT','CCCGTCAATTCCTTTGAGYYY'])
A4F967=json.loads(json.dumps(D967['A4-F'])); A4F967['prim'].append('RTACMCGAARAACCTTACC'); A4F967['ext'].append(0)   # 含支原体补充（未延长）
SB=json.loads(json.dumps(REC)); SB['A4-F']=A4F967
# 池内二聚体
pool={}
for n,s in SB.items(): pool[n]=s['prim']
pool['V5-R']=V5R['prim']
sev,hp=cp.pool_issues(pool); print('S_B 池（含 V5R）寡核苷酸',sum(len(v) for v in pool.values()),'严重二聚体',[(x[0],x[1],x[2],x[3],x[4],x[5]) for x in sev],'发夹',hp)
# 引物对覆盖（GG/SILVA 主要门平均/最差、支原体）
def pair_cov(f,r):
    out={}
    for vn,(ref,rows_) in VAL.items():
        phs=ref.phylum[rows_]; h=site_hit(ref,rows_,f)&site_hit(ref,rows_,r); per=[h[phs==p].mean() for p in pdz.KEY_PHYLA if (phs==p).sum()>=30]; out[vn]=(round(float(np.mean(per)),3),round(float(min(per)),3))
    ten={}
    for db,(ref,rr) in refs.items():
        hf,vf=site_cov_ten(f,ref,rr); hr,vr=site_cov_ten(r,ref,rr); v=vf&vr; ten[db]=round(float((hf&hr)[v].mean()),3)
    return out,ten
rows=[]
for nm,(f,r) in {'V5（A3-F × V5R，371 bp）':(SB['A3-F'],V5R),'V4（A3-F × A3-R，249 bp）':(SB['A3-F'],SB['A3-R']),'V6·V7 967 版（229 bp）':(SB['A4-F'],SB['A4-R']),'V6·V7 906 版（原，290 bp）':(REC['A4-F'],REC['A4-R'])}.items():
    o,t=pair_cov(f,r); rows.append(dict(pair=nm,GG=o['GG'],SILVA=o['SILVA'],Ten_GG=t['GG'],Ten_SILVA=t['SILVA']))
print(pd.DataFrame(rows).to_string()); pd.DataFrame(rows).to_csv(SC+'v5_pair_cov.csv',index=False)
# 属准确率（realistic）：加 V5 窗口，A4 用 967 版
def merge(iv):
    iv=sorted((a,b) for a,b in iv if b>=a); out=[]
    for a,b in iv:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
exec(open(SC+'ext2.py').read().split("rows=[]")[0])
def design_of(sites,extra=None):
    d=build(sites)
    if extra: d+=extra
    return d
hitV5=site_hit(Gh,np.arange(len(held)),SB['A3-F'])&site_hit(Gh,np.arange(len(held)),V5R)
res=[]
for nm,(sites,extra) in {'当前推荐（906 正向，无 V5）':(REC,None),'仅缩短 V6·V7（967 版，无 V5）':({**REC,'A4-F':A4F967},None),'S_B：967 版 + V5（V5 由 A3-F×V5R 提供）':({**REC,'A4-F':A4F967},[((788,906),hitV5)])}.items():
    d=design_of(sites,extra); r=dict(design=nm)
    for lab,a,f in SET: r[lab]=round(realistic2(d,a,f),4)
    ints=[iv for iv,_ in d]; r['ideal']=round(ts.accuracy_direct(xall,merge(ints),gall,vall),4); res.append(r)
print(pd.DataFrame(res).to_string()); pd.DataFrame(res).to_csv(SC+'v5_design_eval.csv',index=False)
pickle.dump((SB,V5R),open(SC+'v5_sb.pkl','wb'))
