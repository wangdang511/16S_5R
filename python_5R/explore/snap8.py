# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/snap_products.csv
import pickle, json, itertools, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore'); pd.set_option('display.width',250); pd.set_option('display.max_colwidth',50)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
exec(open(SC+'ext2.py').read().split("rows=[]")[0])
SNAP=pickle.load(open(SC+'snap_sites.pkl','rb')); P=pd.read_csv(SC+'snap_positions.csv')
Fs=['V1_f','V3_f','V4_f',"V6'_f",'V6_f','V7_f']; Rs=['V2_r','V4_r','V5_r','V8_r','V9_r']
pos={g:(int(min(P[P['name'].isin([n for n in P['name'] if n.startswith(g[:-2]) or n==g])].start)),) for g in []}
# 位置（E. coli）：每组引物的 5′ 位置
ppos={'V1_f':(9,27),'V3_f':(341,357),'V4_f':(517,533),"V6'_f":(1055,1070),'V6_f':(967,985),'V7_f':(1099,1114),'V2_r':(372,391),'V4_r':(785,803),'V5_r':(907,926),'V8_r':(1390,1407),'V9_r':(1492,1507)}
# 每条序列各引物组是否结合（≤1 错配、3′ 3 碱基匹配；另算宽松：≤2 错配）
def hits(site_dict, ref, rows_, mm_rule):
    return site_hit(ref,rows_,site_dict)
res={}
for vn,(ref,rows_) in VAL.items():
    phs=ref.phylum[rows_]; H={g:site_hit(ref,rows_,s) for g,s in SNAP.items()}
    res[vn]=(H,phs)
rows=[]
for f,r in itertools.product(Fs,Rs):
    if ppos[r][0]<=ppos[f][1]: continue
    L=ppos[r][1]-ppos[f][0]+1; row=dict(F=f,R=r,length=L)
    for vn,(H,phs) in res.items():
        both=H[f]&H[r]; per=[both[phs==p].mean() for p in pdz.KEY_PHYLA if (phs==p).sum()>=30]
        row[vn+'_all']=round(float(both.mean()),3); row[vn+'_phylum_mean']=round(float(np.mean(per)),3)
    rows.append(row)
R=pd.DataFrame(rows).sort_values('length'); print(R.to_string()); R.to_csv(SC+'snap_products.csv',index=False)
# 引物池二聚体
pool={n:[s] for n,s in zip(P['name'],P['seq'])}
sev,hp=cp.pool_issues(pool)
print('SNAP 引物池严重二聚体',len(sev)); 
for x in sorted(sev,key=lambda t:t[4])[:25]: print(x[0],x[2],x[4],x[5])
print('发夹',hp)
pd.DataFrame([dict(a=x[0],b=x[2],dG=x[4],dG_3end=x[5]) for x in sev]).to_csv(SC+'snap_dimers.csv',index=False)
