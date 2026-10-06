# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/par_*
import pickle, json, itertools, time, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore'); pd.set_option('display.width',250); pd.set_option('display.max_colwidth',50)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'ten1.py').read().split("rows=[]; masks={}")[0]
exec(src)                       # VAL, refs, site_cov_ten, site_hit 等
SB,V5R=pickle.load(open(SC+'v5_sb.pkl','rb')); SNAP=pickle.load(open(SC+'snap_sites.pkl','rb'))
REC=json.load(open(SC+'final_rec_sites.json'))
# ---- 寡核苷酸池
U=[]   # dict(name, slot, site dict(single))
def add(slot,name,s,i):
    U.append(dict(name=name,slot=slot,s=dict(orient=s['orient'],start=s['start'],L=s['L'],prim=[s['prim'][i]],ext=[s['ext'][i]]),seq=s['prim'][i],nexp=len(pdz.expand(s['prim'][i]))))
labels={'A1-F':['A1-F'],'A1-R':['A1-R.1','A1-R.2','A1-R.3','A1-R.4(补)'],'A2-F':['A2-F.1','A2-F.2(补)'],'A2-R':['A2-R'],'A3-F':['A3-F.1','A3-F.2'],'A3-R':['A3-R'],'A4-F':['A4-F.1','A4-F.2','A4-F.3(补)'],'A4-R':['A4-R.1','A4-R.2(补)'],'A5-F':['A5-F.1','A5-F.2','A5-F.3(补)'],'A5-R':['A5-R']}
for slot,names in labels.items():
    for i,n in enumerate(names): add(slot,n,SB[slot],i)
for i,n in enumerate(['V5-R.1','V5-R.2']): add('V5-R',n,V5R,i)
# 借鉴成熟引物（SNAP）
v1=SNAP['V1_f']; U.append(dict(name='V1_f（SNAP）',slot='A1-F',s=dict(v1,prim=[v1['prim'][0]],ext=[0]),seq=v1['prim'][0],nexp=len(pdz.expand(v1['prim'][0]))))
v5=SNAP['V5_r']; U.append(dict(name='V5_r（SNAP 907R）',slot='V5-R',s=dict(v5,prim=[v5['prim'][0]],ext=[0]),seq=v5['prim'][0],nexp=len(pdz.expand(v5['prim'][0]))))
print(len(U),'oligos'); 
for u in U: print(u['slot'],u['name'],u['seq'],u['nexp'])
# ---- 命中矩阵
rowsH=np.arange(len(held))
t=time.time()
for u in U:
    u['held']=site_hit(Gh,rowsH,u['s'])
    u['val']={vn:site_hit(ref,rows_,u['s']) for vn,(ref,rows_) in VAL.items()}
    u['ten']={db:site_cov_ten(u['s'],ref,rr) for db,(ref,rr) in refs.items()}
print('命中矩阵',round(time.time()-t,1),'s')
pickle.dump(U,open(SC+'par_U.pkl','wb'))
