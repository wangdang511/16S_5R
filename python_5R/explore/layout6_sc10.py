# 注：此脚本在本次会话中从 scratchpad 运行（借用 wangdang511/Olivar_primer 的 SADDLE 打分 PrimerSetBadnessFast），依赖未入库的中间文件；结果见 docs/primer_design/layout6_dimer_*.csv
import sys; sys.path.insert(0,'/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad')
from sc2 import *
newR=[('16S-A5-R.1','ATGMTGAYTTGACGTCRTC'),('16S-A5-R.2','ATRCTGACYTGACGTCRTC')]
ids5=[add(n,'L6A5-R','R',1188,s) for n,s in newR]; byslot['L6A5-R']=ids5
W=list(WINS); W[4]=(986,1187); st=STRUCT['L6']; st['wins']=W; st['S']=build_S(W); st['lens']=np.array([b-a+1 for a,b in W],np.float32)
Tp,res=pickle.load(open(SC+'pools_final.pkl','rb'))
name2id={U[i]['name']:i for i in range(len(U)) if U[i]['slot'].startswith('L6')}
name2id['16S-A1-F.SNAP']=name2id['V1_f(SNAP 9–27)']; name2id['16S-A5-R.1']=ids5[0]; name2id['16S-A5-R.2']=ids5[1]
def pool(n): return [r['名称'] for _,r in Tp.iterrows() if r[f'池{n}']=='✓']
ADD=['16S-A1-R.7s','16S-A1-R.9s','16S-A6-F.5s','16S-A1-R.8s','16S-A6-F.4s']
P28=pool(28); P33=P28+ADD


import json
ids=lambda names:{name2id[n] for n in names}
# 最终池：A2-F 换为起点 334 的两条
i_a2=[add('16S-A2-F.1','L6A2-F','F',334,'CCAKACWCCTACGGGAGGC'),add('16S-A2-F.2','L6A2-F','F',334,'CCCSACTCCTACGGGAGRC')]
byslot['L6A2-F']=list(byslot.get('L6A2-F',[]))+i_a2
W=list(STRUCT['L6']['wins']); W[1]=(353,503)
def setW(W):
    S_=STRUCT['L6']; S_['wins']=W; S_['S']=build_S(W); S_['lens']=np.array([b-a+1 for a,b in W],np.float32)

import json, sys
sys.path.insert(0,"/home/user/wangdang511/olivar_primer/src/olivar")
import design as od
Tt=pd.read_pickle(SC+'pool32_table.pkl')
def vols(T,total=2.5):
    v={}
    for s_,g in T.groupby('位点'):
        sup=g[g['名称'].str.endswith('s')]; main=g[~g['名称'].str.endswith('s')]
        for _,r in main.iterrows(): v[r['名称']]=total*(0.9 if len(sup) else 1)*r['展开数']/main['展开数'].sum()
        for _,r in sup.iterrows(): v[r['名称']]=total*0.1*r['展开数']/sup['展开数'].sum()
    return v
def badness(rows):
    # rows: list of (name, site, seq)
    T=pd.DataFrame([dict(名称=n,位点=s_,展开数=len(pdz.expand(q)),seq=q) for n,s_,q in rows]); V=vols(T)
    F=[];R=[];cf=[];cr=[]
    for _,r in T.iterrows():
        ex=pdz.expand(r['seq']); 
        for e in ex: (F if r['位点'].endswith('F') else R).append(e.lower()); (cf if r['位点'].endswith('F') else cr).append(V[r['名称']]/len(ex))
    m=np.mean(cf+cr); tot,comp=od.PrimerSetBadnessFast(F,R,fP_conc=[c/m for c in cf],rP_conc=[c/m for c in cr]); return tot,max(comp[0]+comp[1])
base_rows=[(r['名称'],r['位点'],r['5′→3′序列']) for _,r in Tt.iterrows()]
name_set={n for n,_,_ in base_rows}
i_a2=[name2id.get(n) for n in ('16S-A2-F.1','16S-A2-F.2')]
# 当前 pool32 的 ids：A2-F 为新加两条（name2id 中以 add 的名字？）
cur_ids=set()
for n,_,_ in base_rows:
    if n in ('16S-A2-F.1','16S-A2-F.2'): continue
    cur_ids.add(name2id[n])
a2f=[add('16S-A2-F.1','L6A2-F','F',334,'CCAKACWCCTACGGGAGGC'),add('16S-A2-F.2','L6A2-F','F',334,'CCCSACTCCTACGGGAGRC')]
byslot['L6A2-F']=list(byslot.get('L6A2-F',[]))+a2f
W0=list(STRUCT['L6']['wins']); W0[1]=(353,503)
def setW(W):
    S_=STRUCT['L6']; S_['wins']=W; S_['S']=build_S(W); S_['lens']=np.array([b-a+1 for a,b in W],np.float32)
def run(label,drop=(),a5r=None):
    rows=[r for r in base_rows if r[0] not in drop]
    ids_=set(cur_ids)|set(a2f)
    for n in drop:
        if n in name2id: ids_.discard(name2id[n])
    W=list(W0)
    if a5r:
        for n in ('16S-A5-R.1','16S-A5-R.2'): ids_.discard(name2id[n])
        rows=[r for r in rows if not r[0].startswith('16S-A5-R')]
        for k,q in enumerate(a5r,1):
            i=add(f'16S-A5-R.{k}x','L6A5-R','R',1192,q); byslot['L6A5-R']=list(byslot.get('L6A5-R',[]))+[i]; ids_.add(i); rows.append((f'16S-A5-R.{k}x','A5-R',q))
        W[4]=(986,1191)
    setW(W); r=evaluate(ids_,'L6'); b,mx=badness(rows)
    print(f'{label:34s} 引物{len(ids_):3d} 展开{r["expansions"]:4d} J {r["J"]:.4f} 准确率 {r["accf"]:.4f} 覆盖 {r["cov"]:.4f} 支原体 {r["ten"]:.4f} | Olivar 总坏度 {b:9.0f} 最大单条 {mx:8.0f}',flush=True)
run('当前 32 条池')
run('去掉 A6-F.3',drop=('16S-A6-F.3',))
run('去掉 A2-R.3',drop=('16S-A2-R.3',))
run('去掉 A6-F.3 + A2-R.3',drop=('16S-A6-F.3','16S-A2-R.3'))
run('A5-R 换起点 1192',a5r=['GGSCATGMTGAYTTGACGT','GGCCATRCTGACTWGACRT'])
run('A5-R 换 1192 + 去掉 A6-F.3',a5r=['GGSCATGMTGAYTTGACGT','GGCCATRCTGACTWGACRT'],drop=('16S-A6-F.3',))
run('A5-R 换 1192 + 去 A6-F.3 + A2-R.3',a5r=['GGSCATGMTGAYTTGACGT','GGCCATRCTGACTWGACRT'],drop=('16S-A6-F.3','16S-A2-R.3'))
print('done')
