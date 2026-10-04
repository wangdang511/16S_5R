# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/extend_*
import pickle, collections, json, sys
import pandas as pd
sys.path.insert(0,'.')
import explore.offtarget as ot, primer_design as pdz
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
sets=json.load(open(SC+'ext2/offtarget_primers.json')); rows=[]; per=[]
for sn,ol in sets.items():
    res,allh=pickle.load(open(SC+f'ext2/offtarget_hits_{sn}.pkl','rb')); H=[h for h in allh[sn] if '_' not in h[0]]
    for mm in (1,2,3):
        P=[p for p in ot.pairs(H,1500,mm) if p[3]>=80]; uniq=set((p[0],p[1]//50,p[2]//50) for p in P)
        loc=set((h[0],h[1]//20,h[2]) for h in H if h[4]<=mm); rows.append(dict(set=sn,mm=mm,oligos=len(ol),loci=len(loc),products=len(uniq)))
    mt=res[(sn,'rCRS')]; print(sn,'rCRS ≤4mm',sorted({(h[1],h[2],h[3].split('#')[0],h[4]) for h in mt if h[4]<=4}))
    byo=collections.defaultdict(set)
    for h in H:
        if h[4]<=2: byo[h[3].split('#')[0]].add((h[0],h[1]//20,h[2]))
    for o in ol: ne=len(pdz.expand(o['seq'])); per.append((sn,o['name'],o['seq'],len(o['seq']),ne,len(byo[o['name']]),round(len(byo[o['name']])/ne,1)))
S=pd.DataFrame(rows); print(S.to_string())
P=pd.DataFrame(per,columns=['set','oligo','seq','nt','expansions','loci_le2mm','loci_per_expansion']); print(P.to_string())
S.to_csv(SC+'ext_offtarget_summary.csv',index=False); P.to_csv(SC+'ext_offtarget_per_oligo.csv',index=False)
