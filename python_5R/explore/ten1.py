# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/ten_*
import pickle, json, numpy as np, pandas as pd, warnings, collections
warnings.filterwarnings('ignore'); pd.set_option('display.width',250); pd.set_option('display.max_colwidth',50)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'robust_realistic.py').read().split("def realistic2")[0]
exec(src)
rf=open(SC+'refine.py').read(); exec(rf[rf.index('def extend'):rf.index('def site_cov')]) if False else None
S_=pickle.load(open(SC+'S.pkl','rb')); B=pickle.load(open(SC+'big_labeled.pkl','rb')); tax=B['tax']; ids=B['ids']
ph0=np.array([tax[i][1] if tax[i][1] else 'unclassified' for i in ids]); gen0=np.array([tax[i][5] if len(tax[i])>5 and tax[i][5] else '' for i in ids]); fam0=np.array([tax[i][4] if len(tax[i])>4 and tax[i][4] else '' for i in ids])
G=pdz.Reference(B['aln'],ids,np.array(['Bacteria']*len(ids)),ph0,S_.col_of_pos,S_.ref_seq,{})
refs={'GG':(G,np.flatnonzero(ph0=='Tenericutes')),'SILVA':(S_,np.flatnonzero((S_.phylum=='Tenericutes')&(S_.domain=='Bacteria')))}
print({k:len(v[1]) for k,v in refs.items()})
def site_cov_ten(s,ref,rows):
    h=np.zeros(len(rows),bool); valid=np.ones(len(rows),bool)
    for p,e in zip(s['prim'],s['ext']):
        st,Ln=window(s['orient'],s['start'],s['L'],e); M=pdz.site_matrix(ref,st,Ln)[rows]; ok=(M!=0).all(1)&(M!=45).all(1); valid&=ok if e==0 else ok
        h|=ohit(ref,rows,s['orient'],st,Ln,top_of(s['orient'],p))
    return h,valid
D1=json.load(open(SC+'short_sites.json'))['D1']; orig=json.load(open(SC+'a5f_final5.json'))
sites={**{k:v for k,v in D1.items()}, 'A4-F(原 906–927)':orig['A4-F'],'A4-R(原)':orig['A4-R']}
five=[("R1-F","TGGCGAACGGGTGAGTAA","F",103),("R1-R","CCGTGTCTCAGTCCCARTG","R",314),("R2-F","ACTCCTACGGGAGGCAGC","F",338),("R2-R","GTATTACCGCGGCTGCTG","R",519),("R3-F","GTGTAGCGGTGRAATGCG","F",685),("R3-R","CCCGTCAATTCMTTTGAGTT","R",908),("R4-F","GGAGCATGTGGWTTAATTCGA","F",944),("R4-R","CGTTGCGGGACTTAACCC","R",1087),("R5-F","GGAGGAAGGTGGGGATGAC","F",1175),("R5-R","AAGGCCCGGGAACGTATT","R",1374)]
sites5R={f'5R {n}':dict(orient=o,start=(p if o=='F' else p-0),L=len(s),prim=[s],ext=[0]) for n,s,o,p in five}
# 5R 的 F 起点：p 为起点；R 的 p 为 3′ 端起点（顶链左端）？按 iterate_5R.FIVE：(name,seq,orient,pos) —— pos 对 F 是起点，对 R 是位点起点
rows=[]; masks={}
for nm,s in {**sites,**sites5R}.items():
    r=dict(site=nm,oligos=len(s['prim']))
    for db,(ref,rr) in refs.items():
        h,valid=site_cov_ten(s,ref,rr); r[db+'_n_valid']=int(valid.sum()); r[db+'_cov']=round(float(h[valid].mean()),3) if valid.sum() else np.nan; masks[(nm,db)]=(h,valid)
    rows.append(r)
R=pd.DataFrame(rows); print(R.to_string()); R.to_csv(SC+'ten_site_cov.csv',index=False)
pickle.dump((masks,gen0,fam0,refs[ 'GG'][1]),open(SC+'ten_masks.pkl','wb'))
