# 注：此脚本在本次会话中从 scratchpad 运行（借用 Olivar 的 dG/GC/复杂度/自比对函数）；结果见 docs/primer_design/layout6_olivar_primer_checks.csv 与 layout6_amplicon_gc.csv
import sys, pickle, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0,'/home/user/16S_5R/python_5R'); sys.path.insert(0,'/home/user/wangdang511/olivar_primer/src/olivar')
import primer_design as pdz, silva_ref as sr, design as od, basic
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
# 1) 不同退火温度下的 Olivar dG 通过率
T=pd.read_pickle(SC+'pool32v2_table_off.pkl'); out=[]
for temp in (50,55,58,60):
    g=od.primer_generator(temperature=temp,salinity=0.18); n=0; tot=0
    for q in T['5′→3′序列']:
        for e in pdz.expand(q): tot+=1; n+= (g.dG_init+g.StacksDG(e.lower()))<=-11.8
    out.append((temp,n,tot)); print(f'退火 {temp} °C：dG ≤ -11.8 的展开序列 {n}/{tot} ({100*n/tot:.0f}%)',flush=True)
# 2) 扩增子 GC 含量
S=pickle.load(open(SC+'S.pkl','rb')); p2c=np.load(SC+'lay_pos2col.npy')
ids,A,total=sr.sample_alignment(SC+'cand/silva128/share/sepp/ref/99_otus_aligned_masked1977.fasta',60000,1)
AMPS={'A1':(9,264),'A2':(334,522),'A3':(518,792),'A4':(788,972),'A5':(967,1210),'A6':(1222,1504)}
isgc=(A==ord('G'))|(A==ord('C'))|(A==ord('g'))|(A==ord('c')); isb=np.isin(A,np.frombuffer(b'ACGTUacgtu',np.uint8))
ph=np.array(S.phylum,dtype=object); dom=np.array(S.domain)
rows=[]
for a,(f,r) in AMPS.items():
    c0,c1=p2c[f-1],p2c[r-1]; gc=isgc[:,c0:c1+1].sum(1); nb=isb[:,c0:c1+1].sum(1); ok=(nb>=0.9*(r-f+1)); frac=gc/np.maximum(nb,1)
    x=frac[ok&(dom=='Bacteria')]
    rows.append(dict(扩增子=a,n=int(len(x)),GC中位=round(float(np.median(x)),3),p5=round(float(np.percentile(x,5)),3),p95=round(float(np.percentile(x,95)),3),低于35=round(float((x<0.35).mean()*100),1),高于65=round(float((x>0.65).mean()*100),1)))
    for p in ('Actinobacteria','Tenericutes','Proteobacteria','Firmicutes'):
        m=ok&(ph==p)&(dom=='Bacteria')
        if m.sum()>=30: rows[-1][p+'中位GC']=round(float(np.median(frac[m])),3)
D=pd.DataFrame(rows); print(D.to_string(index=False)); D.to_csv(SC+'gc_amp.csv',index=False)
