# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/v67_* 与 order_list_5amp_v67short.csv
import pickle, json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore'); pd.set_option('display.width',250); pd.set_option('display.max_colwidth',60)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'ext2.py').read().split("rows=[]")[0]
exec(src)
S=json.load(open(SC+'short_sites.json')); orig=json.load(open(SC+'a5f_final5.json')); D1=S['D1']
def keycov(s):
    out={}
    for vn,(ref,rows_) in VAL.items():
        phs=ref.phylum[rows_]; h=site_hit(ref,rows_,s)
        out[vn]=(round(float(np.mean([h[phs==g].mean() for g in pdz.KEY_PHYLA if (phs==g).sum()>=30])),3),round(float(min([h[phs==g].mean() for g in pdz.KEY_PHYLA if (phs==g).sum()>=30])),3))
    return out
rows=[]
for nm,s in (('A4-F 原（906–927，22 nt）',orig['A4-F']),('A4-F 新（967–985，19 nt ×2）',D1['A4-F']),('A4-R 原（1177–1195）',orig['A4-R']),('A4-R 新（1177–1195）',D1['A4-R'])):
    kc=keycov(s); tm=[round(cp.oligo_tm(p)[0],1) for p in s['prim']]
    rows.append(dict(site=nm,oligos=len(s['prim']),prim=' | '.join(s['prim']),expansions=sum(len(pdz.expand(p)) for p in s['prim']),Tm=tm,SILVA_mean=kc['SILVA'][0],SILVA_min=kc['SILVA'][1],GG_mean=kc['GG'][0],GG_min=kc['GG'][1]))
R=pd.DataFrame(rows); print(R.to_string()); R.to_csv(SC+'v67_primer.csv',index=False)
# 位点保守度（主要门平均 / 最差）以及在设计集上的单条引物覆盖：同一预算下 旧 vs 新
D=pd.read_pickle(SC+'short_scan.pkl')
for kind,pos,Ls in (('A4F',927,(19,20)),('A4F',985,(19,)),('A4F',986,(19,)),('A4F',987,(19,))):
    g=D[(D.kind==kind)&(D.pos==pos)&(D.L.isin(Ls))&(D.severe==0)]
    for k in (1,2):
        h=g[g.k==k].sort_values('cov_min',ascending=False).head(1)
        for _,r in h.iterrows(): print(kind,pos,'L',r.L,'k',k,'cov_min %.3f'%r.cov_min,'nexp',r.nexp,r.tm,r.prim)
