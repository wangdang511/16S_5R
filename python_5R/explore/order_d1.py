# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/v67_* 与 order_list_5amp_v67short.csv
import pickle, json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore'); pd.set_option('display.width',250); pd.set_option('display.max_colwidth',60)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'ext2.py').read().split("rows=[]")[0]
exec(src)
sites=json.load(open(SC+'short_sites.json'))['D1']
AMP={'A1':'V1·V2','A2':'V3','A3':'V4','A4':'V6·V7','A5':'V8·V9'}
rows=[]; GC=lambda s:sum(c in 'GCS' for c in s)/len(s)
for sn in sorted(sites):
    s=sites[sn]; o=s['orient']
    # 该位点各寡核苷酸命中（SILVA+GG 验证集，主要门等权）
    hits=[]
    for p,e in zip(s['prim'],s['ext']):
        one=dict(s,prim=[p],ext=[e]); per=[]
        for vn,(ref,rows_) in VAL.items():
            phs=ref.phylum[rows_]; h=site_hit(ref,rows_,one)
            per.append(np.mean([h[phs==g].mean() for g in pdz.KEY_PHYLA if (phs==g).sum()>=30]))
        hits.append(float(np.mean(per)))
    tot=sum(hits)
    for i,(p,e) in enumerate(zip(s['prim'],s['ext'])):
        a,b=(s['start']-e,s['start']+s['L']-1) if o=='F' else (s['start'],s['start']+s['L']-1+e)
        tm=cp.oligo_tm(p); ts=[cp.oligo_tm(x)[0] for x in pdz.expand(p)]; lo,hi=min(ts),max(ts)
        rows.append(dict(site=sn,amplicon=AMP[sn[:2]],orient=o,idx=i,n_at_site=len(s['prim']),seq=p,nt=len(p),pos=f'{a}–{b}',Tm=round(tm[0],1),Tm_lo=round(lo,1),Tm_hi=round(hi,1),GC=round(GC(p),2),expansions=len(pdz.expand(p)),cov_alone=round(hits[i],3),share_raw=hits[i]/tot))
df=pd.DataFrame(rows)
# 建议比例：同一位点内按各条单独覆盖率归一化，取整到 5%，最后一条补足
sh=[]
for sn,g in df.groupby('site',sort=False):
    if len(g)==1: sh+=[1.0]; continue
    r=(g.share_raw*20).round()/20; r.iloc[-1]=round(1-r.iloc[:-1].sum(),2); sh+=list(r)
df['share']=sh
# 扩增子范围
amps={}
for k in AMP:
    f=df[(df.site==k+'-F')].pos.str.split('–').str[0].astype(int).min(); r=df[(df.site==k+'-R')].pos.str.split('–').str[1].astype(int).max(); amps[k]=(f,r)
df['amplicon_span']=[f'{amps[x[:2]][0]}–{amps[x[:2]][1]} ({amps[x[:2]][1]-amps[x[:2]][0]+1} bp)' for x in df.site]
# 10× 混合液：每个位点 2.5 µM（终浓度 250 nM），100 µM 母液
df['uL_100uM_per_100uL_10x']=(df.share*2.5).round(2)
df['name']=[f"SMURF5-{r.site}{'' if r.n_at_site==1 else '.'+str(r.idx+1)}" for r in df.itertuples()]
df=df.drop(columns=['share_raw'])
df.to_csv(SC+'order_list_5amp_v67short.csv',index=False)
print(df[['name','amplicon','orient','seq','nt','pos','Tm','Tm_lo','Tm_hi','GC','expansions','share','uL_100uM_per_100uL_10x','amplicon_span']].to_string())
sev,hp=cp.pool_issues({r.name:[r.seq] for r in df.itertuples()}); print('pool severe',sev,'hairpins',hp)
print('total stock µL',df.uL_100uM_per_100uL_10x.sum())
