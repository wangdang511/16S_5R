# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/
import pickle, json, numpy as np, pandas as pd, warnings, itertools, collections
warnings.filterwarnings('ignore'); pd.set_option('display.width',250)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'ten1.py').read().split("rows=[]; masks={}")[0]
exec(src)
S_big=pickle.load(open(SC+'S.pkl','rb'))
core=pdz.expand('YCTACGGGWGGCAGCA')
# 各主要门在 338–340 三个位置的碱基分布（该 16 nt 核心命中的序列）
for vn,(ref,rows_) in VAL.items():
    phs=ref.phylum[rows_]; s16=dict(orient='F',start=341,L=16,prim=['YCTACGGGWGGCAGCA'],ext=[0])
    h=site_hit(ref,rows_,s16)
    print(vn)
    for p in pdz.KEY_PHYLA:
        m=h&(phs==p)
        if m.sum()<10: continue
        cols=[]
        for c in (338,339,340):
            col=ref.aln[rows_[m],c-1]; cnt=collections.Counter(chr(x) if x in b'ACGT' else '-' for x in col); tot=sum(cnt.values()); cols.append(''.join(f'{b}{cnt[b]/tot:.2f} ' for b,_ in cnt.most_common(2)))
        print('  %-16s n=%4d'%(p,m.sum()),' | '.join(cols))
