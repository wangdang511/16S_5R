# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/snap_*
import pickle, json, numpy as np, pandas as pd, warnings, collections, re
warnings.filterwarnings('ignore'); pd.set_option('display.width',250); pd.set_option('display.max_colwidth',60)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
exec(open(SC+'ext2.py').read().split("rows=[]")[0])
S_=pickle.load(open(SC+'S.pkl','rb')); ref=S_.ref_seq
if isinstance(ref,bytes): ref=ref.decode()
ref=ref.upper().replace('U','T'); print('ref len',len(ref))
prim={}; name=None
for l in open('/root/.claude/uploads/0a157df1-6424-5d79-bd19-c9f07612faab/aa48595b-Swift_16S_SNAP_primers_v2.fasta'):
    l=l.strip()
    if l.startswith('>'): name=l[1:]
    elif l: prim[name]=l.lstrip('^').upper()
IU={'A':'A','C':'C','G':'G','T':'T','R':'AG','Y':'CT','S':'GC','W':'AT','K':'GT','M':'AC','B':'CGT','D':'AGT','H':'ACT','V':'ACG','N':'ACGT'}
def mm(p,w): return sum(1 for a,b in zip(p,w) if b not in IU[a])
rows=[]
for n,p in prim.items():
    f=n.endswith('_f') or '_f' in n
    best=None
    for strand,q in (('+',p),('-',pdz.revcomp(p))):
        for i in range(len(ref)-len(q)+1):
            w=ref[i:i+len(q)]; m=sum(1 for a,b in zip(q,w) if b not in IU[a]) if strand=='+' else sum(1 for a,b in zip(q,w) if b not in IU[a])
            if best is None or m<best[0]: best=(m,strand,i+1,i+len(q))
    rows.append(dict(name=n,seq=p,nt=len(p),mismatch=best[0],strand=best[1],start=best[2],end=best[3]))
R=pd.DataFrame(rows); print(R.to_string()); R.to_csv(SC+'snap_positions.csv',index=False)
