# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/layout6_*
from v5d_scan import *
if __name__=='__main__':
    done=set()
    tasks=[(o,st,19) for o in 'FR' for st in range(8,1500,2)]
    res=[]
    with Pool(4) as p:
        for r in p.imap_unordered(work2,tasks): res.extend(r)
    pd.DataFrame(res).to_pickle(SC+'v5d_scan3.pkl'); print(len(res))
