# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/v3_*
import sys; sys.path.insert(0,'/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad')
from v2r_scan import *
from short_scan import work
if __name__=='__main__':
    tasks=[]
    for L in (19,20):
        tasks+=[('A4F',p,L) for p in range(304,347)]   # 复用 work：F 用 3′ 端位置
        tasks+=[('A4R',p,L) for p in range(505,517)]
    res=[]
    with Pool(4) as p:
        for r in p.imap_unordered(work,tasks): res.extend(r)
    df=pd.DataFrame(res); df['tmmax']=df.tm.apply(max); df['tmmin']=df.tm.apply(min); df.to_pickle(SC+'v3_scan.pkl'); print(len(df))
