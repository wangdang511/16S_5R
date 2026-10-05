import pickle, json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'ext2.py').read().split("rows=[]")[0]
exec(src)
def ent(c):
    cnt=np.stack([(c==b).sum(0) for b in range(4)]).astype(float); n=cnt.sum(0); p=cnt/np.maximum(n,1)
    with np.errstate(all='ignore'): return -(np.where(p>0,p*np.log2(p),0)).sum(0)
H=ent(xall)
def s(a,b): h=H[a-1:b]; return 'n=%d 平均熵 %.2f bit，变异位点（熵>0.5）%.0f%%，最大 %.2f'%(b-a+1,h.mean(),100*(h>0.5).mean(),h.max())
print('缺口 1046–1055:',s(1046,1055)); print('V6 986–1043:',s(986,1043)); print('间隔 1044–1116:',s(1044,1116)); print('V7 1117–1173:',s(1117,1173)); print('V8V9 缺口 1363–1370:',s(1363,1370)); print('V8 1243–1294:',s(1243,1294)); print('V9 1435–1465:',s(1435,1465))
print('1044–1060 逐位熵:',[round(float(x),2) for x in H[1043:1060]])
print('1046–1055 逐位熵:',[round(float(x),2) for x in H[1045:1055]])
