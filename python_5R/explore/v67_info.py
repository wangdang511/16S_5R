# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/v67_* 与 order_list_5amp_v67short.csv
import pickle, json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore'); pd.set_option('display.width',250); pd.set_option('display.max_colwidth',50)
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
src=open(SC+'ext2.py').read().split("rows=[]")[0]
exec(src)
S=json.load(open(SC+'short_sites.json')); orig=json.load(open(SC+'a5f_final5.json')); D1=S['D1']
ints0=[(orig[f'A{k}-F']['start']+orig[f'A{k}-F']['L'],orig[f'A{k}-R']['start']-1) for k in range(1,6)]
print('原 A4 内部',ints0[3])
base=ts.accuracy_direct(xall,ints0,gall,vall)
def acc(ints): return ts.accuracy_direct(xall,ints,gall,vall)
rows=[]
for nm,iv in (('原 V6V7 内部 928–1176（249 bp）',(928,1176)),('缩短后 986–1176（191 bp）',(986,1176)),('被去掉的 928–985（58 bp）',(928,985)),('V6（986–1043）',(986,1043)),('V7（1117–1173）',(1117,1173)),('间隔区 1044–1116',(1044,1116))):
    rows.append(dict(region=nm,bp=iv[1]-iv[0]+1,alone=round(acc([iv]),4)))
# 整体设计中的边际贡献
full=acc(ints0); noA4=acc([i for k,i in enumerate(ints0) if k!=3]); short=acc([ints0[0],ints0[1],ints0[2],(986,1176),ints0[4]])
print('整体 5 扩增子',round(full,4),'去掉 V6V7',round(noA4,4),'V6V7 缩短',round(short,4))
R=pd.DataFrame(rows); print(R.to_string())
# 位点熵（带标签序列）
X=xall
def ent(c):
    cnt=np.stack([(c==b).sum(0) for b in range(4)]).astype(float); n=cnt.sum(0); p=cnt/np.maximum(n,1); 
    with np.errstate(all='ignore'): h=-(np.where(p>0,p*np.log2(p),0)).sum(0)
    return h
H=ent(X)
for nm,(a,b) in (('928–985',(928,985)),('986–1176',(986,1176)),('V6',(986,1043)),('V7',(1117,1173))):
    h=H[a-1:b]; print(nm,'平均熵 %.3f bit，变异位点（熵>0.5）占 %.0f%%，熵之和 %.1f bit'%(h.mean(),100*(h>0.5).mean(),h.sum()))
# 单位长度：准确率增益 per bp
print('每 bp 的边际准确率贡献（整体设计）：原 %.4f 个百分点/bp，缩短后 %.4f'%((full-noA4)*100/249,(short-noA4)*100/191))
pickle.dump(dict(full=full,noA4=noA4,short=short,R=R),open(SC+'v67_info.pkl','wb'))
R.to_csv(SC+'v67_info.csv',index=False)
