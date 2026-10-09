"""用已有的 480 个管的精确结果，找出决定分级的“热点引物对”（按引物名，不是展开序列）。
目的：第一阶段/第三阶段只算热点对，使 960 条序列的评估量下降一个数量级；本脚本同时给出这种简化与精确分级的一致率。
用法：python3 udp960_hotpairs.py <输出 json> <pkl...>"""
import sys, pickle, json, itertools, collections
import numpy as np
out = sys.argv[1]; pkls = sys.argv[2:]
tubes = {}
for p in pkls:
    r = pickle.load(open(p, 'rb'))
    for d in r['tubes'].values(): tubes[int(d['name'][2:])] = d
EX = r['EX']; N = len(EX); names = sorted({x[0] for x in EX}); nid = {n: i for i, n in enumerate(names)}
grp = np.array([nid[x[0]] for x in EX]); G = len(names)
ks = sorted(tubes)
def tier(h, e): return 0 if (h >= -7 and e >= -4.5) else 1 if (h >= -8.5 and e >= -5.5) else 2
# 每个管、每个名称对（含自身）的最小 H、E
pairs = [(a, b) for a in range(G) for b in range(a, G)]
PH = np.zeros((len(ks), len(pairs))); PE = np.zeros((len(ks), len(pairs)))
idx = [np.where(grp == g)[0] for g in range(G)]
for t, k in enumerate(ks):
    H, E = tubes[k]['H'], tubes[k]['E']; Es = np.minimum(E, E.T)
    for q, (a, b) in enumerate(pairs):
        ia, ib = idx[a], idx[b]; PH[t, q] = H[np.ix_(ia, ib)].min(); PE[t, q] = Es[np.ix_(ia, ib)].min()
full = np.array([tier(PH[t].min(), PE[t].min()) for t in range(len(ks))])
print('精确分级（A,B,C）', [(full == i).sum() for i in range(3)])
# 贪心：逐个加入名称对，使“只看已选名称对的分级”与精确分级一致的管最多
sel = []; cur_h = np.zeros(len(ks)); cur_e = np.zeros(len(ks)); best_acc = 0
for step in range(60):
    bestq, bestscore = None, -1
    for q in range(len(pairs)):
        if q in sel: continue
        h = np.minimum(cur_h, PH[:, q]); e = np.minimum(cur_e, PE[:, q])
        tt = np.array([tier(h[t], e[t]) for t in range(len(ks))]); score = (tt == full).mean() - 1e-4 * (tt > full).sum() / len(ks)
        if score > bestscore: bestq, bestscore, besttt = q, score, tt
    sel.append(bestq); cur_h = np.minimum(cur_h, PH[:, bestq]); cur_e = np.minimum(cur_e, PE[:, bestq])
    acc = (besttt == full).mean(); print(step + 1, names[pairs[bestq][0]], names[pairs[bestq][1]], '一致率 %.3f' % acc, flush=True)
    if acc >= 0.995: break
hot = [(names[pairs[q][0]], names[pairs[q][1]]) for q in sel]
json.dump(dict(hot=hot, acc=float(acc), n=len(sel), total_pairs=len(pairs)), open(out, 'w'), ensure_ascii=False, indent=1)
