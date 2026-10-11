"""把自建 UDP 从 384 对扩展到 480 对：保留已有的 384 对（已评估），追加 96 个新 i5 和 96 个新 i7（分 2 个 48 组，每组每个周期配色平衡）。
新序列与同类已有序列：编辑距离 ≥4、与反向互补 ≥3（池级约束；选 48/96 组时再要求所有序列两两 ≥4 且反向互补 ≥4）。
用法：python3 udp_denovo_480.py <dn384_final.json> <ext i5 pkl> <ext i7 pkl> <输出 json> <seed>"""
import sys, pickle, json, random, math
import numpy as np
from rapidfuzz.process import cdist
from rapidfuzz.distance import Levenshtein
fin, p5, p7, out, seed = sys.argv[1:6]; seed = int(seed)
rc = lambda s: s[::-1].translate(str.maketrans('ACGT', 'TGCA'))
d = json.load(open(fin)); old5, old7 = d['i5'], d['i7']
P5 = pickle.load(open(p5, 'rb')); P7 = pickle.load(open(p7, 'rb'))
B = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
def conf_count(P, others):
    """P 中每条序列与 others 的冲突数（编辑距离 <4 或 与反向互补 <4）"""
    D = cdist(P, others, scorer=Levenshtein.distance, dtype=np.uint8, workers=2); R = cdist(P, [rc(s) for s in others], scorer=Levenshtein.distance, dtype=np.uint8, workers=2)
    R2 = cdist([rc(s) for s in P], others, scorer=Levenshtein.distance, dtype=np.uint8, workers=2)
    return ((D < 4) | (R < 4) | (R2 < 4)).sum(1)
def pen(cnt, N=48):
    gc = cnt[:, 1] + cnt[:, 2]; red = cnt[:, 0] + cnt[:, 1]; green = cnt[:, 0] + cnt[:, 3]; v = 0.0
    for x in (gc, red, green): v += (np.maximum(0, 0.4 * N - x) ** 2 + np.maximum(0, x - 0.6 * N) ** 2).sum()
    return v + (np.maximum(0, 0.15 * N - cnt) ** 2).sum()
def select(P, cost, sd, iters=400000):
    n = len(P); M = np.array([[B[c] for c in s] for s in P]); rnd = random.Random(sd); order = list(range(n)); rnd.shuffle(order)
    slot = order[:96]; unused = order[96:]; cnt = np.zeros((2, 10, 4), int)
    for k, s in enumerate(slot):
        for pos in range(10): cnt[k // 48, pos, M[s, pos]] += 1
    gp = np.array([pen(cnt[g]) for g in range(2)]); tot = 1000 * gp.sum() + cost[slot].sum(); best = (tot, list(slot))
    for it in range(iters):
        T = 40.0 * (1 - it / iters) + 0.02; i = rnd.randrange(96); gi = i // 48
        if rnd.random() < 0.5 and unused:
            j = rnd.randrange(len(unused)); a, b = slot[i], unused[j]
            for pos in range(10): cnt[gi, pos, M[a, pos]] -= 1; cnt[gi, pos, M[b, pos]] += 1
            ng = pen(cnt[gi]); dE = (ng - gp[gi]) * 1000 + cost[b] - cost[a]
            if dE <= 0 or rnd.random() < math.exp(-dE / T): slot[i] = b; unused[j] = a; gp[gi] = ng; tot += dE
            else:
                for pos in range(10): cnt[gi, pos, M[a, pos]] += 1; cnt[gi, pos, M[b, pos]] -= 1
        else:
            k = rnd.randrange(96); gk = k // 48
            if gk == gi: continue
            a, b = slot[i], slot[k]
            for pos in range(10): cnt[gi, pos, M[a, pos]] -= 1; cnt[gi, pos, M[b, pos]] += 1; cnt[gk, pos, M[b, pos]] -= 1; cnt[gk, pos, M[a, pos]] += 1
            n1, n2 = pen(cnt[gi]), pen(cnt[gk]); dE = (n1 + n2 - gp[gi] - gp[gk]) * 1000
            if dE <= 0 or rnd.random() < math.exp(-dE / T): slot[i], slot[k] = b, a; gp[gi], gp[gk] = n1, n2; tot += dE
            else:
                for pos in range(10): cnt[gi, pos, M[a, pos]] += 1; cnt[gi, pos, M[b, pos]] -= 1; cnt[gk, pos, M[b, pos]] += 1; cnt[gk, pos, M[a, pos]] -= 1
        if tot < best[0]: best = (tot, list(slot))
    s = best[1]; cnt = np.zeros((2, 10, 4), int)
    for k, x in enumerate(s):
        for pos in range(10): cnt[k // 48, pos, M[x, pos]] += 1
    return [P[x] for x in s], float(sum(pen(cnt[g]) for g in range(2)))
c5 = conf_count(P5, old7).astype(float); n5, p5v = select(P5, c5, seed)
print('新 i5：周期越界', p5v, '与已有 i7 的冲突数合计', int(conf_count(n5, old7).sum()), flush=True)
used = set(old5 + n5); P7 = [x for x in P7 if x not in used and rc(x) not in used]
c7 = conf_count(P7, old5 + n5).astype(float); n7, p7v = select(P7, c7, seed + 1)
print('新 i7：周期越界', p7v, '与已有 i5+新 i5 的冲突数合计', int(conf_count(n7, old5 + n5).sum()), flush=True)
json.dump(dict(i5=old5 + n5, i7=old7 + n7, new5=n5, new7=n7, pen5=p5v, pen7=p7v), open(out, 'w'))
