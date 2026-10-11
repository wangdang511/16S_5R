"""第四阶段（a）：在 100×100 里同时选 96 个 i5、96 个 i7 并配对（前 48 对是核心 48 组，是 96 组的子集）。
管的预测分级用热点代理：H = min(FF(i5), RR(i7), FR(i5,i7))，E 同理（FF、RR 来自第一阶段，FR 来自第三阶段）。
硬约束：96 个 i5、96 个 i7 各自每个周期配色达标；核心 48 个 i5、48 个 i7 也各自达标（N=48）。池内序列已满足编辑距离 ≥4、跨类反向互补 ≥3。
目标：A 级尽量多（核心 48 组的权重 ×2）。
用法：python3 udp960_stage4_pair.py <pool144 json> <subset100 json> <stage1 pkl> <stage3 pkl> <tails tsv> <out json> <seed> <iters>"""
import sys, json, pickle, random, math
import numpy as np
poolj, subj, s1p, s3p, tsv, out, seed, iters = sys.argv[1:9]; seed, iters = int(seed), int(iters)
pool = json.load(open(poolj)); sub = json.load(open(subj)); S1 = pickle.load(open(s1p, 'rb')); st = pickle.load(open(s3p, 'rb'))
seqs = {l.split()[0]: l.split()[1] for l in open(tsv)}
R = sub['rows']; C = sub['cols']; nr, nc = len(R), len(C)
names5 = [pool['i5'][i] for i in R]; names7 = [pool['i7'][j] for j in C]
HM = np.zeros((nr, nc)); EM = np.zeros((nr, nc))
for a, i in enumerate(R):
    hs, es, early = st['rows'][i]
    for b, j in enumerate(C): HM[a, b] = hs[j]; EM[a, b] = es[j]
assert not np.isnan(HM).any()
H5 = np.array([S1[n]['FF'][0] for n in names5]); E5 = np.array([S1[n]['FF'][1] for n in names5])
H7 = np.array([S1[n]['RR'][0] for n in names7]); E7 = np.array([S1[n]['RR'][1] for n in names7])
TH = np.minimum(HM, np.minimum(H5[:, None], H7[None, :])); TE = np.minimum(EM, np.minimum(E5[:, None], E7[None, :]))
def tier_cost(h, e):
    t = np.where((h >= -7) & (e >= -4.5), 0.0, np.where((h >= -8.5) & (e >= -5.5), 3.0, 30.0))
    return t + 0.2 * np.maximum(0, -(h + 7)) + 0.2 * np.maximum(0, -(e + 4.5))
COST = tier_cost(TH, TE); TIER = np.where((TH >= -7) & (TE >= -4.5), 'A', np.where((TH >= -8.5) & (TE >= -5.5), 'B', 'C'))
print('100×100 组合分级（预测）：', {t: int((TIER == t).sum()) for t in 'ABC'}, '占比 A %.1f%%' % (100 * (TIER == 'A').mean()), flush=True)
Bm = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
M5 = np.array([[Bm[c] for c in seqs[n]] for n in names5]); M7 = np.array([[Bm[c] for c in seqs[n]] for n in names7])
def pen(M, idx, N):
    p = 0.0
    for pos in range(10):
        c = np.bincount(M[idx, pos], minlength=4); gc = c[1] + c[2]; red = c[0] + c[1]; green = c[0] + c[3]
        for x in (gc, red, green): p += max(0, 0.4 * N - x) ** 2 + max(0, x - 0.6 * N) ** 2
        for x in c: p += max(0, 0.15 * N - x) ** 2
    return p
rnd = random.Random(seed)
p5 = list(range(nr)); rnd.shuffle(p5); p7 = list(range(nc)); rnd.shuffle(p7)
s5 = p5[:96]; u5 = p5[96:]; s7 = p7[:96]; u7 = p7[96:]
W = np.array([2.0] * 48 + [1.0] * 48)
def energy(s5, s7):
    pc = (COST[s5, s7] * W).sum()
    return pc + 1000 * (pen(M5, s5, 96) + pen(M7, s7, 96) + pen(M5, s5[:48], 48) + pen(M7, s7[:48], 48))
E = energy(s5, s7); best = (E, list(s5), list(s7)); print('起始', E, flush=True)
for it in range(iters):
    Tm = 40.0 * (1 - it / iters) + 0.05; mv = rnd.random(); n5, n7 = list(s5), list(s7); a5, a7 = list(u5), list(u7)
    if mv < 0.3:                                   # 交换两个位置的 i7（改配对）
        i, j = rnd.sample(range(96), 2); n7[i], n7[j] = n7[j], n7[i]
    elif mv < 0.5:                                 # 交换两个位置的 i5
        i, j = rnd.sample(range(96), 2); n5[i], n5[j] = n5[j], n5[i]
    elif mv < 0.75:                                # 用未选 i5 替换
        i = rnd.randrange(96); k = rnd.randrange(len(a5)); n5[i], a5[k] = a5[k], n5[i]
    else:                                          # 用未选 i7 替换
        i = rnd.randrange(96); k = rnd.randrange(len(a7)); n7[i], a7[k] = a7[k], n7[i]
    E2 = energy(n5, n7)
    if E2 <= E or rnd.random() < math.exp((E - E2) / Tm):
        s5, s7, u5, u7, E = n5, n7, a5, a7, E2
        if E < best[0]: best = (E, list(s5), list(s7))
    if it % 100000 == 0 and it: print(it, round(E, 1), round(best[0], 1), flush=True)
E, s5, s7 = best
tiers = [TIER[a, b] for a, b in zip(s5, s7)]
print('最终：96 对分级（预测）', {t: tiers.count(t) for t in 'ABC'}, '；核心 48 对', {t: tiers[:48].count(t) for t in 'ABC'},
      '；周期越界 96:', pen(M5, s5, 96) + pen(M7, s7, 96), '核心 48:', pen(M5, s5[:48], 48) + pen(M7, s7[:48], 48))
json.dump(dict(pairs=[dict(i5=names5[a], i7=names7[b], tier=TIER[a, b], H=float(TH[a, b]), E=float(TE[a, b]), seq5=seqs[names5[a]], seq7=seqs[names7[b]]) for a, b in zip(s5, s7)]), open(out, 'w'), ensure_ascii=False)
