"""第四阶段（b）：在完整 144×144 里同时选 96 个 i5、96 个 i7 并配对（前 48 对 = 核心 48 组，是 96 组的子集）。
管的预测：H = min(FF(i5), RR(i7), FR(i5,i7))，E 同理；FF/RR 取第一阶段与增量（*_extra）的最小值，FR 取第三阶段与增量的最小值。
硬约束：96 个 i5、96 个 i7 各自每个周期配色达标；核心 48 个也各自达标；池内序列已满足编辑距离 ≥4、跨类反向互补 ≥3。
目标：A 级尽量多（核心 48 的权重 ×2）。
用法：python3 udp960_stage4b_pair.py <pool144> <stage1> <stage1_extra> <stage3> <stage3_extra> <tails tsv> <out json> <seed> <iters>"""
import sys, json, pickle, random, math
import numpy as np
poolj, s1p, s1xp, s3p, s3xp, tsv, out, seed, iters = sys.argv[1:10]; seed, iters = int(seed), int(iters)
pool = json.load(open(poolj)); S1 = pickle.load(open(s1p, 'rb')); S1x = pickle.load(open(s1xp, 'rb')); st = pickle.load(open(s3p, 'rb')); stx = pickle.load(open(s3xp, 'rb'))
seqs = {l.split()[0]: l.split()[1] for l in open(tsv)}
names5, names7 = pool['i5'], pool['i7']; nr, nc = len(names5), len(names7)
HM = np.zeros((nr, nc)); EM = np.zeros((nr, nc))
for i in range(nr):
    h0, e0, _ = st['rows'][i]; h1, e1, _ = stx['rows'][i]
    HM[i] = np.minimum(h0, h1); EM[i] = np.minimum(e0, e1)
assert not np.isnan(HM).any()
mn = lambda a, b: (min(a[0], b[0]), min(a[1], b[1]))
FF = [mn(S1[n]['FF'], S1x[n]['FF']) for n in names5]; RR = [mn(S1[n]['RR'], S1x[n]['RR']) for n in names7]
H5 = np.array([x[0] for x in FF]); E5 = np.array([x[1] for x in FF]); H7 = np.array([x[0] for x in RR]); E7 = np.array([x[1] for x in RR])
TH = np.minimum(HM, np.minimum(H5[:, None], H7[None, :])); TE = np.minimum(EM, np.minimum(E5[:, None], E7[None, :]))
def tier_cost(h, e):
    t = np.where((h >= -7) & (e >= -4.5), 0.0, np.where((h >= -8.5) & (e >= -5.5), 3.0, 30.0))
    return t + 0.2 * np.maximum(0, -(h + 7)) + 0.2 * np.maximum(0, -(e + 4.5))
COST = tier_cost(TH, TE); TIER = np.where((TH >= -7) & (TE >= -4.5), 'A', np.where((TH >= -8.5) & (TE >= -5.5), 'B', 'C'))
print('144×144 组合分级（预测）：', {t: int((TIER == t).sum()) for t in 'ABC'}, '占比 A %.1f%%' % (100 * (TIER == 'A').mean()), flush=True)
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
def energy(s5, s7): return (COST[s5, s7] * W).sum() + 1000 * (pen(M5, s5, 96) + pen(M7, s7, 96) + pen(M5, s5[:48], 48) + pen(M7, s7[:48], 48))
E = energy(s5, s7); best = (E, list(s5), list(s7))
for it in range(iters):
    Tm = 40.0 * (1 - it / iters) + 0.05; mv = rnd.random(); n5, n7 = list(s5), list(s7); a5, a7 = list(u5), list(u7)
    if mv < 0.3: i, j = rnd.sample(range(96), 2); n7[i], n7[j] = n7[j], n7[i]
    elif mv < 0.5: i, j = rnd.sample(range(96), 2); n5[i], n5[j] = n5[j], n5[i]
    elif mv < 0.75: i = rnd.randrange(96); k = rnd.randrange(len(a5)); n5[i], a5[k] = a5[k], n5[i]
    else: i = rnd.randrange(96); k = rnd.randrange(len(a7)); n7[i], a7[k] = a7[k], n7[i]
    E2 = energy(n5, n7)
    if E2 <= E or rnd.random() < math.exp((E - E2) / Tm):
        s5, s7, u5, u7, E = n5, n7, a5, a7, E2
        if E < best[0]: best = (E, list(s5), list(s7))
    if it % 100000 == 0 and it: print(it, round(E, 1), round(best[0], 1), flush=True)
E, s5, s7 = best; tiers = [TIER[a, b] for a, b in zip(s5, s7)]
print('最终：96 对分级（预测）', {t: tiers.count(t) for t in 'ABC'}, '；核心 48 对', {t: tiers[:48].count(t) for t in 'ABC'}, '；周期越界 96:', pen(M5, s5, 96) + pen(M7, s7, 96), '核心 48:', pen(M5, s5[:48], 48) + pen(M7, s7[:48], 48))
json.dump(dict(pairs=[dict(i5=names5[a], i7=names7[b], tier=TIER[a, b], H=float(TH[a, b]), E=float(TE[a, b]), seq5=seqs[names5[a]], seq7=seqs[names7[b]]) for a, b in zip(s5, s7)]), open(out, 'w'), ensure_ascii=False)
