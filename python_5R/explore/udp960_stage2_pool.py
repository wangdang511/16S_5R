"""第二阶段：从 960 条尾巴（U001–U960，每条都可当 i5 或 i7）中选候选池：M 个 i5 + M 个 i7。
硬约束：所选 2M 条序列两两 Levenshtein ≥4；i5 与 i7 之间（跨类）与反向互补的编辑距离 ≥3（x vs rc(y) 与 rc(x) vs y 都要）；
        i5 池、i7 池各自每个周期 GC、A+C、A+T 占比 40–60% 且 A/C/G/T 各 ≥15%。
代价：第一阶段单个序列的分级（i5 角色看 FF 与 FRu，i7 角色看 RR 与 RFu；A=0、B=3、C=30，另加超出阈值的量和发夹）。
用法：python3 udp960_stage2_pool.py <stage1 pkl> <tails tsv> <out json> <M> <seed> <iters>"""
import sys, os, json, pickle, random, math
import numpy as np
from rapidfuzz.process import cdist
from rapidfuzz.distance import Levenshtein
s1, tsv, out, M, seed, iters = sys.argv[1:7]; M, seed, iters = int(M), int(seed), int(iters)
rc = lambda s: s[::-1].translate(str.maketrans('ACGT', 'TGCA'))
rows = [l.split() for l in open(tsv)]; names = [r[0] for r in rows]; seqs = [r[1] for r in rows]; n = len(seqs)
S1 = pickle.load(open(s1, 'rb'))
D = cdist(seqs, seqs, scorer=Levenshtein.distance, dtype=np.uint8, workers=2); Dr = cdist(seqs, [rc(s) for s in seqs], scorer=Levenshtein.distance, dtype=np.uint8, workers=2)
CD = (D < 4).astype(np.int16); np.fill_diagonal(CD, 0)
CR = ((Dr < 3) | (Dr.T < 3)).astype(np.int16); np.fill_diagonal(CR, 0)
print('冲突密度：编辑距离 %.4f，跨类反向互补 %.4f' % (CD.mean(), CR.mean()), flush=True)
def tier_cost(h, e):
    t = 0.0 if (h >= -7 and e >= -4.5) else 3.0 if (h >= -8.5 and e >= -5.5) else 30.0
    return t + 0.2 * max(0, -(h + 7)) + 0.2 * max(0, -(e + 4.5))
c5 = np.zeros(n); c7 = np.zeros(n)
for i, nm in enumerate(names):
    r = S1[nm]
    c5[i] = tier_cost(min(r['FF'][0], r['FRu'][0]), min(r['FF'][1], r['FRu'][1])) + 0.5 * max(0, -(r['hpF'] + 2))
    c7[i] = tier_cost(min(r['RR'][0], r['RFu'][0]), min(r['RR'][1], r['RFu'][1])) + 0.5 * max(0, -(r['hpR'] + 2))
B = {'A': 0, 'C': 1, 'G': 2, 'T': 3}; MM = np.array([[B[c] for c in s] for s in seqs])
def pen_cnt(cnt, N):
    gc = cnt[:, 1] + cnt[:, 2]; red = cnt[:, 0] + cnt[:, 1]; green = cnt[:, 0] + cnt[:, 3]; v = 0.0
    for x in (gc, red, green): v += (np.maximum(0, 0.4 * N - x) ** 2 + np.maximum(0, x - 0.6 * N) ** 2).sum()
    return v + (np.maximum(0, 0.15 * N - cnt) ** 2).sum()
rnd = random.Random(seed)
order = list(range(n)); rnd.shuffle(order)
slot5 = order[:M]; slot7 = order[M:2 * M]; used = np.zeros(n, bool); used[slot5 + slot7] = True
cnt5 = np.zeros((10, 4), int); cnt7 = np.zeros((10, 4), int)
for s in slot5:
    for p in range(10): cnt5[p, MM[s, p]] += 1
for s in slot7:
    for p in range(10): cnt7[p, MM[s, p]] += 1
cD = CD[:, slot5 + slot7].sum(1).astype(int); cR5 = CR[:, slot5].sum(1).astype(int); cR7 = CR[:, slot7].sum(1).astype(int)
allsel = slot5 + slot7
conf = int(sum(cD[s] for s in allsel) / 2 + sum(cR7[s] for s in slot5))
pen5 = pen_cnt(cnt5, M); pen7 = pen_cnt(cnt7, M); cost = c5[slot5].sum() + c7[slot7].sum()
E = 1000 * conf + 1000 * (pen5 + pen7) + cost; best = (E, list(slot5), list(slot7)); print('起始', conf, pen5 + pen7, flush=True)
for it in range(iters):
    if it % 300000 == 0 and it: json.dump(dict(it=it, E=best[0], i5=[names[x] for x in best[1]], i7=[names[x] for x in best[2]]), open(out + '.part', 'w'))
    Tm = 80.0 * (1 - it / iters) + 0.05; role = rnd.randrange(2); b = rnd.randrange(n)
    if used[b]: continue
    sl = slot5 if role == 0 else slot7; i = rnd.randrange(M); a = sl[i]
    cO = cR7 if role == 0 else cR5            # 另一角色的反向互补冲突计数（b、a 与对方角色所选序列）
    dconf = (cD[b] - CD[b, a] + cO[b]) - (cD[a] + cO[a])
    cnt = cnt5 if role == 0 else cnt7
    for p in range(10): cnt[p, MM[a, p]] -= 1; cnt[p, MM[b, p]] += 1
    cc = c5 if role == 0 else c7
    newpen = pen_cnt(cnt, M); oldpen = pen5 if role == 0 else pen7
    dE = 1000 * dconf + 1000 * (newpen - oldpen) + cc[b] - cc[a]
    if dE <= 0 or rnd.random() < math.exp(-dE / Tm):
        sl[i] = b; used[a] = False; used[b] = True; cD += CD[:, b] - CD[:, a]
        if role == 0: cR5 += CR[:, b] - CR[:, a]; pen5 = newpen
        else: cR7 += CR[:, b] - CR[:, a]; pen7 = newpen
        conf += dconf; E += dE
        if E < best[0]: best = (E, list(slot5), list(slot7))
    else:
        for p in range(10): cnt[p, MM[a, p]] += 1; cnt[p, MM[b, p]] -= 1
E_, s5, s7 = best
sel = s5 + s7; cdd = int(CD[np.ix_(sel, sel)].sum() / 2); crr = int(CR[np.ix_(s5, s7)].sum())
cn5 = np.zeros((10, 4), int); cn7 = np.zeros((10, 4), int)
for s in s5:
    for p in range(10): cn5[p, MM[s, p]] += 1
for s in s7:
    for p in range(10): cn7[p, MM[s, p]] += 1
print('最终：编辑距离冲突', cdd, '跨类反向互补冲突', crr, '周期越界', pen_cnt(cn5, M) + pen_cnt(cn7, M), 'i5 平均代价 %.2f i7 平均代价 %.2f' % (c5[s5].mean(), c7[s7].mean()), flush=True)
json.dump(dict(i5=[names[x] for x in s5], i7=[names[x] for x in s7], seq={names[x]: seqs[x] for x in sel}, c5=[float(c5[x]) for x in s5], c7=[float(c7[x]) for x in s7]), open(out, 'w'))
