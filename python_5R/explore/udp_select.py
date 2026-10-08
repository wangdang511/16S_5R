"""从 UDP0001–0384 中挑 48 / 96 组 i5/i7 内联索引。
硬约束：i5 之间、i7 之间的 Levenshtein 编辑距离 ≥4；每个索引 GC 4–6/10；无同聚 ≥3；无近回文；每个周期（第 1–10 位）上 GC 占比 40–60%，A+C 与 A+T（两通道信号）也在 40–60%，且 A/C/G/T 各 ≥15%；
        管内二聚体分级 A/B（C 级不用）。软目标：A 级尽量多、每位点碱基更均衡。
用法：python3 udp_select.py <评估 pkl> <候选 tsv> <输出 json> [cross_min]   cross_min = i5 与 i7 之间的最小编辑距离（0 表示不约束）"""
import sys, pickle, json, random, math, collections, itertools
import numpy as np
pkl, tsv, out = sys.argv[1:4]; CROSS = int(sys.argv[4]) if len(sys.argv) > 4 else 0
R = pickle.load(open(pkl, 'rb')); EX = R['EX']; N0 = len(EX)
def lev(a, b):
    p = list(range(len(b) + 1))
    for i, x in enumerate(a, 1):
        c = [i]
        for j, y in enumerate(b, 1): c.append(min(p[j] + 1, c[-1] + 1, p[j - 1] + (x != y)))
        p = c
    return p[-1]
iu = np.triu_indices(N0)
cand = []
for k, d in sorted(R['tubes'].items()):
    minH = float(d['H'][iu].min()); minE = float(d['E'].min())
    tier = 'A' if (minH >= -7 and minE >= -4.5) else 'B' if (minH >= -8.5 and minE >= -5.5) else 'C'
    cand.append(dict(name=d['name'], i5=d['i5'], i7=d['i7'], minH=minH, minE=minE, tier=tier, bad=float(d['bad'])))
cand = [c for c in cand if c['tier'] != 'C']; n = len(cand)
print('可用候选', n, collections.Counter(c['tier'] for c in cand))
L5 = np.zeros((n, n), int); L7 = np.zeros((n, n), int); X = np.zeros((n, n), int)
for a in range(n):
    for b in range(n):
        if a < b: L5[a, b] = L5[b, a] = lev(cand[a]['i5'], cand[b]['i5']); L7[a, b] = L7[b, a] = lev(cand[a]['i7'], cand[b]['i7'])
        X[a, b] = lev(cand[a]['i5'], cand[b]['i7'])
ok = (L5 >= 4) & (L7 >= 4)
if CROSS: ok &= (X >= CROSS) & (X.T >= CROSS)
np.fill_diagonal(ok, True)
B = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
M5 = np.array([[B[ch] for ch in c['i5']] for c in cand]); M7 = np.array([[B[ch] for ch in c['i7']] for c in cand])
def balance(S, N):
    """每个周期（位点）：GC 占比须在 40–60%，且 A/C/G/T 每种都 ≥15%（两通道化学里 A+C 与 A+T 也都在 40–60%）。pen = 越界量的平方和。"""
    lo_gc, hi_gc = 0.40 * N, 0.60 * N; lo_b = 0.15 * N; pen = 0.0; sse = 0.0
    for M in (M5, M7):
        for pos in range(10):
            cnt = np.bincount(M[S, pos], minlength=4)      # A C G T
            gc = cnt[1] + cnt[2]; red = cnt[0] + cnt[1]; green = cnt[0] + cnt[3]
            for v in (gc, red, green): pen += max(0, lo_gc - v) ** 2 + max(0, v - hi_gc) ** 2
            for v in cnt: pen += max(0, lo_b - v) ** 2
            sse += ((cnt - N / 4) ** 2).sum()
    return pen, sse
isB = np.array([c['tier'] == 'B' for c in cand])
def energy(S, N):
    S = list(S); sub = ok[np.ix_(S, S)]; conf = (~sub).sum() / 2
    pen, sse = balance(S, N)
    return 1000 * conf + 200 * pen + 3 * isB[S].sum() + 0.5 * sse, conf, pen, sse
def anneal(N, pool, seed, iters=60000):
    rnd = random.Random(seed); pool = list(pool); S = rnd.sample(pool, N); rest = [x for x in pool if x not in S]
    e = energy(S, N)[0]; T0 = 30.0
    best = (e, list(S))
    for it in range(iters):
        T = T0 * (1 - it / iters) + 0.05
        i = rnd.randrange(N); j = rnd.randrange(len(rest)); old = S[i]; S[i] = rest[j]
        e2 = energy(S, N)[0]
        if e2 <= e or rnd.random() < math.exp((e - e2) / T): rest[j] = old; e = e2
        else: S[i] = old
        if e < best[0]: best = (e, list(S))
    return best
res = {}
for N in (96, 48):
    pool = range(n); allb = []
    for seed in range(6):
        allb.append(anneal(N, pool, seed))
    e, S = min(allb, key=lambda x: x[0]); en, conf, pen, sse = energy(S, N)
    res[N] = dict(sel=[cand[i]['name'] for i in S], energy=e, conflicts=int(conf), balance_excess=float(pen), sse=float(sse), nB=int(isB[S].sum()))
    print(N, '能量 %.1f 冲突 %d 超容差 %.1f SSE %.1f B 级 %d' % (e, conf, pen, sse, isB[S].sum()), flush=True)
json.dump(dict(res=res, cand=cand), open(out, 'w'), ensure_ascii=False)
