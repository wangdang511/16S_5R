"""从头设计 UDP：从预筛候选中选 96 个 i5 + 96 个 i7（前 48 个是 48 组），满足
 - 192 条序列两两 Levenshtein ≥ DMIN，且任意两条序列与对方反向互补的编辑距离 ≥ RCMIN；
 - 每个周期 GC、A+C、A+T 占比 40–60%，A/C/G/T 各 ≥15%：48 组（前 48 个）与 96 组（全部 96 个）的 i5、i7 各自满足；
 - 预筛二聚体指标（H、E）越好越优。
用法：python3 udp_denovo_select.py <cand pkl> <输出 json> <DMIN> <RCMIN> <Hmin> <Emin> <seed> [iters]"""
import sys, json, pickle, random, math
import numpy as np
from rapidfuzz.process import cdist
from rapidfuzz.distance import Levenshtein
pkl, out, DMIN, RCMIN, HMIN, EMIN, seed = sys.argv[1:8]; DMIN, RCMIN, HMIN, EMIN, seed = int(DMIN), int(RCMIN), float(HMIN), float(EMIN), int(seed)
ITERS = int(sys.argv[8]) if len(sys.argv) > 8 else 300000
rc = lambda s: s[::-1].translate(str.maketrans('ACGT', 'TGCA'))
cand = [c for c in pickle.load(open(pkl, 'rb')) if c[1] >= HMIN and c[2] >= EMIN]
seqs = [c[0] for c in cand]; n = len(seqs); print('候选', n, flush=True)
D = cdist(seqs, seqs, scorer=Levenshtein.distance, dtype=np.uint8, workers=4)
Drc = cdist(seqs, [rc(s) for s in seqs], scorer=Levenshtein.distance, dtype=np.uint8, workers=4)
bad = (D < DMIN) | (Drc < RCMIN) | (cdist([rc(s) for s in seqs], seqs, scorer=Levenshtein.distance, dtype=np.uint8, workers=4) < RCMIN)
np.fill_diagonal(bad, False)
B = {'A': 0, 'C': 1, 'G': 2, 'T': 3}; M = np.array([[B[c] for c in s] for s in seqs])
cost = np.array([-(c[1] + 3) * 0.5 - (c[2] + 2) * 0.5 for c in cand]); cost = np.maximum(cost, 0)
def bal(idx, N):
    lo_gc, hi_gc, lo_b = 0.4 * N, 0.6 * N, 0.15 * N; pen = 0.0; sse = 0.0
    Mi = M[idx]
    for pos in range(10):
        cnt = np.bincount(Mi[:, pos], minlength=4); gc = cnt[1] + cnt[2]; red = cnt[0] + cnt[1]; green = cnt[0] + cnt[3]
        for v in (gc, red, green): pen += max(0, lo_gc - v) ** 2 + max(0, v - hi_gc) ** 2
        for v in cnt: pen += max(0, lo_b - v) ** 2
        sse += ((cnt - N / 4) ** 2).sum()
    return pen, sse
# 槽位：0-47 i5(48组)，48-95 i5(后48)，96-143 i7(48组)，144-191 i7(后48)
def energy(S):
    S = np.array(S); conf = bad[np.ix_(S, S)].sum() / 2
    if len(set(S.tolist())) < len(S): conf += 1000
    p = 0.0; q = 0.0
    for a, b in ((0, 48), (0, 96)):
        pa, sa = bal(S[a:b], b - a); p += pa; q += sa
        pa, sa = bal(S[96 + a:96 + b], b - a); p += pa; q += sa
    return 1000 * conf + 200 * p + 0.3 * q + cost[S].sum(), conf, p
rnd = random.Random(seed); S = rnd.sample(range(n), 192); e = energy(S)[0]; best = (e, list(S))
for it in range(ITERS):
    T = 20.0 * (1 - it / ITERS) + 0.05; i = rnd.randrange(192); j = rnd.randrange(n)
    if j in S: continue
    old = S[i]; S[i] = j; e2 = energy(S)[0]
    if e2 <= e or rnd.random() < math.exp((e - e2) / T): e = e2
    else: S[i] = old
    if e < best[0]: best = (e, list(S))
    if it % 20000 == 0: print(it, round(e, 1), round(best[0], 1), flush=True)
e, S = best; en, conf, pen = energy(S); print('最终 能量 %.1f 冲突 %d 周期越界 %.1f' % (en, conf, pen))
json.dump(dict(i5=[seqs[k] for k in S[:96]], i7=[seqs[k] for k in S[96:]], conflicts=int(conf), pen=float(pen), H=[cand[k][1] for k in S], E=[cand[k][2] for k in S]), open(out, 'w'))
