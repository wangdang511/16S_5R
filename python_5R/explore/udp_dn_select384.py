"""从自建 384 对里选 48 / 96 组：管内二聚体尽量 A 级；所选集合可一起上机 pooling：
 i5∪i7 两两编辑距离 ≥4、与对方反向互补编辑距离 ≥4；每个周期 GC、A+C、A+T 占比 40–60% 且 A/C/G/T 各 ≥15%（i5、i7 各自）。
用法：python3 udp_dn_select384.py <out json> <conflict npy> <seed> <pkl1> [pkl2 ...]"""
import sys, pickle, json, random, math
import numpy as np
out, confl, seed = sys.argv[1:4]; seed = int(seed); pkls = sys.argv[4:]
T = np.load(confl); n = T.shape[0]
tubes = {}
for p in pkls:
    r = pickle.load(open(p, 'rb'))
    for d in r['tubes'].values(): tubes[int(d['name'][2:])] = d
    N0 = len(r['EX'])
assert len(tubes) == n, len(tubes)
iu = np.triu_indices(N0)
tier = []; mh = []; me = []; i5 = []; i7 = []
for k in range(1, n + 1):
    d = tubes[k]; h = float(d['H'][iu].min()); e = float(d['E'].min())
    tier.append('A' if (h >= -7 and e >= -4.5) else 'B' if (h >= -8.5 and e >= -5.5) else 'C'); mh.append(h); me.append(e); i5.append(d['i5']); i7.append(d['i7'])
tier = np.array(tier); mh = np.array(mh); me = np.array(me)
B = {'A': 0, 'C': 1, 'G': 2, 'T': 3}; M5 = np.array([[B[c] for c in s] for s in i5]); M7 = np.array([[B[c] for c in s] for s in i7])
cst = np.where(tier == 'A', 0.0, np.where(tier == 'B', 3.0, 30.0)) + 0.2 * np.maximum(0, -(mh + 7)) + 0.2 * np.maximum(0, -(me + 4.5))
def pen(S, N):
    p = 0.0
    for M in (M5, M7):
        Ms = M[S]
        for pos in range(10):
            c = np.bincount(Ms[:, pos], minlength=4); gc = c[1] + c[2]; red = c[0] + c[1]; green = c[0] + c[3]
            for v in (gc, red, green): p += max(0, 0.4 * N - v) ** 2 + max(0, v - 0.6 * N) ** 2
            for v in c: p += max(0, 0.15 * N - v) ** 2
    return p
def energy(S, N):
    conf = T[np.ix_(S, S)].sum() / 2
    return 1000 * conf + 1000 * pen(S, N) + cst[S].sum(), conf
def anneal(N, pool, sd, iters):
    rnd = random.Random(sd); pool = list(pool); S = rnd.sample(pool, N); rest = [x for x in pool if x not in S]; e = energy(S, N)[0]; best = (e, list(S))
    for it in range(iters):
        Tm = 60.0 * (1 - it / iters) + 0.05; i = rnd.randrange(N); j = rnd.randrange(len(rest)); a = S[i]; S[i] = rest[j]; e2 = energy(S, N)[0]
        if e2 <= e or rnd.random() < math.exp((e - e2) / Tm): rest[j] = a; e = e2
        else: S[i] = a
        if e < best[0]: best = (e, list(S))
    return best
res = {}
for N, nested in ((96, False), (48, True)):
    pool = range(n) if not nested else res[96]['idx']
    allb = [anneal(N, pool, seed + s, 250000 if N == 96 else 120000) for s in range(4)]
    e, S = min(allb, key=lambda x: x[0]); en, conf = energy(S, N)
    res[N] = dict(idx=S, energy=e, conflicts=int(conf), pen=float(pen(S, N)), tiers={t: int((tier[S] == t).sum()) for t in 'ABC'})
    print(N, res[N]['tiers'], '冲突', conf, '周期越界', res[N]['pen'], flush=True)
json.dump(dict(res=res, tier=tier.tolist(), mh=mh.tolist(), me=me.tolist(), i5=i5, i7=i7), open(out, 'w'))
