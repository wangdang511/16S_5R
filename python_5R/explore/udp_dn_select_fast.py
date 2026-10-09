"""快速版：从自建候选管里选 N 组（48/96），增量计算冲突和每周期配色。目标：A 级尽量多；所选序列两两编辑距离 ≥4、反向互补 ≥4；每个周期 GC、A+C、A+T 占比 40–60% 且 A/C/G/T 各 ≥15%（i5、i7 各自）。
用法：python3 udp_dn_select_fast.py <out json> <conflict npy> <seed> <iters> <N 列表，如 96,48> <pkl...>   （管名 DNxxx 的编号 = 行号）"""
import sys, pickle, json, random, math
import numpy as np
out, confl, seed, iters, Ns = sys.argv[1:6]; seed = int(seed); iters = int(iters); Ns = [int(x) for x in Ns.split(',')]; pkls = sys.argv[6:]
T = np.load(confl).astype(np.int16); n = T.shape[0]
tubes = {}
for p in pkls:
    r = pickle.load(open(p, 'rb')); N0 = len(r['EX'])
    for d in r['tubes'].values(): tubes[int(d['name'][2:])] = d
assert len(tubes) >= n, len(tubes)
iu = np.triu_indices(N0); tier = []; mh = []; me = []; i5 = []; i7 = []
for k in range(1, n + 1):
    d = tubes[k]; h = float(d['H'][iu].min()); e = float(d['E'].min())
    tier.append('A' if (h >= -7 and e >= -4.5) else 'B' if (h >= -8.5 and e >= -5.5) else 'C'); mh.append(h); me.append(e); i5.append(d['i5']); i7.append(d['i7'])
tier = np.array(tier); mh = np.array(mh); me = np.array(me)
import os
own = np.load(os.environ['OWN']) if os.environ.get('OWN') else np.zeros(n, bool)
print('分级', {t: int((tier == t).sum()) for t in 'ABC'}, '管内自冲突', int(own.sum()), flush=True)
import os
cst = 1e6 * own + np.where(tier == 'A', 0.0, np.where(tier == 'B', 3.0, float(os.environ.get('CC', 30)))) + 0.2 * np.maximum(0, -(mh + 7)) + 0.2 * np.maximum(0, -(me + 4.5))
B = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
M = np.zeros((2, n, 10), int)
for k in range(n):
    for p in range(10): M[0, k, p] = B[i5[k][p]]; M[1, k, p] = B[i7[k][p]]
LO_, HI_, BM_ = float(os.environ.get('PLO', 0.4)), float(os.environ.get('PHI', 0.6)), float(os.environ.get('PBASE', 0.15))
def pen_cnt(cnt, N):
    gc = cnt[..., 1] + cnt[..., 2]; red = cnt[..., 0] + cnt[..., 1]; green = cnt[..., 0] + cnt[..., 3]; v = 0.0
    for x in (gc, red, green): v += (np.maximum(0, LO_ * N - x) ** 2 + np.maximum(0, x - HI_ * N) ** 2).sum()
    return v + (np.maximum(0, BM_ * N - cnt) ** 2).sum()
def run(N, pool, sd):
    rnd = random.Random(sd); pool = list(pool); S = rnd.sample(pool, N); inS = np.zeros(n, bool); inS[S] = True
    cnt = np.zeros((2, 10, 4), int)
    for s in S:
        for t in range(2):
            for p in range(10): cnt[t, p, M[t, s, p]] += 1
    cc = T[:, S].sum(1).astype(int)             # 每个管与已选管的冲突数
    conf = int(cc[S].sum() / 2); pen = pen_cnt(cnt, N); csum = cst[S].sum(); E = 1000 * conf + 1000 * pen + csum; best = (E, list(S))
    poolarr = np.array(pool)
    for it in range(iters):
        if it % 300000 == 0 and it: json.dump(dict(N=N, it=it, E=best[0], idx=best[1]), open(out + f'.part{N}_{sd}', 'w'))
        Tm = 60.0 * (1 - it / iters) + 0.05; i = rnd.randrange(N); a = S[i]; b = int(poolarr[rnd.randrange(len(poolarr))])
        if inS[b]: continue
        dconf = int(cc[b] - T[b, a] - cc[a])
        for t in range(2):
            for p in range(10): cnt[t, p, M[t, a, p]] -= 1; cnt[t, p, M[t, b, p]] += 1
        npen = pen_cnt(cnt, N); dE = 1000 * dconf + 1000 * (npen - pen) + cst[b] - cst[a]
        if dE <= 0 or rnd.random() < math.exp(-dE / Tm):
            S[i] = b; inS[a] = False; inS[b] = True; cc += T[:, b] - T[:, a]; conf += dconf; pen = npen; E += dE
            if E < best[0]: best = (E, list(S))
        else:
            for t in range(2):
                for p in range(10): cnt[t, p, M[t, a, p]] += 1; cnt[t, p, M[t, b, p]] -= 1
    return best
res = {}
for N in Ns:
    import os
    pool = range(n) if not res else res[max(res)]['idx']
    if os.environ.get('POOL'): pool = [k for k in pool if tier[k] in os.environ['POOL']]
    allb = [run(N, pool, seed + s) for s in range(3)]
    e, S = min(allb, key=lambda x: x[0])
    conf = int(T[np.ix_(S, S)].sum() / 2); cnt = np.zeros((2, 10, 4), int)
    for s in S:
        for t in range(2):
            for p in range(10): cnt[t, p, M[t, s, p]] += 1
    res[N] = dict(idx=S, conflicts=conf, pen=float(pen_cnt(cnt, N)), tiers={t: int((tier[S] == t).sum()) for t in 'ABC'})
    print(N, res[N]['tiers'], '冲突', conf, '周期越界', res[N]['pen'], flush=True)
json.dump(dict(res=res, tier=tier.tolist(), mh=mh.tolist(), me=me.tolist(), i5=i5, i7=i7), open(out, 'w'))
