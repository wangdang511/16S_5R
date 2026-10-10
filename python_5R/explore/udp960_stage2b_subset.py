"""从 144+144 候选池里再选 100 个 i5 + 100 个 i7 子集（先算这个子集的第三阶段）：
代价（第二阶段的单序列代价）最低，且 i5 子集、i7 子集各自每个周期配色达标（N=100）；已经算完第三阶段的行（i5）给予奖励以便复用。
用法：python3 udp960_stage2b_subset.py <pool json> <stage1 pkl> <tails tsv> <stage3 pkl|-> <out json> <seed>"""
import sys, json, pickle, os, random, math
import numpy as np
poolj, s1, tsv, s3, out, seed = sys.argv[1:7]; seed = int(seed)
pool = json.load(open(poolj)); seqs = {l.split()[0]: l.split()[1] for l in open(tsv)}
done_rows = set(pickle.load(open(s3, 'rb'))['rows']) if s3 != '-' and os.path.exists(s3) else set()
B = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
def pen(cnt, N):
    gc = cnt[:, 1] + cnt[:, 2]; red = cnt[:, 0] + cnt[:, 1]; green = cnt[:, 0] + cnt[:, 3]; v = 0.0
    for x in (gc, red, green): v += (np.maximum(0, 0.4 * N - x) ** 2 + np.maximum(0, x - 0.6 * N) ** 2).sum()
    return v + (np.maximum(0, 0.15 * N - cnt) ** 2).sum()
def pick(names, cost, N=100, iters=200000, bonus=None):
    M = np.array([[B[c] for c in seqs[n]] for n in names]); n = len(names); rnd = random.Random(seed)
    cost = np.array(cost) - (np.array(bonus) if bonus is not None else 0)
    S = rnd.sample(range(n), N); rest = [x for x in range(n) if x not in S]
    def cnt_of(S):
        c = np.zeros((10, 4), int)
        for s in S:
            for p in range(10): c[p, M[s, p]] += 1
        return c
    cnt = cnt_of(S); E = 1000 * pen(cnt, N) + cost[S].sum(); best = (E, list(S))
    for it in range(iters):
        T = 30 * (1 - it / iters) + 0.02; i = rnd.randrange(N); j = rnd.randrange(len(rest)); a, b = S[i], rest[j]
        for p in range(10): cnt[p, M[a, p]] -= 1; cnt[p, M[b, p]] += 1
        E2 = 1000 * pen(cnt, N) + cost[S].sum() - cost[a] + cost[b]
        if E2 <= E or rnd.random() < math.exp((E - E2) / T): S[i] = b; rest[j] = a; E = E2
        else:
            for p in range(10): cnt[p, M[a, p]] += 1; cnt[p, M[b, p]] -= 1
        if E < best[0]: best = (E, list(S))
    S = best[1]; return [names[x] for x in S], pen(cnt_of(S), N)
i5, p5 = pick(pool['i5'], pool['c5'], bonus=[0.6 if k in done_rows else 0 for k in range(len(pool['i5']))])
i7, p7 = pick(pool['i7'], pool['c7'])
print('i5 子集越界', p5, '平均代价 %.2f' % np.mean([pool['c5'][pool['i5'].index(x)] for x in i5]), '，其中已算完的行', sum(pool['i5'].index(x) in done_rows for x in i5))
print('i7 子集越界', p7, '平均代价 %.2f' % np.mean([pool['c7'][pool['i7'].index(x)] for x in i7]))
json.dump(dict(rows=[pool['i5'].index(x) for x in i5], cols=[pool['i7'].index(x) for x in i7], i5=i5, i7=i7), open(out, 'w'))
