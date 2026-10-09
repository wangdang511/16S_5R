"""阶段 2：从独立集（两两 Levenshtein ≥4、与反向互补 ≥3）挑 768 条，分成 16 组（8 块 × {i5, i7}），每组 48 条，每组每个周期 GC、A+C、A+T 占比 40–60% 且 A/C/G/T 各 ≥15%。
二聚体预筛指标作为软代价。用法：python3 udp_denovo384_stage2.py <订购清单 xlsx> <indep.pkl> <输出 json> <seed> [iters]"""
import sys, pickle, json, random, math
import numpy as np, openpyxl, primer3
from multiprocessing import Pool
import primer_design as pdz
xlsx, indep, out, seed = sys.argv[1:5]; seed = int(seed); ITERS = int(sys.argv[5]) if len(sys.argv) > 5 else 600000
KW = dict(mv_conc=50, dv_conc=2, dntp_conc=0.2, dna_conc=250); T = 58.0
ws = openpyxl.load_workbook(xlsx)['订购清单']
BODY = sorted({e for r in ws.iter_rows(min_row=4, max_row=36, values_only=True) for e in pdz.expand(r[4])})
def proxy(s):
    H = min(primer3.calc_heterodimer(s, b, temp_c=T, **KW).dg for b in BODY) / 1000
    E = min(primer3.calc_end_stability(b, s, temp_c=T, **KW).dg for b in BODY) / 1000
    return H, E
if __name__ == '__main__':
    S = pickle.load(open(indep, 'rb')); n = len(S)
    with Pool(4) as p: HE = p.map(proxy, S, chunksize=10)
    cost = np.array([max(0, -h - 3.0) + max(0, -e - 2.5) for h, e in HE]) * 5
    B = {'A': 0, 'C': 1, 'G': 2, 'T': 3}; M = np.array([[B[c] for c in s] for s in S])
    G, NG = 16, 48
    def pen(cnt):
        gc = cnt[:, 1] + cnt[:, 2]; red = cnt[:, 0] + cnt[:, 1]; green = cnt[:, 0] + cnt[:, 3]; v = 0.0
        for x in (gc, red, green): v += (np.maximum(0, 0.4 * NG - x) ** 2 + np.maximum(0, x - 0.6 * NG) ** 2).sum()
        v += (np.maximum(0, 0.15 * NG - cnt) ** 2).sum()
        return v
    rnd = random.Random(seed); order = list(range(n)); rnd.shuffle(order)
    slot = order[:G * NG]; unused = order[G * NG:]
    cnt = np.zeros((G, 10, 4), int)
    for k, s in enumerate(slot):
        g = k // NG
        for pos in range(10): cnt[g, pos, M[s, pos]] += 1
    gp = np.array([pen(cnt[g]) for g in range(G)]); tot = gp.sum() * 100 + cost[slot].sum()
    best = (tot, list(slot), list(unused)); print('起始', gp.sum(), flush=True)
    for it in range(ITERS):
        Tm = 15.0 * (1 - it / ITERS) + 0.02
        i = rnd.randrange(G * NG); gi = i // NG
        if rnd.random() < 0.5 and unused:                       # 与未用序列交换
            j = rnd.randrange(len(unused)); a, b = slot[i], unused[j]
            for pos in range(10): cnt[gi, pos, M[a, pos]] -= 1; cnt[gi, pos, M[b, pos]] += 1
            ng = pen(cnt[gi]); dE = (ng - gp[gi]) * 100 + cost[b] - cost[a]
            if dE <= 0 or rnd.random() < math.exp(-dE / Tm): slot[i] = b; unused[j] = a; gp[gi] = ng; tot += dE
            else:
                for pos in range(10): cnt[gi, pos, M[a, pos]] += 1; cnt[gi, pos, M[b, pos]] -= 1
        else:                                                   # 与另一组的槽位交换
            k = rnd.randrange(G * NG); gk = k // NG
            if gk == gi: continue
            a, b = slot[i], slot[k]
            for pos in range(10):
                cnt[gi, pos, M[a, pos]] -= 1; cnt[gi, pos, M[b, pos]] += 1; cnt[gk, pos, M[b, pos]] -= 1; cnt[gk, pos, M[a, pos]] += 1
            n1, n2 = pen(cnt[gi]), pen(cnt[gk]); dE = (n1 + n2 - gp[gi] - gp[gk]) * 100
            if dE <= 0 or rnd.random() < math.exp(-dE / Tm): slot[i], slot[k] = b, a; gp[gi], gp[gk] = n1, n2; tot += dE
            else:
                for pos in range(10):
                    cnt[gi, pos, M[a, pos]] += 1; cnt[gi, pos, M[b, pos]] -= 1; cnt[gk, pos, M[b, pos]] += 1; cnt[gk, pos, M[a, pos]] -= 1
        if tot < best[0]: best = (tot, list(slot), list(unused))
        if it % 50000 == 0: print(it, round(gp.sum(), 1), round(tot, 1), flush=True)
    tot, slot, unused = best
    cnt = np.zeros((G, 10, 4), int)
    for k, s in enumerate(slot):
        for pos in range(10): cnt[k // NG, pos, M[s, pos]] += 1
    gp = np.array([pen(cnt[g]) for g in range(G)]); print('最终：周期越界合计', gp.sum(), '软代价', cost[slot].sum())
    json.dump(dict(slots=[S[k] for k in slot], H=[HE[k][0] for k in slot], E=[HE[k][1] for k in slot], pen=gp.tolist()), open(out, 'w'))
