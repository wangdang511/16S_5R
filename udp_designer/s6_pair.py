"""第 6 步：模拟退火，从 M×M 候选池里选 n_pairs 对并配对（前 n_core 对 = 核心子集，权重 ×2，自己也满足每周期配色）。
一个 i5×i7 组合的预测 = min(FR 精确值, i5 的 FF, i7 的 RR)；等级代价 A 0 / B 3 / C 30；每周期越界重罚。
修复循环用的两个开关：
  --meas meas.json   [[i5名, i7名, "A"|"B"|"C"], ...] 已整管实测过的组合，用实测等级覆盖代理预测
  --ban5 U001,U002 / --ban7 ...   禁用在几乎所有搭配里都是 C 的序列
用法：python3 s6_pair.py --work W [--cfg config.json] [--meas meas.json] [--ban7 U864] [--out pairs.json]
输出：W/pairs.json（预测等级通常比整管实测乐观，必须再跑第 7 步）"""
import glob
from udp_common import *

if __name__ == "__main__":
    ap = common_args(__doc__); ap.add_argument("--meas"); ap.add_argument("--ban5", default=""); ap.add_argument("--ban7", default=""); ap.add_argument("--out", default="pairs.json")
    args = ap.parse_args(); cfg = setup(args); th = load_thresholds(args.work, cfg); pc = cfg["per_cycle"]; pr = cfg["pairs"]; NP, NC = pr["n_pairs"], pr["n_core"]
    pool = json.load(open(wp(args.work, "pool.json"))); seq = {r[0]: r[1] for r in read_tsv(wp(args.work, "universe.tsv"))}; S1 = pickle.load(open(wp(args.work, "single.pkl"), "rb"))
    n5, n7 = pool["i5"], pool["i7"]; nr, nc = len(n5), len(n7); HM = np.zeros((nr, nc)); EM = np.zeros((nr, nc)); have = np.zeros(nr, bool)
    for f in sorted(glob.glob(wp(args.work, "fr_*of*.pkl"))):
        for i, (hs, es, early) in pickle.load(open(f, "rb"))["rows"].items(): HM[i] = hs; EM[i] = es; have[i] = True
    assert have.all(), f"FR 还有 {int((~have).sum())} 行没算完（见 s5_fr.py 的各分片）"
    FF = np.array([S1[x]["FF"] for x in n5]); RR = np.array([S1[x]["RR"] for x in n7])
    TH = np.minimum(HM, np.minimum(FF[:, 0][:, None], RR[:, 0][None, :])); TE = np.minimum(EM, np.minimum(FF[:, 1][:, None], RR[:, 1][None, :]))
    TIER = np.vectorize(lambda h, e: tier_of(h, e, th))(TH, TE); COST = np.vectorize(lambda h, e: tier_cost(h, e, th))(TH, TE)
    if args.meas:
        i5x = {x: i for i, x in enumerate(n5)}; i7x = {x: i for i, x in enumerate(n7)}
        for m in json.load(open(args.meas)):
            a, b = m[0], m[1]
            if a in i5x and b in i7x:      # 有整管实测值的按当前阈值重算（阈值可能已从绝对改为相对）
                if len(m) >= 5: TIER[i5x[a], i7x[b]] = tier_of(m[3], m[4], th); COST[i5x[a], i7x[b]] = tier_cost(m[3], m[4], th)
                else: TIER[i5x[a], i7x[b]] = m[2]; COST[i5x[a], i7x[b]] = {"A": 0.0, "B": 3.0, "C": 30.0}[m[2]]
    for x in filter(None, args.ban5.split(",")): COST[n5.index(x), :] += 1e4; TIER[n5.index(x), :] = "C"
    for x in filter(None, args.ban7.split(",")): COST[:, n7.index(x)] += 1e4; TIER[:, n7.index(x)] = "C"
    print(f"{nr}×{nc} 组合预测：", {t: int((TIER == t).sum()) for t in "ABC"}, flush=True)
    M5 = to_mat([seq[x] for x in n5]); M7 = to_mat([seq[x] for x in n7]); L = M5.shape[1]
    def pen(M, idx, N):
        c = np.zeros((L, 4), int)
        for s in idx:
            for p in range(L): c[p, M[s, p]] += 1
        return cycle_pen_counts(c, N, pc)
    rnd = random.Random(pr["seed"]); p5 = list(range(nr)); rnd.shuffle(p5); p7 = list(range(nc)); rnd.shuffle(p7)
    s5 = p5[:NP]; u5 = p5[NP:]; s7 = p7[:NP]; u7 = p7[NP:]; W = np.array([2.0] * NC + [1.0] * (NP - NC))
    def energy(a, b):
        e = (COST[a, b] * W).sum() + 1000 * (pen(M5, a, NP) + pen(M7, b, NP))
        if NC: e += 1000 * (pen(M5, a[:NC], NC) + pen(M7, b[:NC], NC))
        return e
    Ecur = energy(s5, s7); best = (Ecur, list(s5), list(s7)); iters = pr["iters"]
    for it in range(iters):
        T = 40.0 * (1 - it / iters) + 0.05; mv = rnd.random(); a5, a7, b5, b7 = list(s5), list(s7), list(u5), list(u7)
        if mv < 0.3: i, j = rnd.sample(range(NP), 2); a7[i], a7[j] = a7[j], a7[i]
        elif mv < 0.5: i, j = rnd.sample(range(NP), 2); a5[i], a5[j] = a5[j], a5[i]
        elif mv < 0.75: i = rnd.randrange(NP); k = rnd.randrange(len(b5)); a5[i], b5[k] = b5[k], a5[i]
        else: i = rnd.randrange(NP); k = rnd.randrange(len(b7)); a7[i], b7[k] = b7[k], a7[i]
        E2 = energy(a5, a7)
        if E2 <= Ecur or rnd.random() < math.exp((Ecur - E2) / T):
            s5, s7, u5, u7, Ecur = a5, a7, b5, b7, E2
            if Ecur < best[0]: best = (Ecur, list(s5), list(s7))
        if it and it % 100000 == 0: print(f"  迭代 {it}/{iters} 能量 {Ecur:.1f} 最好 {best[0]:.1f}", flush=True)
    _, s5, s7 = best; tiers = [TIER[a, b] for a, b in zip(s5, s7)]
    print("最终（预测）：", NP, "对", {t: tiers.count(t) for t in "ABC"}, "核心", {t: tiers[:NC].count(t) for t in "ABC"}, "周期越界", pen(M5, s5, NP) + pen(M7, s7, NP))
    json.dump(dict(pairs=[dict(i5=n5[a], i7=n7[b], tier=str(TIER[a, b]), H=float(TH[a, b]), E=float(TE[a, b]), seq5=seq[n5[a]], seq7=seq[n7[b]]) for a, b in zip(s5, s7)]), open(wp(args.work, args.out), "w"), ensure_ascii=False)
