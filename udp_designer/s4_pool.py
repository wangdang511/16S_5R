"""第 4 步：从 universe 里选候选池：M 个 i5 + M 个 i7。
硬约束（以很大的惩罚实现，模拟退火）：入选的 2M 条两两 Levenshtein ≥ tail.min_edit；i5 与 i7 之间 x 对 rc(y) 与 rc(x) 对 y 的编辑距离 ≥ tail.cross_rc_min；
          i5 池、i7 池各自每个周期 GC、A+C、A+T 占比在 [lo,hi]，A/C/G/T 各 ≥ base_min。
代价：第 3 步单条尾巴的分级（i5 角色看 FF、FRu，i7 角色看 RR、RFu，另加发夹与“尾巴×引物”代理）。
M 选多大：第 5 步要算 M×M 个 i5×i7 组合，计算量 ∝ M²。留冗余（M ≥ 1.5×要选的对数）让配对有余地，但 M 太大第 5 步会很慢。
用法：python3 s4_pool.py --work W [--cfg config.json]     输出：W/pool.json"""
from udp_common import *

if __name__ == "__main__":
    ap = common_args(__doc__); args = ap.parse_args(); cfg = setup(args); th = load_thresholds(args.work, cfg); t = cfg["tail"]; pc = cfg["per_cycle"]
    M = cfg["pool"]["n_each"]; iters = cfg["pool"]["iters"]; seed = cfg["pool"]["seed"]; out = wp(args.work, "pool.json")
    if os.path.exists(out): print("已存在，跳过", out); sys.exit()
    rows = read_tsv(wp(args.work, "universe.tsv")); names = [r[0] for r in rows]; seqs = [r[1] for r in rows]; n = len(seqs); assert 2 * M <= n, "universe 不够大"
    S1 = pickle.load(open(wp(args.work, "single.pkl"), "rb"))
    D = lev_cdist(seqs, seqs); Dr = lev_cdist(seqs, [revcomp(s) for s in seqs])
    CD = (D < t["min_edit"]).astype(np.int16); np.fill_diagonal(CD, 0)
    CR = ((Dr < t["cross_rc_min"]) | (Dr.T < t["cross_rc_min"])).astype(np.int16); np.fill_diagonal(CR, 0)
    print("冲突密度：编辑距离 %.4f，跨类反向互补 %.4f" % (CD.mean(), CR.mean()), flush=True)
    c5 = np.zeros(n); c7 = np.zeros(n)
    for i, nm in enumerate(names):
        r = S1[nm]; px = tier_cost(r["ph"], r["pe"], th)
        c5[i] = tier_cost(min(r["FF"][0], r["FRu"][0]), min(r["FF"][1], r["FRu"][1]), th) + 0.5 * max(0, -(r["hpF"] + 2)) + px
        c7[i] = tier_cost(min(r["RR"][0], r["RFu"][0]), min(r["RR"][1], r["RFu"][1]), th) + 0.5 * max(0, -(r["hpR"] + 2)) + px
    MM = to_mat(seqs); L = MM.shape[1]
    def counts(sel):
        c = np.zeros((L, 4), int)
        for s in sel:
            for p in range(L): c[p, MM[s, p]] += 1
        return c
    rnd = random.Random(seed); order = list(range(n)); rnd.shuffle(order)
    s5 = order[:M]; s7 = order[M:2 * M]; used = np.zeros(n, bool); used[s5 + s7] = True
    cnt5, cnt7 = counts(s5), counts(s7)
    cD = CD[:, s5 + s7].sum(1).astype(int); cR5 = CR[:, s5].sum(1).astype(int); cR7 = CR[:, s7].sum(1).astype(int)
    conf = int(sum(cD[s] for s in s5 + s7) / 2 + sum(cR7[s] for s in s5))
    pen5 = cycle_pen_counts(cnt5, M, pc); pen7 = cycle_pen_counts(cnt7, M, pc); cost = c5[s5].sum() + c7[s7].sum()
    E_ = 1000 * conf + 1000 * (pen5 + pen7) + cost; best = (E_, list(s5), list(s7)); print("起始 冲突", conf, "周期越界", pen5 + pen7, flush=True)
    for it in range(iters):
        T = 80.0 * (1 - it / iters) + 0.05; role = rnd.randrange(2); b = rnd.randrange(n)
        if used[b]: continue
        sl = s5 if role == 0 else s7; i = rnd.randrange(M); a = sl[i]
        cO = cR7 if role == 0 else cR5; cnt = cnt5 if role == 0 else cnt7; cc = c5 if role == 0 else c7
        dconf = (cD[b] - CD[b, a] + cO[b]) - (cD[a] + cO[a])
        for p in range(L): cnt[p, MM[a, p]] -= 1; cnt[p, MM[b, p]] += 1
        newpen = cycle_pen_counts(cnt, M, pc); oldpen = pen5 if role == 0 else pen7
        dE = 1000 * dconf + 1000 * (newpen - oldpen) + cc[b] - cc[a]
        if dE <= 0 or rnd.random() < math.exp(-dE / T):
            sl[i] = b; used[a] = False; used[b] = True; cD += CD[:, b] - CD[:, a]
            if role == 0: cR5 += CR[:, b] - CR[:, a]; pen5 = newpen
            else: cR7 += CR[:, b] - CR[:, a]; pen7 = newpen
            conf += dconf; E_ += dE
            if E_ < best[0]: best = (E_, list(s5), list(s7))
        else:
            for p in range(L): cnt[p, MM[a, p]] += 1; cnt[p, MM[b, p]] -= 1
        if it and it % 300000 == 0: print(f"  迭代 {it}/{iters} 能量 {E_:.1f} 最好 {best[0]:.1f}", flush=True)
    _, s5, s7 = best; sel = s5 + s7
    cdd = int(CD[np.ix_(sel, sel)].sum() / 2); crr = int(CR[np.ix_(s5, s7)].sum()); pen = cycle_pen_counts(counts(s5), M, pc) + cycle_pen_counts(counts(s7), M, pc)
    print(f"最终：编辑距离冲突 {cdd}，跨类反向互补冲突 {crr}，周期越界 {pen}，i5 平均代价 {c5[s5].mean():.2f}，i7 平均代价 {c7[s7].mean():.2f}")
    if cdd or crr or pen: print("警告：硬约束没有全部满足；增大 universe、增加迭代数或降低 M")
    json.dump(dict(i5=[names[x] for x in s5], i7=[names[x] for x in s7], seq={names[x]: seqs[x] for x in sel}, c5=[float(c5[x]) for x in s5], c7=[float(c7[x]) for x in s7]), open(out, "w"))
