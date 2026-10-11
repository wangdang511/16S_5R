"""第 3b 步（相对分级）：按“初步测试”的数值分布定等级阈值，而不是用固定的经验阈值。
做法：从 universe 随机抽 K 个 (i5, i7) 组合（固定种子，不经任何优选），按第 7 步的方法整管实测（全部引物两两，不用热点代理），
得到这个 panel 上“随便加一对尾巴”时整管最差全局 ΔG（H）和 3′端 ΔG（E）的分布；
A = 分布里最好的 rel_A（默认 20%），B = 最好的 rel_B（默认 70%），其余 C。
两个维度联合取分位：找 s，使 “H ≥ q_H(s) 且 E ≥ q_E(s)” 的组合恰好占 rel_A，阈值就是 (q_H(s), q_E(s))。

为什么必须整管实测、不能用第 5/6 步的组合预测：热点代理只看一部分引物对，对大 panel 系统性偏乐观
（panel2_1 第 1 轮：代理预测 96 管 C 60，整管实测 64/64 都比代理差）；用代理分布定阈值会把阈值定得过严。
为什么放在第 4 步之前：第 4 步（候选池）、第 5 步（提前终止）、第 6 步（配对代价）、第 7/8 步（分级）都读阈值。

panel 越大（靶标越多），一个管里至少有一对差的概率越大，整管分布整体左移；相对分级让 A/B/C 始终表示“在这个 panel 上相对好/中/差”，
但**不表示绝对安全**：thresholds.json 同时保留绝对阈值，导出时两种等级都给。
已有 thresholds.json 且 grading=relative、参数相同 → 跳过（保证各轮阈值一致）。
用法：python3 s3b_calibrate.py --work W --panel panel.csv [--cfg config.json]     输出：W/thresholds.json、W/tubes_calib/、W/calib.tsv"""
import random, time
from multiprocessing import Pool
from udp_common import *
import s7_tubes as S7

def joint_q(H, E, frac):
    lo, hi = 0.0, 1.0
    for _ in range(50):
        s = (lo + hi) / 2; f = float(((H >= np.quantile(H, s)) & (E >= np.quantile(E, s))).mean())
        lo, hi = (s, hi) if f > frac else (lo, s)
    return float(np.quantile(H, s)), float(np.quantile(E, s))

if __name__ == "__main__":
    ap = common_args(__doc__); args = ap.parse_args(); cfg = setup(args); tc = cfg["tier"]
    if tc.get("grading", "relative") != "relative": print("grading=absolute，跳过"); sys.exit()
    K, pA, pB, seed = tc.get("calib_n", 48), tc.get("rel_A", 0.2), tc.get("rel_B", 0.7), tc.get("calib_seed", 7)
    tp = wp(args.work, "thresholds.json"); old = json.load(open(tp)) if os.path.exists(tp) else {}
    if old.get("grading") == "relative" and old.get("calib", {}).get("params") == [K, pA, pB, seed]: print("已校准，跳过：", {k: old[k] for k in ("A_H", "A_E", "B_H", "B_E")}); sys.exit()
    rows = read_tsv(wp(args.work, "universe.tsv")); rnd = random.Random(seed); combos = []
    while len(combos) < K:
        a, b = rnd.sample(rows, 2); combos.append((f"C{len(combos)+1:03d}", a[1], b[1]))
    panel = load_panel(args.panel); ex, gi, ng, names = S7.prep(panel); d = wp(args.work, "tubes_calib"); os.makedirs(d, exist_ok=True)
    def done(nm, a, b):
        f = os.path.join(d, nm + ".pkl")
        try: r = pickle.load(open(f, "rb")); return r["i5"] == a and r["i7"] == b
        except Exception: return False
    todo = [c for c in combos if not done(*c)]
    print(f"校准：随机 {K} 个组合整管实测，待算 {len(todo)}（每管约 {int(1.5*len(ex)**2)} 次调用）", flush=True); t0 = time.time()
    with Pool(cfg["procs"], S7._init, (cfg, ex, gi, ng)) as p:
        for q, r in enumerate(p.imap_unordered(S7.tube_task, todo), 1):
            pickle.dump(r, open(os.path.join(d, r["name"] + ".pkl"), "wb")); print(f"  完成 {q}/{len(todo)} ({time.time()-t0:.0f}s)", flush=True)
    R = [pickle.load(open(os.path.join(d, c[0] + ".pkl"), "rb")) for c in combos]; H = np.array([r["minH"] for r in R]); E = np.array([r["minE"] for r in R])
    aH, aE = joint_q(H, E, pA); bH, bE = joint_q(H, E, pB)
    ab = {k: tc[k] for k in ("A_H", "A_E", "B_H", "B_E")}
    th = dict(A_H=round(aH, 3), A_E=round(aE, 3), B_H=round(bH, 3), B_E=round(bE, 3), grading="relative", abs=ab,
              calib=dict(params=[K, pA, pB, seed], H_q=dict(zip(["p10", "p25", "p50", "p75", "p90"], np.quantile(H, [.1, .25, .5, .75, .9]).round(2).tolist())),
                         E_q=dict(zip(["p10", "p25", "p50", "p75", "p90"], np.quantile(E, [.1, .25, .5, .75, .9]).round(2).tolist()))))
    json.dump(th, open(tp, "w"), ensure_ascii=False, indent=1)
    # 校准管里每管最差的引物名对若不在热点里，补进 hot_extra.json（之后的 s2/s3 增量补算）
    hp = wp(args.work, "hot_pairs.json"); hot0 = {tuple(x) for x in json.load(open(hp))["hot"]} if os.path.exists(hp) else set(); miss = {}
    for r in R:
        for M in (r["PH"], r["PE"]):
            i, j = np.unravel_index(np.argmin(M), M.shape); k = tuple(sorted((names[i], names[j])))
            if k not in hot0 and (k[1], k[0]) not in hot0: miss[k] = miss.get(k, 0) + 1
    if miss:
        f = wp(args.work, "hot_extra.json"); ex0 = json.load(open(f)) if os.path.exists(f) else []
        json.dump(ex0 + [list(k) for k in miss if list(k) not in ex0], open(f, "w")); print("校准管里热点漏掉的引物对：", sorted(miss.items(), key=lambda x: -x[1])[:8], "→ hot_extra.json")
    with open(wp(args.work, "calib.tsv"), "w") as fh:
        fh.write("#组合\ti5\ti7\t整管最差全局ΔG\t整管最差3′端ΔG\t相对等级\t绝对等级\n")
        for c, h, e in zip(combos, H, E): fh.write(f"{c[0]}\t{c[1]}\t{c[2]}\t{h:.2f}\t{e:.2f}\t{tier_of(h, e, th)}\t{tier_of(h, e, ab)}\n")
    print(f"整管分布（随机 {K} 组合）：H {th['calib']['H_q']}  E {th['calib']['E_q']}")
    print(f"相对阈值：A（前 {pA:.0%}）H ≥ {aH:.2f} 且 E ≥ {aE:.2f}；B（前 {pB:.0%}）H ≥ {bH:.2f} 且 E ≥ {bE:.2f}；绝对阈值 {ab} → thresholds.json")
