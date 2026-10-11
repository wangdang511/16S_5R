"""第 8 步：按整管实测修复。对实测等级差于 --target 的管，和同组（核心 / 非核心）里的其他管互换 i7，
把互换后的两个管整管实测，挑总代价下降的交换。互换不改变 i5 / i7 各自的序列集合，所以编辑距离、反向互补、每周期配色都不受影响。
如果所有交换都修不好，说明问题在这个 i5 或 i7 本身：用 s6_pair.py 的 --ban5 / --ban7 禁用它们后重选。
用法：python3 s8_repair.py --work W --panel panel.csv --pairs pairs.json --tag base [--target B|C] [--trials 8]
输出：W/pairs_repaired.json 与 W/tubes_<tag>_rep/（更新后的管结果，之后可用 s9_export.py）"""
import shutil, time
from multiprocessing import Pool
from udp_common import *
from s7_tubes import _init, tube_task, prep, tube_names

if __name__ == "__main__":
    ap = common_args(__doc__); ap.add_argument("--pairs", default="pairs.json"); ap.add_argument("--tag", default="base"); ap.add_argument("--target", default="B"); ap.add_argument("--trials", type=int, default=8)
    args = ap.parse_args(); cfg = setup(args); th = load_thresholds(args.work, cfg); NC = cfg["pairs"]["n_core"]
    panel = load_panel(args.panel); ex, gi, ng, names = prep(panel); pairs = json.load(open(wp(args.work, args.pairs)))["pairs"]; tn = tube_names(pairs); NP = len(pairs)
    src = wp(args.work, f"tubes_{args.tag}"); dst = wp(args.work, f"tubes_{args.tag}_rep"); shutil.rmtree(dst, ignore_errors=True); shutil.copytree(src, dst)
    res = {n: pickle.load(open(os.path.join(dst, n + ".pkl"), "rb")) for n in tn}
    tier = lambda r: tier_of(r["minH"], r["minE"], th); cost = lambda r: tier_cost(r["minH"], r["minE"], th)
    bad = [k for k in range(NP) if tier(res[tn[k]]) >= args.target]          # "A"<"B"<"C"
    print("实测：", {t: sum(tier(res[n]) == t for n in tn) for t in "ABC"}, f"；需要修复（≥{args.target}）的管 {len(bad)}", flush=True)
    rnd = random.Random(cfg["pairs"]["seed"]); rep = 0
    with Pool(cfg["procs"], _init, (cfg, ex, gi, ng)) as pool:
        for k in bad:
            if tier(res[tn[k]]) < args.target: continue          # 前面的交换可能已经顺带修好
            grp = range(0, NC) if k < NC else range(NC, NP)
            ys = rnd.sample([y for y in grp if y != k], min(args.trials, len(grp) - 1)) if NC else rnd.sample([y for y in range(NP) if y != k], min(args.trials, NP - 1))
            jobs = []
            for y in ys: jobs += [(f"X{y}a", pairs[k]["seq5"], pairs[y]["seq7"]), (f"X{y}b", pairs[y]["seq5"], pairs[k]["seq7"])]
            out = {r["name"]: r for r in pool.map(tube_task, jobs)}
            old = cost(res[tn[k]]) 
            best = None
            for y in ys:
                a, b = out[f"X{y}a"], out[f"X{y}b"]; c_new = cost(a) + cost(b); c_old = old + cost(res[tn[y]])
                if tier(a) != "C" and tier(b) != "C" and c_new < c_old - 1e-6 and (best is None or c_new < best[0]): best = (c_new, y, a, b)
            if best is None: print(f"  管 {tn[k]}（{pairs[k]['i5']}/{pairs[k]['i7']}）{tier(res[tn[k]])}：没有可行的交换 → 考虑 --ban5 {pairs[k]['i5']} 或 --ban7 {pairs[k]['i7']}"); continue
            _, y, a, b = best; print(f"  管 {tn[k]} 与 {tn[y]} 互换 i7：{tier(res[tn[k]])}+{tier(res[tn[y]])} → {tier(a)}+{tier(b)}")
            pairs[k]["i7"], pairs[y]["i7"] = pairs[y]["i7"], pairs[k]["i7"]; pairs[k]["seq7"], pairs[y]["seq7"] = pairs[y]["seq7"], pairs[k]["seq7"]
            a["name"], b["name"] = tn[k], tn[y]; res[tn[k]], res[tn[y]] = a, b
            pickle.dump(a, open(os.path.join(dst, tn[k] + ".pkl"), "wb")); pickle.dump(b, open(os.path.join(dst, tn[y] + ".pkl"), "wb")); rep += 1
    for k, p in enumerate(pairs):
        r = res[tn[k]]; p["tier"] = tier(r); p["H"] = r["minH"]; p["E"] = r["minE"]
    json.dump(dict(pairs=pairs), open(wp(args.work, "pairs_repaired.json"), "w"), ensure_ascii=False)
    print("修复后：", {t: sum(p['tier'] == t for p in pairs) for t in "ABC"}, "；核心", {t: sum(p['tier'] == t for p in pairs[:NC]) for t in "ABC"}, f"；交换 {rep} 次")
