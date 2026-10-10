"""第 7 步：整管实测。对每一对 (i5,i7)，把尾巴加到 panel 的全部展开引物上（F 加 i5，R 加 i7），算全部两两
（含自身）的全局二聚体 H 和 3′端锚定 E，取最差值分级。这是最终判据，不依赖热点代理。
同时把“实测比预测差”的管里最差的引物名对写进 W/hot_extra.json（--update-hot），下次重跑第 2、3、5、6 步时作为热点。
计算量：每管 ≈ 1.5·N² 次 primer3 调用（N=展开引物数）；每个任务 = 一个管，可断点续算，--shard k/n 可多机并行。
用法：python3 s7_tubes.py --work W --panel panel.csv --pairs pairs.json --tag base [--shard 0/1] [--update-hot]
输出：W/tubes_<tag>/<管名>.pkl 与 W/tubes_<tag>.summary.json"""
import time, glob
from multiprocessing import Pool
from udp_common import *

G = {}
def _init(cfg, ex, gi, ng):
    init_thermo(cfg); G.update(ex=ex, gi=gi, ng=ng)
def eval_seqs(seqs, gi, ng):
    """seqs：带尾巴的全部展开引物；返回 minH, minE, PH(G×G), PE(G×G), minHP"""
    n = len(seqs); PH = np.zeros((ng, ng), np.float32); PE = np.zeros((ng, ng), np.float32)
    for i in range(n):
        a = seqs[i]
        for j in range(i, n):
            b = seqs[j]; h = H(a, b); e = Edir(a, b)
            if j != i: e = min(e, Edir(b, a))
            x, y = gi[i], gi[j]
            if h < PH[x, y]: PH[x, y] = PH[y, x] = h
            if e < PE[x, y]: PE[x, y] = PE[y, x] = e
    return float(PH.min()), float(PE.min()), PH, PE, min(HP(s) for s in seqs)
def tube_task(args):
    name, t5, t7 = args; ex = G["ex"]; seqs = [(t5 if o == "F" else t7) + s for _, o, s in ex]
    h, e, PH, PE, hp = eval_seqs(seqs, G["gi"], G["ng"]); return dict(name=name, i5=t5, i7=t7, minH=h, minE=e, PH=PH, PE=PE, hp=hp)

def tube_names(pairs): return [f"T{k+1:03d}" for k in range(len(pairs))]
def prep(panel):
    ex = build_ex(panel); names = [o["name"] for o in panel]; nid = {n: i for i, n in enumerate(names)}
    return ex, np.array([nid[x[0]] for x in ex]), len(names), names

if __name__ == "__main__":
    ap = common_args(__doc__); ap.add_argument("--pairs", default="pairs.json"); ap.add_argument("--tag", default="base"); ap.add_argument("--shard", default="0/1"); ap.add_argument("--update-hot", action="store_true")
    args = ap.parse_args(); cfg = setup(args); th = load_thresholds(args.work, cfg); k0, n0 = [int(x) for x in args.shard.split("/")]
    panel = load_panel(args.panel); ex, gi, ng, names = prep(panel); pairs = json.load(open(wp(args.work, args.pairs)))["pairs"]; tn = tube_names(pairs)
    d = wp(args.work, f"tubes_{args.tag}"); os.makedirs(d, exist_ok=True)
    todo = [(tn[k], p["seq5"], p["seq7"]) for k, p in enumerate(pairs) if k % n0 == k0 and not os.path.exists(os.path.join(d, tn[k] + ".pkl"))]
    print(f"{len(pairs)} 个管，本分片待算 {len(todo)}；展开引物 {len(ex)} 条，每管约 {int(1.5*len(ex)**2)} 次调用", flush=True); t0 = time.time()
    with Pool(cfg["procs"], _init, (cfg, ex, gi, ng)) as p:
        for q, r in enumerate(p.imap_unordered(tube_task, todo), 1):
            pickle.dump(r, open(os.path.join(d, r["name"] + ".pkl"), "wb")); print(f"  完成 {q}/{len(todo)} {r['name']} ({time.time()-t0:.0f}s)", flush=True)
    files = [os.path.join(d, x + ".pkl") for x in tn]
    if not all(os.path.exists(f) for f in files): print("还有管没算完（其他分片），先不汇总"); sys.exit()
    res = [pickle.load(open(f, "rb")) for f in files]; act = [tier_of(r["minH"], r["minE"], th) for r in res]; pred = [p.get("tier") for p in pairs]
    print("整管实测：", {t: act.count(t) for t in "ABC"}, "前 n_core：", {t: act[:cfg['pairs']['n_core']].count(t) for t in "ABC"})
    if all(pred): print("预测：", {t: pred.count(t) for t in "ABC"}, "一致率 %.3f" % np.mean([a == b for a, b in zip(act, pred)]), "实测比预测差的管", sum(a > b for a, b in zip(act, pred)))
    miss = {}
    for r, a, p_ in zip(res, act, pred):
        if p_ and a > p_:                        # 字母越靠后越差
            for M, key in ((r["PH"], "H"), (r["PE"], "E")):
                i, j = np.unravel_index(np.argmin(M), M.shape); miss[tuple(sorted((names[i], names[j])))] = miss.get(tuple(sorted((names[i], names[j]))), 0) + 1
    hot0 = {tuple(x) for x in json.load(open(wp(args.work, "hot_pairs.json")))["hot"]} if os.path.exists(wp(args.work, "hot_pairs.json")) else set()
    new = {k: v for k, v in miss.items() if k not in hot0 and (k[1], k[0]) not in hot0}
    json.dump(dict(tiers=dict(zip(tn, act)), counts={t: act.count(t) for t in "ABC"}, missed_pairs=[[a, b, v] for (a, b), v in sorted(new.items(), key=lambda x: -x[1])]), open(wp(args.work, f"tubes_{args.tag}.summary.json"), "w"), ensure_ascii=False, indent=1)
    print("代理漏掉的引物对（不在热点里）：", sorted(new.items(), key=lambda x: -x[1])[:8])
    if args.update_hot and new:
        f = wp(args.work, "hot_extra.json"); old = json.load(open(f)) if os.path.exists(f) else []
        json.dump(old + [list(k) for k in new if list(k) not in old], open(f, "w")); print("已写入", f, "；请重跑 s2_hotpairs、s3_single、s5_fr、s6_pair（可加 --meas）")
