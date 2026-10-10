"""第 5 步：候选池里每个 i5 × 每个 i7 的“正向 × 反向”热点相互作用（两个尾巴都加上）。
只算热点引物对里 正向×反向 的展开序列对；热点按基线严重程度排序，一旦已经比 B 级还差（全局 < B_H 或 3′端 < B_E）就提前终止（记为 C，early=True）。
这是整个流程里最贵的一步：计算量 ∝ M² × 热点 FR 对数 × 展开数乘积。可用 --shard k/n 在多台机器上分行并行，结果各自一个 pkl，再用 s6_pair.py 一起读。
可断点续算（每行一存）。用法：python3 s5_fr.py --work W --panel panel.csv [--cfg config.json] [--shard 0/1]
输出：W/fr_<k>of<n>.pkl"""
import time
from multiprocessing import Pool
from udp_common import *

G = {}
def _init(cfg, th, fr, T5, T7):
    init_thermo(cfg); G.update(th=th, fr=fr, T5=T5, T7=T7); G["order"] = [(f, (None, fr[f][1])) for f in fr]
def Esym(a, b): return min(Edir(a, b), Edir(b, a))
def row(i):
    th = G["th"]; t5 = G["T5"][i]; Fl = {f: [t5 + e for e in G["fr"][f][0]] for f in G["fr"]}; nj = len(G["T7"])
    hs = np.zeros(nj, np.float32); es = np.zeros(nj, np.float32); early = np.zeros(nj, bool)
    for j in range(nj):
        t7 = G["T7"][j]; h = e = 0.0
        for f, (_, rlist) in G["order"]:
            for r in rlist:
                Rl = [t7 + x for x in r[1]]
                for a in Fl[f]:
                    for b in Rl:
                        h = min(h, H(a, b)); e = min(e, Esym(a, b))
            if h < th["B_H"] or e < th["B_E"]: early[j] = True; break
        hs[j], es[j] = h, e
    return i, hs, es, early

if __name__ == "__main__":
    ap = common_args(__doc__); ap.add_argument("--shard", default="0/1"); args = ap.parse_args(); cfg = setup(args); th = load_thresholds(args.work, cfg)
    k, n = [int(x) for x in args.shard.split("/")]
    panel = load_panel(args.panel); EXP = {o["name"]: expand(o["seq"]) for o in panel}; ORI = {o["name"]: o["orient"] for o in panel}
    z = np.load(wp(args.work, "baseline_names.npz")); PH = z["PH"]; names = list(z["names"]); nid = {x: i for i, x in enumerate(names)}
    hot = [tuple(x) for x in json.load(open(wp(args.work, "hot_pairs.json")))["hot"]]
    pairs = sorted([((a, b) if ORI[a] == "F" else (b, a)) for a, b in hot if ORI[a] != ORI[b]], key=lambda p: PH[nid[p[0]], nid[p[1]]])
    # 以 F 引物名为键：[(F 展开列表, [(R 名, R 展开列表), ...])]，F 的顺序按其最严重的热点对
    byF = {}
    for f, r in pairs: byF.setdefault(f, []).append((r, EXP[r]))
    fr = {f: (EXP[f], rl) for f, rl in byF.items()}; order = [(f, (None, rl)) for f, rl in byF.items()]
    pool = json.load(open(wp(args.work, "pool.json"))); seq = {r[0]: r[1] for r in read_tsv(wp(args.work, "universe.tsv"))}
    T5 = [seq[x] for x in pool["i5"]]; T7 = [seq[x] for x in pool["i7"]]
    out = wp(args.work, f"fr_{k}of{n}.pkl"); st = pickle.load(open(out, "rb")) if os.path.exists(out) else dict(rows={})
    todo = [i for i in range(len(T5)) if i % n == k and i not in st["rows"]]
    print(f"热点 FR 引物对 {len(pairs)}；本分片 {sum(1 for i in range(len(T5)) if i % n == k)} 行，待算 {len(todo)}", flush=True); t0 = time.time()
    with Pool(cfg["procs"], _init, (cfg, th, fr, T5, T7)) as p:
        for q, (i, hs, es, early) in enumerate(p.imap_unordered(row, todo), 1):
            st["rows"][i] = (hs, es, early); pickle.dump(st, open(out, "wb"))
            if q % 5 == 0 or q == len(todo): print(f"  完成行 {len(st['rows'])} ({time.time()-t0:.0f}s)", flush=True)
    st["i5"] = pool["i5"]; st["i7"] = pool["i7"]; pickle.dump(st, open(out, "wb")); print("分片完成")
