"""第 5 步：候选池里每个 i5 × 每个 i7 的“正向 × 反向”热点相互作用（两个尾巴都加上）。
只算热点引物对里 正向×反向 的展开序列对；热点按基线严重程度排序，一旦已经比 B 级还差就提前终止（early=True，记为 C）。
**增量**：结果文件记录已算过的热点对；第 7 步发现漏掉的引物对写进 hot_extra.json、重跑第 2 步后，本步只补算新增的热点对并合并（取最差值），不会从头算。
这是整个流程里最贵的一步：计算量 ∝ M² × 热点 FR 对数 × 展开数乘积。可用 --shard k/n 多机并行（每个分片一个 pkl，第 6 步一起读）。
用法：python3 s5_fr.py --work W --panel panel.csv [--cfg config.json] [--shard 0/1]        输出：W/fr_<k>of<n>.pkl"""
import time
from multiprocessing import Pool
from udp_common import *

G = {}
def _init(cfg, th, EXP, T5, T7):
    init_thermo(cfg); G.update(th=th, EXP=EXP, T5=T5, T7=T7)
def Esym(a, b): return min(Edir(a, b), Edir(b, a))
def row(args):
    i, plist, h0, e0 = args; th, EXP = G["th"], G["EXP"]; t5 = G["T5"][i]; nj = len(G["T7"])
    hs = h0.copy() if h0 is not None else np.zeros(nj, np.float32); es = e0.copy() if e0 is not None else np.zeros(nj, np.float32); early = np.zeros(nj, bool)
    Fl = {f: [t5 + e for e in EXP[f]] for f, _ in plist}
    for j in range(nj):
        t7 = G["T7"][j]; h, e = float(hs[j]), float(es[j])
        for f, r in plist:
            Rl = [t7 + x for x in EXP[r]]
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
    pairs = sorted({((a, b) if ORI[a] == "F" else (b, a)) for a, b in hot if ORI[a] != ORI[b]}, key=lambda p: PH[nid[p[0]], nid[p[1]]])
    pool = json.load(open(wp(args.work, "pool.json"))); seq = {r[0]: r[1] for r in read_tsv(wp(args.work, "universe.tsv"))}
    T5 = [seq[x] for x in pool["i5"]]; T7 = [seq[x] for x in pool["i7"]]
    out = wp(args.work, f"fr_{k}of{n}.pkl"); st = pickle.load(open(out, "rb")) if os.path.exists(out) else dict(rows={}, done_pairs=[])
    done = {tuple(x) for x in st.get("done_pairs", [])}; new = [p for p in pairs if p not in done]
    mine = [i for i in range(len(T5)) if i % n == k]; t0 = time.time()
    # 任务：没算过的行 → 全部热点对；已算过的行 → 只补新增热点对
    todo = []
    for i in mine:
        old = st["rows"].get(i)
        if old is None: todo.append((i, pairs, None, None))
        elif new: todo.append((i, new, old[0], old[1]))
    print(f"热点 FR 引物对 {len(pairs)}（新增 {len(new)}）；本分片 {len(mine)} 行，待算 {len(todo)}", flush=True)
    if todo:
        with Pool(cfg["procs"], _init, (cfg, th, EXP, T5, T7)) as p:
            for q, (i, hs, es, early) in enumerate(p.imap_unordered(row, todo), 1):
                st["rows"][i] = (hs, es, early)
                if q % 3 == 0 or q == len(todo): pickle.dump(st, open(out, "wb")); print(f"  完成 {q}/{len(todo)} ({time.time()-t0:.0f}s)", flush=True)
    st["done_pairs"] = [list(p) for p in pairs]; st["i5"] = pool["i5"]; st["i7"] = pool["i7"]; pickle.dump(st, open(out, "wb")); print("分片完成")
