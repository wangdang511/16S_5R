"""第 3 步：单条尾巴评估（universe 里每条尾巴，分别按 i5 角色加在正向引物、按 i7 角色加在反向引物）。
对热点引物对计算（每项 = 全局异源二聚体最差 ΔG H，3′端锚定最差 ΔG E）：
  FF  两条正向引物都带这条尾巴        RR  两条反向引物都带这条尾巴
  FRu 带尾巴的正向 × 不带尾巴的反向    RFu 带尾巴的反向 × 不带尾巴的正向
  hpF / hpR  带尾巴的正向 / 反向引物的发夹最差 ΔG
另外对 panel 里全部展开引物做“尾巴本身 × 引物”的代理（不依赖热点）：
  ph  尾巴与每条不带尾巴引物的全局最差 H；pe  每条不带尾巴引物的 3′端锚在尾巴上的最差 E
  （抓“尾巴自己与某条引物 3′端互补”这类热点漏掉的情形）
可断点续算。用法：python3 s3_single.py --work W --panel panel.csv [--cfg config.json]
输出：W/single.pkl"""
import time
from multiprocessing import Pool
from udp_common import *
import hashlib

G = {}
def _init(cfg, panel, hot, ex):                # hot 在此仅占位
    init_thermo(cfg); G["EXP"] = {o["name"]: expand(o["seq"]) for o in panel}; G["ORI"] = {o["name"]: o["orient"] for o in panel}
    G["HOT"] = hot; G["BODY"] = sorted({x[2] for x in ex}); G["F"] = [e for o in panel if o["orient"] == "F" for e in G["EXP"][o["name"]]]
    G["R"] = [e for o in panel if o["orient"] == "R" for e in G["EXP"][o["name"]]]
def Esym(a, b): return min(Edir(a, b), Edir(b, a))
def pm(A, B, same):
    """A、B：[(带尾巴序列, 不带尾巴序列), ...]"""
    h = e = 0.0
    for i, (a, xa) in enumerate(A):
        for j, (b, xb) in enumerate(B):
            if same and j < i: continue
            hh, ee = adj(H(a, b), Esym(a, b), xa, xb); h = min(h, hh); e = min(e, ee)
    return h, e
def do(arg):
    tail, old, plist = arg; name, t = tail[0], tail[1]; EXP, ORI = G["EXP"], G["ORI"]
    res = old if old is not None else dict(name=name, seq=t, FF=(0.0, 0.0), RR=(0.0, 0.0), FRu=(0.0, 0.0), RFu=(0.0, 0.0))
    for a, b in plist:
        oa, ob = ORI[a], ORI[b]
        if oa == ob:
            cl = "FF" if oa == "F" else "RR"; h, e = pm([(t + x, x) for x in EXP[a]], [(t + x, x) for x in EXP[b]], a == b)
            res[cl] = (min(res[cl][0], h), min(res[cl][1], e))
        else:
            f, r = (a, b) if oa == "F" else (b, a)
            h, e = pm([(t + x, x) for x in EXP[f]], [(x, x) for x in EXP[r]], False); res["FRu"] = (min(res["FRu"][0], h), min(res["FRu"][1], e))
            h, e = pm([(t + x, x) for x in EXP[r]], [(x, x) for x in EXP[f]], False); res["RFu"] = (min(res["RFu"][0], h), min(res["RFu"][1], e))
    if old is None:      # 发夹和“尾巴×引物”代理与热点无关，只算一次
        res["hpF"] = min([HP(t + e) for e in G["F"]] or [0.0]); res["hpR"] = min([HP(t + e) for e in G["R"]] or [0.0])
        res["ph"] = min(H(t, b) for b in G["BODY"]); res["pe"] = min(Edir(b, t) for b in G["BODY"])
    return res

if __name__ == "__main__":
    ap = common_args(__doc__); args = ap.parse_args(); cfg = setup(args)
    panel = load_panel(args.panel); ex = build_ex(panel); hot = [tuple(x) for x in json.load(open(wp(args.work, "hot_pairs.json")))["hot"]]
    tails = read_tsv(wp(args.work, "universe.tsv")); out = wp(args.work, "single.pkl")
    done = pickle.load(open(out, "rb")) if os.path.exists(out) else {}
    hp = wp(args.work, "single_hot_done.json"); hdone = frozenset(tuple(x) for x in json.load(open(hp))) if os.path.exists(hp) else frozenset()
    # 每条尾巴记录自己已覆盖的热点集合（hk → 热点列表存在 single_hot_sets.json），中途掉线重跑时只补缺的部分
    sp = wp(args.work, "single_hot_sets.json"); sets = {k: frozenset(tuple(x) for x in v) for k, v in json.load(open(sp)).items()} if os.path.exists(sp) else {}
    hotset = frozenset(hot); hk = hashlib.md5(repr(sorted(hotset)).encode()).hexdigest()[:12]
    if hk not in sets:
        sets[hk] = hotset; json.dump({k: [list(x) for x in sorted(v)] for k, v in sets.items()}, open(sp, "w"))
    todo = []
    for t in tails:
        if t[0] not in done: todo.append((t, None, hot)); continue
        cov = sets.get(done[t[0]].get("hk"), hdone); miss = [h for h in hot if h not in cov]
        if miss: todo.append((t, done[t[0]], miss))
    print(f"待算 {len(todo)} / {len(tails)} 条尾巴；热点引物对 {len(hot)}（已有记录 {sum(1 for t in tails if t[0] in done)} 条）", flush=True); t0 = time.time()
    with Pool(cfg["procs"], _init, (cfg, panel, hot, ex)) as p:
        for k, r in enumerate(p.imap_unordered(do, todo, chunksize=2), 1):
            r["hk"] = hk; done[r["name"]] = r
            if k % 40 == 0: pickle.dump(done, open(out + ".tmp", "wb")); os.replace(out + ".tmp", out); print(f"  完成 {k}/{len(todo)} ({time.time()-t0:.0f}s)", flush=True)
    pickle.dump(done, open(out + ".tmp", "wb")); os.replace(out + ".tmp", out); json.dump([list(h) for h in hot], open(hp, "w")); print("全部完成", len(done))
