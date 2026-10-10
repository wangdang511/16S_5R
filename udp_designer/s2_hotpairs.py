"""第 2 步：选“热点引物对”（按引物名）。只有这些对会在后面的步骤里带尾巴精确计算，计算量因此降低一到两个数量级。
热点 = (a) 基线（不带尾巴）全局 ΔG 比 A 阈值只差 margin_H 以内，或 3′端 ΔG 在 margin_E 以内；
       (b) 第 1 步 lint 标出的“跨尾巴交界处互补”引物对；
       (c) 用户补充文件（hot.extra_file，或 W/hot_extra.json）——第 7 步发现代理漏掉的引物对会自动写进去。
热点只是代理，会漏（本项目第一轮代理预测 66A，整管实测 61A）；所以第 6 步一定要整管实测，并把漏掉的对补回这里重跑。
用法：python3 s2_hotpairs.py --work W --panel panel.csv [--cfg config.json]
输出：W/hot_pairs.json"""
from udp_common import *

if __name__ == "__main__":
    ap = common_args(__doc__); args = ap.parse_args(); cfg = setup(args); th = load_thresholds(args.work, cfg); hc = cfg["hot"]
    z = np.load(wp(args.work, "baseline_names.npz")); PH, PE, names = z["PH"], z["PE"], list(z["names"]); G = len(names)
    hot = {}
    for a in range(G):
        for b in range(a, G):
            why = []
            if PH[a, b] <= th["A_H"] + hc["margin_H"]: why.append(f"基线全局 {PH[a,b]:.2f}")
            if PE[a, b] <= th["A_E"] + hc["margin_E"]: why.append(f"基线3′端 {PE[a,b]:.2f}")
            if why: hot[(names[a], names[b])] = why
    n_base = len(hot)
    for l in read_tsv(wp(args.work, "lint.tsv")):
        if len(l) >= 3 and l[0] in names and l[1] in names:
            key = tuple(sorted((l[0], l[1]), key=names.index)); hot.setdefault(key, []).append(f"lint k={l[2]}")
    ex_files = [hc.get("extra_file"), wp(args.work, "hot_extra.json")]
    for f in ex_files:
        if f and os.path.exists(f):
            for a, b in json.load(open(f)):
                key = tuple(sorted((a, b), key=names.index)); hot.setdefault(key, []).append("补充")
    items = [dict(a=k[0], b=k[1], why=v) for k, v in hot.items()]
    json.dump(dict(hot=[[i["a"], i["b"]] for i in items], detail=items, total_pairs=G * (G + 1) // 2, thresholds=th), open(wp(args.work, "hot_pairs.json"), "w"), ensure_ascii=False, indent=1)
    print(f"热点引物对 {len(items)} / {G*(G+1)//2}（基线 {n_base}，加 lint 与补充后 {len(items)}）")
