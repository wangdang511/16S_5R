"""第 10 步（可选）：对某一条引物做最小改动并重新评估整管等级。用于“上了 UDP 后发现某条引物与别的引物在尾巴下互补，能不能微调一下”。
变体：3′端修剪 1–2 nt、5′端修剪 1–3 nt、3′端最后 --last 位的单碱基替换（长度不变）。
两种模式：
  筛选：--partners 引物名,引物名   只算变体与这些伙伴引物（及自身）的相互作用，几秒到几分钟，用来淘汰明显没用的变体
  整管：不给 --partners             变体换进每个管，对 panel 全部展开引物重算，给出整管 A/B/C 计数（慢，只对筛选后的少数变体做）
本脚本只算热力学。覆盖率（目标物种匹配率）和基因组脱靶必须另外算：
  覆盖率：用你自己的参考序列库，口径为“≤1 错配且 3′端 3 nt 无错配”；脱靶：改 3′种子区(最后 8 nt)的变体必须对基因组重新全长扫描，不能用旧命中列表（会低估）。
用法：python3 s10_tweak.py --work W --panel panel.csv --pairs pairs_repaired.json --tubes tubes_base_rep --oligo 16S-A2-R.1 [--partners 16S-A6-F.2] [--last 7] [--max-tubes 96]"""
import time
from multiprocessing import Pool
from udp_common import *

G = {}
def _init(cfg, ex, names, gi, th):
    init_thermo(cfg); G.update(ex=ex, names=names, gi=gi, th=th, delta=cfg["tier"].get("mode", "delta") == "delta", grace=cfg["tier"].get("delta_grace", 0.5))
def task(args):
    """args: (tube 序号, i5, i7, 变体序列列表(展开后), 目标引物名, 伙伴名集合或 None) → (minH, minE) 变体 × (其他引物 + 变体自身)"""
    k, t5, t7, vexp, target, partners = args; ex = G["ex"]; ori = [o for n, o, s in ex if n == target][0]; tv = t5 if ori == "F" else t7
    V = [(tv + s, s) for s in vexp]; h = e = 0.0; th = G["th"]; delta = G["delta"]
    others = [((t5 if o == "F" else t7) + s, s) for n, o, s in ex if n != target and (partners is None or n in partners)]
    for a, xa in V:
        for b, xb in others + V:
            hh = H(a, b); ee = min(Edir(a, b), Edir(b, a))
            if delta:      # 变体不在基线里：现算不带尾巴的值，扣除基线超标量
                h0 = H(xa, xb); e0 = min(Edir(xa, xb), Edir(xb, xa))
                g = G["grace"]; hh -= min(0.0, h0 - th["A_H"] - g); ee -= min(0.0, e0 - th["A_E"] - g)     # 与 s1 相同的连续规则
            h = min(h, hh); e = min(e, ee)
    return k, h, e

def variants(seq, last, orient):
    v = {"原始": seq}
    for k in (1, 2): v[f"3′修剪{k}"] = seq[:-k]
    for k in (1, 2, 3): v[f"5′修剪{k}"] = seq[k:]
    for pos in range(len(seq) - last, len(seq)):
        alts = [b for b in "ACGT" if b != seq[pos]] if seq[pos] in "ACGT" else list("ACGT")
        for b in alts: v[f"位{pos+1}:{seq[pos]}→{b}"] = seq[:pos] + b + seq[pos + 1:]
    return v

if __name__ == "__main__":
    ap = common_args(__doc__); ap.add_argument("--pairs", default="pairs_repaired.json"); ap.add_argument("--tubes", default="tubes_base_rep"); ap.add_argument("--oligo", required=True)
    ap.add_argument("--partners", default=None); ap.add_argument("--last", type=int, default=7); ap.add_argument("--max-tubes", type=int, default=0); ap.add_argument("--only", default=None, help="只评估这些变体（逗号分隔）")
    args = ap.parse_args(); cfg = setup(args); th = load_thresholds(args.work, cfg)
    panel = load_panel(args.panel); ex = build_ex(panel); names = [o["name"] for o in panel]; nid = {n: i for i, n in enumerate(names)}; gi = np.array([nid[x[0]] for x in ex])
    pairs = json.load(open(wp(args.work, args.pairs)))["pairs"]; NP = len(pairs) if not args.max_tubes else min(args.max_tubes, len(pairs)); tn = [f"T{k+1:03d}" for k in range(NP)]
    tubes = [pickle.load(open(os.path.join(wp(args.work, args.tubes), n + ".pkl"), "rb")) for n in tn]
    me = nid[args.oligo]; keep = [i for i in range(len(names)) if i != me]
    base_rest = [(float(t["PH"][np.ix_(keep, keep)].min()), float(t["PE"][np.ix_(keep, keep)].min())) for t in tubes]        # 去掉目标引物后每个管剩下的最差值
    seq0 = [o for o in panel if o["name"] == args.oligo][0]["seq"]; V = variants(seq0, args.last, None)
    if args.only: V = {k: v for k, v in V.items() if k in args.only.split(",") or k == "原始"}
    partners = set(args.partners.split(",")) if args.partners else None
    print(f"{args.oligo}：{len(V)} 个变体，{NP} 个管，" + ("筛选模式，伙伴 " + str(sorted(partners)) if partners else "整管模式"), flush=True); out = []; t0 = time.time()
    with Pool(cfg["procs"], _init, (cfg, ex, names, gi, th)) as pool:
        for vn, vs in V.items():
            vexp = expand(vs); jobs = [(k, pairs[k]["seq5"], pairs[k]["seq7"], vexp, args.oligo, partners) for k in range(NP)]
            r = sorted(pool.map(task, jobs)); h = np.array([x[1] for x in r]); e = np.array([x[2] for x in r])
            if partners:
                row = dict(variant=vn, seq=vs, 展开数=len(vexp), 管数_3端低于A=int((e < th["A_E"]).sum()), 最差3端=round(float(e.min()), 2), 管数_全局低于A=int((h < th["A_H"]).sum()), 最差全局=round(float(h.min()), 2))
            else:
                tiers = [tier_of(min(bh, hh), min(be, ee), th) for (bh, be), hh, ee in zip(base_rest, h, e)]
                row = dict(variant=vn, seq=vs, 展开数=len(vexp), A=tiers.count("A"), B=tiers.count("B"), C=tiers.count("C"))
            out.append(row); print(" ", row, f"({time.time()-t0:.0f}s)", flush=True)
    json.dump(out, open(wp(args.work, f"tweak_{args.oligo}_{'screen' if partners else 'full'}.json"), "w"), ensure_ascii=False, indent=1)
