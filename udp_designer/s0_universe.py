"""第 0 步：生成候选尾巴库（universe）。
枚举全部 4^L 条序列 → 单条规则（GC、同聚物、二核苷酸重复、GGGG、接头 6-mer、非回文、排除名单）→
贪心取两个“类型”各 n_each 条（默认 480+480=960，每条之后都可以当 i5 或 i7 用）：
  类型 A 内两两 Levenshtein ≥ min_edit；类型 B 内同样；A 与 B 之间 x 对 rc(y) 的编辑距离 ≥ cross_rc_min。
（两类之间不要求 Levenshtein ≥ min_edit：全体两两 ≥4 且反向互补 ≥3 的最大独立集只有约 840 条，不够 960；
 第 3 步选候选池时再要求入选的全部序列两两 ≥ min_edit。）
用法：python3 s0_universe.py --work W [--cfg config.json] [--exclude 现有条形码.txt]
输出：W/universe.tsv（名称 序列 类型）"""
import itertools, time
import numpy as np
from udp_common import *

def greedy(cands, n, D, R, other=None, log=""):
    """cands 已按优先级排好。同类内 Levenshtein ≥ D；other（另一类已选序列）要求 x 对 rc(o) 与 rc(x) 对 o 的编辑距离 ≥ R。"""
    S = []; t0 = time.time()
    for b0 in range(0, len(cands), 4000):
        B = cands[b0:b0 + 4000]
        if S:
            B = [x for x, bad in zip(B, lev_cdist(B, S, 4).min(1) < D) if not bad]
        if other and B:
            ro = [revcomp(o) for o in other]
            bad = (lev_cdist(B, ro, 4).min(1) < R) | (lev_cdist([revcomp(x) for x in B], other, 4).min(1) < R)
            B = [x for x, bd in zip(B, bad) if not bd]
        from rapidfuzz.distance import Levenshtein
        acc = []
        for x in B:
            if all(Levenshtein.distance(x, y) >= D for y in acc): acc.append(x)
            if len(S) + len(acc) >= n: break
        S += acc
        if (b0 // 4000) % 10 == 0: print(f"  {log} 候选 {b0}/{len(cands)} 已选 {len(S)} ({time.time()-t0:.0f}s)", flush=True)
        if len(S) >= n: break
    return S[:n]

if __name__ == "__main__":
    ap = common_args(__doc__); ap.add_argument("--exclude", default=None, help="要避开的现有条形码（每行一条，含其反向互补）")
    args = ap.parse_args(); cfg = setup(args); t = cfg["tail"]; u = cfg["universe"]
    out = wp(args.work, "universe.tsv")
    if os.path.exists(out): print("已存在，跳过", out); sys.exit()
    bad = set()
    if args.exclude: bad = {l.split()[0].upper() for l in open(args.exclude) if l.strip()}
    K6 = kmers6(); L = t["len"]
    allc = [s for s in map("".join, itertools.product("ACGT", repeat=L)) if single_ok(s, cfg, K6, bad)]
    from rapidfuzz.distance import Levenshtein
    allc = [s for s in allc if Levenshtein.distance(s, revcomp(s)) >= t["palin_min"]]
    print("单条规则通过", len(allc), flush=True)
    def order(sym):
        if sym == "random": r = random.Random(u["seed"]); x = list(allc); r.shuffle(x); return x      # 随机顺序：组成更分散，但最大可取条数更少
        tb = str.maketrans(sym, "ACGT"); return sorted(allc, key=lambda x: x.translate(tb))
    A = greedy(order(u["order_a"]), u["n_each"], t["min_edit"], t["cross_rc_min"], None, "类型A")
    Aset = set(A); B = greedy([x for x in order(u["order_b"]) if x not in Aset], u["n_each"], t["min_edit"], t["cross_rc_min"], A, "类型B")
    print(f"类型A {len(A)} 条，类型B {len(B)} 条（目标 {u['n_each']}）")
    if len(A) < u["n_each"] or len(B) < u["n_each"]: print("警告：不足目标条数；可降低 universe.n_each 或放宽 tail.min_edit")
    with open(out, "w") as fh:
        for i, s in enumerate(A + B, 1): fh.write(f"U{i:03d}\t{s}\t{'A' if i <= len(A) else 'B'}\n")
    print("写出", out)
