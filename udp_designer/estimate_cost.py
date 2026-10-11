"""按本机实测的 primer3 速度，估算一个 panel 跑完整流程的 CPU 小时（和 N 核下的墙钟时间）。
用法：python3 estimate_cost.py panel.csv [--pairs 96] [--pool 144] [--tails 960] [--hot-frac 0.077] [--cores 16] [--bench 2000]
说明：热点比例 hot-frac 取自示例 panel（43 / 561 个引物名对 = 7.7%）；大 panel 的热点比例未知，只能先按这个估，
第 2 步跑完后用 hot_pairs.json 里的真实数量替换（--hot-pairs N）。s5 没有计入“一旦比 B 级差就提前终止”，是上界。"""
import sys, time, argparse, random
from udp_common import *
ap = argparse.ArgumentParser(); ap.add_argument("panel"); ap.add_argument("--pairs", type=int, default=96); ap.add_argument("--pool", type=int, default=144); ap.add_argument("--tails", type=int, default=960)
ap.add_argument("--hot-frac", type=float, default=0.077); ap.add_argument("--hot-pairs", type=int, default=0); ap.add_argument("--cores", type=int, default=16); ap.add_argument("--bench", type=int, default=2000); a = ap.parse_args()
init_thermo(load_cfg()); panel = load_panel(a.panel); ex = build_ex(panel); N = len(ex); G = len(panel); rnd = random.Random(1); seqs = [x[2] for x in ex]
def bench(tail):
    t0 = time.time()
    for _ in range(a.bench):
        x, y = rnd.choice(seqs), rnd.choice(seqs)
        if tail: x = "".join(rnd.choice("ACGT") for _ in range(10)) + x; y = "".join(rnd.choice("ACGT") for _ in range(10)) + y
        H(x, y); Edir(x, y); Edir(y, x)
    return 3 * a.bench / (time.time() - t0)
cps0 = bench(False); cps = bench(True)             # 每个 CPU 核每秒 primer3 调用数（H + 两个方向的 E）：不带尾巴 / 带 10 nt 尾巴（长序列慢约 3–4 倍）
avg_e = N / G; nF = sum(1 for o in panel if o["orient"] == "F"); nR = G - nF
hot_names = a.hot_pairs or int(a.hot_frac * G * (G + 1) / 2); hot_fr = hot_names * (2 * nF * nR) / (G * G)        # 热点里 F×R 对的期望数
def line(name, calls, tailed=True):
    cpuh = calls / (cps if tailed else cps0) / 3600; print(f"{name:34s} {calls:12.3e} 次调用  {cpuh:9.1f} CPU·小时  → {a.cores} 核约 {cpuh / a.cores:8.1f} 小时")
print(f"panel：{G} 条寡核苷酸 → {N} 条展开引物（平均 {avg_e:.1f} 个/条；F {nF}，R {nR}）；本机 primer3 约 {cps0:,.0f} 次调用/秒/核（不带尾巴），{cps:,.0f}（带 10 nt 尾巴）")
print(f"热点引物名对（估）{hot_names}，其中 F×R {hot_fr:.0f}\n")
line("s1 基线(不带尾巴) 1.5·N²", 1.5 * N * N, False)
line(f"s3 单条尾巴 × {a.tails}", a.tails * (3 * hot_names * avg_e * avg_e * 2.0 + 3 * N))
line(f"s5 FR {a.pool}×{a.pool}（上界）", a.pool ** 2 * hot_fr * avg_e * avg_e * 3)
line(f"s7 整管实测 × {a.pairs}", a.pairs * 1.5 * N * N)
line("s7 每轮（再测一遍）", a.pairs * 1.5 * N * N)
