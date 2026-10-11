"""生成随机的“合成 panel”，只用于测速和估算规模（不用于设计）：n 对引物（每对一条正向 + 一条反向寡核苷酸，19–22 nt，GC 40–60%，含若干简并碱基）。
用法：python3 make_synthetic_panel.py 200 out.csv [seed] [每条平均展开数目标，默认 6]"""
import sys, random, csv
n, out = int(sys.argv[1]), sys.argv[2]; seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1; target = float(sys.argv[4]) if len(sys.argv) > 4 else 6.0
rnd = random.Random(seed); DEG = ["R", "Y", "K", "M", "S", "W"]
def one():
    while True:
        L = rnd.randint(19, 22); s = [rnd.choice("ACGT") for _ in range(L)]; gc = sum(c in "GC" for c in s) / L
        if 0.4 <= gc <= 0.6 and s[-1] in "ACGT" and s[-2:] != ["G", "G"]: break
    k = rnd.choice([1, 2, 2, 3])                       # 简并位置数（1–3 个二倍简并 → 2/4/8 个展开）
    for p in rnd.sample(range(2, L - 3), k): s[p] = rnd.choice(DEG)
    return "".join(s)
with open(out, "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["name", "site", "orient", "seq"])
    for i in range(1, n + 1):
        w.writerow([f"S{i:03d}-F", f"A{i}-F", "F", one()]); w.writerow([f"S{i:03d}-R", f"A{i}-R", "R", one()])
print("写出", out, f"{n} 对 = {2*n} 条寡核苷酸")
