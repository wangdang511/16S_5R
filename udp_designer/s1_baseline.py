"""第 1 步：panel 基线（不带尾巴）+ 设计缺陷预检。
  1) 所有展开引物两两（含自身）的全局二聚体 H0 和 3′端锚定 E0（不带尾巴）→ W/baseline.npz
  2) 自动阈值：基线最差值若已经比绝对阈值差，A/B 阈值按基线放宽 → W/thresholds.json（避免“所有管都是 C”）
  3) 接合处预检 lint：X 的 3′端 与 Y 的 5′端“跨尾巴交界处”互补（本项目里 16S-A2-R.1 × 16S-A6-F.2 就是这种，
     不带尾巴时几乎没有相互作用，加尾巴后尾巴末位补上一个配对）→ W/lint.tsv
用法：python3 s1_baseline.py --work W --panel panel.csv [--cfg config.json]"""
import time
from multiprocessing import Pool
from udp_common import *

EXV = None
def _init(cfg, ex):
    global EXV; init_thermo(cfg); EXV = ex
def _row(i):
    n = len(EXV); a = EXV[i][2]; h = np.zeros(n, np.float32); e = np.zeros(n, np.float32)
    for j in range(i, n):
        b = EXV[j][2]; h[j] = H(a, b); e[j] = Edir(a, b)
    er = np.zeros(n, np.float32)
    for j in range(i, n): er[j] = Edir(EXV[j][2], a)       # E0[j,i]
    return i, h, e, er

def lint(ex, kmin):
    """返回 {(nameX,nameY): (k, m, j, gc, X3, Y5)}：X 的 3′端 m-mer 的反向互补 w，Y 的 5′端 k=m-j 个碱基等于 w[j:]（j=1,2 个碱基由尾巴末位提供）"""
    pref = {}
    for nm, o, s in ex:
        for k in range(kmin, 8): pref.setdefault(s[:k], set()).add(nm)
    res = {}
    for nx, ox, x in ex:
        for m in range(kmin, 8):
            w = revcomp(x[-m:])
            for j in (1, 2):
                k = m - j
                if k < kmin: continue
                for ny in pref.get(w[j:], ()):
                    gc = sum(c in "GC" for c in w[j:]); key = (nx, ny)
                    if key not in res or (k, gc) > res[key][:1] + res[key][3:4]:
                        res[key] = (k, m, j, gc, x[-m:], w[j:])
    return res

if __name__ == "__main__":
    ap = common_args(__doc__); args = ap.parse_args(); cfg = setup(args)
    panel = load_panel(args.panel); ex = build_ex(panel); N = len(ex); names = [o["name"] for o in panel]; nid = {n: i for i, n in enumerate(names)}
    print(f"panel：{len(panel)} 条寡核苷酸 → {N} 条展开引物（F {sum(1 for x in ex if x[1]=='F')}，R {sum(1 for x in ex if x[1]=='R')}）", flush=True)
    out = wp(args.work, "baseline.npz")
    if not os.path.exists(out):
        H0 = np.zeros((N, N), np.float32); E0 = np.zeros((N, N), np.float32); t0 = time.time()
        with Pool(cfg["procs"], _init, (cfg, ex)) as p:
            for k, (i, h, e, er) in enumerate(p.imap_unordered(_row, range(N), chunksize=4), 1):
                H0[i, i:] = h[i:]; H0[i:, i] = h[i:]; E0[i, i:] = e[i:]; E0[i:, i] = er[i:]
                if k % 200 == 0: print(f"  基线 {k}/{N} ({time.time()-t0:.0f}s)", flush=True)
        np.savez(out, H0=H0, E0=E0)
    z = np.load(out); H0, E0 = z["H0"], z["E0"]; Es = np.minimum(E0, E0.T)
    h0, e0 = float(H0.min()), float(Es.min()); th = dict(cfg["tier"])
    if th.get("auto"):
        th["A_H"] = min(th["A_H"], round(h0 - th["margin_A_H"], 1)); th["A_E"] = min(th["A_E"], round(e0 - th["margin_A_E"], 1))
        th["B_H"] = min(th["B_H"], round(th["A_H"] - th["gap_B_H"], 1)); th["B_E"] = min(th["B_E"], round(th["A_E"] - th["gap_B_E"], 1))
    json.dump({k: th[k] for k in ("A_H", "A_E", "B_H", "B_E")}, open(wp(args.work, "thresholds.json"), "w"))
    print(f"无尾巴基线最差：全局 {h0:.2f}，3′端 {e0:.2f} kcal/mol")
    print("采用的阈值：", {k: th[k] for k in ("A_H", "A_E", "B_H", "B_E")}, "（绝对阈值" + (" 已按基线放宽）" if (th['A_H'], th['A_E']) != (cfg['tier']['A_H'], cfg['tier']['A_E']) else " 未改动）"))
    # 名称级矩阵（每对引物名的最差值），后面选热点引物对用
    G = len(names); gi = np.array([nid[x[0]] for x in ex]); PH = np.zeros((G, G), np.float32); PE = np.zeros((G, G), np.float32)
    idx = [np.where(gi == g)[0] for g in range(G)]
    for a in range(G):
        for b in range(a, G):
            PH[a, b] = PH[b, a] = H0[np.ix_(idx[a], idx[b])].min(); PE[a, b] = PE[b, a] = Es[np.ix_(idx[a], idx[b])].min()
    np.savez(wp(args.work, "baseline_names.npz"), PH=PH, PE=PE, names=np.array(names))
    L = lint(ex, cfg["hot"]["lint_min_k"])
    rows = sorted(((v[0], v[3], k[0], k[1], v) for k, v in L.items()), reverse=True)
    with open(wp(args.work, "lint.tsv"), "w") as fh:
        fh.write("# X的3′端与Y的5′端跨尾巴交界处互补：k=Y 5′端匹配碱基数（尾巴末位 j 个碱基补齐其余），gc=这 k 个碱基里 GC 数\n#X\tY\tk\tj\tgc\tX_3′_mer\tY_5′_match\n")
        for k, gc, a, b, v in rows: fh.write(f"{a}\t{b}\t{v[0]}\t{v[2]}\t{v[3]}\t{v[4]}\t{v[5]}\n")
    print(f"接合处预检：{len(rows)} 对（前 5）：", [(r[2], r[3], r[0]) for r in rows[:5]])
    worst = np.dstack(np.unravel_index(np.argsort(PH, axis=None)[:10 * 2], PH.shape))[0]
    print("基线最差的引物对（全局）：", sorted({(names[a], names[b], round(float(PH[a, b]), 2)) for a, b in worst if a <= b}, key=lambda x: x[2])[:5])
