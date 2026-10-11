"""第 9 步：导出最终方案：Excel（说明、每对 UDP、每周期统计公式、独立验证）+ 完整引物 csv（尾巴 + 特异部分，可直接订购）。
分级来自第 7 / 8 步的整管实测。同时导出完整 n_pairs 对和核心 n_core 对两套。
用法：python3 s9_export.py --work W --panel panel.csv --pairs pairs_repaired.json --tubes tubes_base_rep --out OUTDIR [--name 方案名] [--exclude 现有条形码.txt]"""
import collections
from udp_common import *

def verify(pairs, cfg, bad):
    from rapidfuzz.distance import Levenshtein
    a = [p["seq5"] for p in pairs]; b = [p["seq7"] for p in pairs]; S = a + b; N = len(pairs); pc = cfg["per_cycle"]
    Dm = lev_cdist(S, S).astype(int); np.fill_diagonal(Dm, 99)
    cross = min(lev_cdist(a, [revcomp(s) for s in b]).min(), lev_cdist([revcomp(s) for s in a], b).min())
    own = min(Levenshtein.distance(s, revcomp(s)) for s in S)
    return {"序列数": len(set(S)), "最小编辑距离（全部序列）": int(Dm.min()), "跨类 i5 对 rc(i7) 最小编辑距离": int(cross), "单条与自身反向互补最小": int(own),
            "每周期全部达标": bool(cycle_ok(a, pc) and cycle_ok(b, pc)), "GC 范围": (min(sum(c in "GC" for c in s) for s in S), max(sum(c in "GC" for c in s) for s in S)),
            "同聚物（≥3）条数": sum(bool(re.search(r"(.)\1\1", s)) for s in S), "与现有名单相同的条数": sum(s in bad or revcomp(s) in bad for s in S)}

def write_xlsx(path, pairs, res, cfg, name, panel, ver, prefix):
    import openpyxl
    from openpyxl import Workbook
    from openpyxl.styles import Font, PatternFill, Alignment
    from openpyxl.utils import get_column_letter
    th = cfg["_th"]; N = len(pairs); tc = collections.Counter(tier_of(r["minH"], r["minE"], th) for r in res)
    wb = Workbook(); s = wb.active; s.title = "说明"
    L = [f"{name}：{N} 对内联 UDP", f"分级（整管实测）：A {tc['A']} / B {tc['B']} / C {tc['C']}。A：全局二聚体 ≥ {th['A_H']} 且 3′端 ≥ {th['A_E']} kcal/mol；B：≥ {th['B_H']} 且 ≥ {th['B_E']}；C：其余（经验阈值，无实验标定）。",
         "独立验证：" + "；".join(f"{k}={v}" for k, v in ver.items()), f"命名：第 n 对的 i5 = {prefix}5nnn（加在所有正向引物 5′端），i7 = {prefix}7nnn（加在所有反向引物 5′端）。",
         "没有评估：真实接头；实验验证（建议同一样品、多条形码对照）。"]
    for i, t in enumerate(L, 1): c = s.cell(i, 1, t); c.alignment = Alignment(wrap_text=True, vertical="top"); c.font = Font(name="Arial", size=12 if i == 1 else 10, bold=i == 1)
    s.column_dimensions["A"].width = 150
    w = wb.create_sheet(f"{N}对"); heads = ["序号", "i5 名称", "i7 名称", "i5 序列", "i7 序列", "分级", "全局二聚体最差 ΔG", "3′端最差 ΔG", "发夹最差 ΔG"] + [f"i5 第{q+1}位" for q in range(10)] + [f"i7 第{q+1}位" for q in range(10)]
    for j, h in enumerate(heads, 1): w.cell(1, j, h).font = Font(name="Arial", bold=True)
    fills = {"A": "E2EFDA", "B": "FFF2CC", "C": "F8CBAD"}
    for k, (p, r) in enumerate(zip(pairs, res)):
        t = tier_of(r["minH"], r["minE"], th)
        vals = [k + 1, f"{prefix}5{k+1:03d}", f"{prefix}7{k+1:03d}", p["seq5"], p["seq7"], t, round(r["minH"], 2), round(r["minE"], 2), round(r.get("hp", 0), 2)]
        for j, v in enumerate(vals, 1): w.cell(k + 2, j, v).font = Font(name="Arial", size=10)
        w.cell(k + 2, 6).fill = PatternFill("solid", fgColor=fills[t])
        for q in range(10): w.cell(k + 2, 10 + q, f"=MID($D{k+2},{q+1},1)"); w.cell(k + 2, 20 + q, f"=MID($E{k+2},{q+1},1)")
    c1 = wb.create_sheet("每周期统计")
    for j, h in enumerate(["序列", "周期", "A", "C", "G", "T", "GC%", "A+C%", "A+T%", "是否达标"], 1): c1.cell(1, j, h).font = Font(name="Arial", bold=True)
    row = 2; pc = cfg["per_cycle"]
    for nm, base in (("i5", 10), ("i7", 20)):
        for p in range(10):
            hc = get_column_letter(base + p); rng = f"'{N}对'!${hc}$2:${hc}${N+1}"; c1.cell(row, 1, nm); c1.cell(row, 2, p + 1)
            for j, b in enumerate("ACGT", 3): c1.cell(row, j, f'=COUNTIF({rng},"{b}")')
            c1.cell(row, 7, f"=(E{row}+F{row})/{N}"); c1.cell(row, 8, f"=(C{row}+D{row})/{N}"); c1.cell(row, 9, f"=(C{row}+F{row})/{N}")
            c1.cell(row, 10, f'=IF(AND(G{row}>={pc["lo"]},G{row}<={pc["hi"]},H{row}>={pc["lo"]},H{row}<={pc["hi"]},I{row}>={pc["lo"]},I{row}<={pc["hi"]},MIN(C{row}:F{row})>={round(pc["base_min"]*N,2)}),"是","否")'); row += 1
    wb.save(path)

if __name__ == "__main__":
    ap = common_args(__doc__); ap.add_argument("--pairs", default="pairs_repaired.json"); ap.add_argument("--tubes", default="tubes_base_rep"); ap.add_argument("--out", required=True)
    ap.add_argument("--name", default="自建 UDP"); ap.add_argument("--exclude", default=None); ap.add_argument("--prefix", default="DN")
    args = ap.parse_args(); cfg = setup(args); cfg["_th"] = load_thresholds(args.work, cfg); os.makedirs(args.out, exist_ok=True)
    panel = load_panel(args.panel); pairs = json.load(open(wp(args.work, args.pairs)))["pairs"]; tn = tube_names = [f"T{k+1:03d}" for k in range(len(pairs))]
    res = [pickle.load(open(os.path.join(wp(args.work, args.tubes), n + ".pkl"), "rb")) for n in tn]
    bad = {l.split()[0].upper() for l in open(args.exclude) if l.strip()} if args.exclude else set()
    for N in sorted({len(pairs), cfg["pairs"]["n_core"]} - {0}):
        P, R = pairs[:N], res[:N]; ver = verify(P, cfg, bad); stem = os.path.join(args.out, f"udp_{N}pairs")
        write_xlsx(stem + ".xlsx", P, R, cfg, args.name, panel, ver, args.prefix)
        with open(stem + "_full_oligos.csv", "w", newline="", encoding="utf-8-sig") as fh:
            cw = csv.writer(fh); cw.writerow(["对", "i5 名称", "i7 名称", "分级", "引物名", "位点", "尾巴", "特异部分", "完整序列 5′→3′", "长度"])
            for k, (p, r) in enumerate(zip(P, R)):
                t5, t7 = p["seq5"], p["seq7"]; tr = tier_of(r["minH"], r["minE"], cfg["_th"])
                for o in panel:
                    t = t5 if o["orient"] == "F" else t7; tag = f"{args.prefix}{'5' if o['orient']=='F' else '7'}{k+1:03d}"
                    cw.writerow([k + 1, f"{args.prefix}5{k+1:03d}", f"{args.prefix}7{k+1:03d}", tr, f"{o['name']}-{tag}", o["site"], t, o["seq"], t + o["seq"], len(t + o["seq"])])
        print(f"{N} 对：", {t: sum(tier_of(r['minH'], r['minE'], cfg['_th']) == t for r in R) for t in 'ABC'}, ver)
