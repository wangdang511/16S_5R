"""
design_primers.py -- 从参考序列集端到端重现 5R 引物重设计分析

    python design_primers.py GG_DIR OUT_DIR [--skip-scan]

GG_DIR  = qiime_default_reference/gg_13_8_otus（pip download qiime-default-reference 后解压）
OUT_DIR = 结果目录（sites_key.csv, frontier.csv, designs_final.json）
步骤: 投影到 E. coli 坐标 -> 保守性剖面 -> 扫描全部候选位点 -> 精确求解 1/2 管平铺 -> 调整 Tm -> 管内二聚体检查
"""
import argparse, itertools, json, os, pickle, warnings
import numpy as np, pandas as pd
import primer_design as pdz

warnings.filterwarnings("ignore")
V_WEIGHTS = {v: 1.0 for v in pdz.V_REGIONS}
V_WEIGHTS["V9"] = 0.6          # V9 常被截短，权重略低
AMP = (180, 250)               # 含引物的扩增子长度
V4_MAX = 292                   # V4 全长 107 nt + 两端引物，放宽到 292
MIN_GAP = 20                   # 同管内相邻扩增子的最小间隔
TARGET_TM = 62.0
PLANS = {"A_two_pool": dict(kmin=0.75, kmean=0.85, n=2),
         "B_single_pool": dict(kmin=0.75, kmean=0.85, n=1),
         "C_two_pool_full": dict(kmin=0.50, kmean=0.80, n=2)}


def main(gg_dir, out_dir, skip_scan=False):
    os.makedirs(out_dir, exist_ok=True)
    R = pdz.load_reference(gg_dir, cache=os.path.join(out_dir, "reference.pkl"))
    pdz.conservation_profile(R, "Bacteria").to_csv(os.path.join(out_dir, "profile_bact.csv"))
    pdz.conservation_profile(R, "Archaea").to_csv(os.path.join(out_dir, "profile_arch.csv"))
    sites_csv = os.path.join(out_dir, "sites_key.csv")
    if skip_scan and os.path.exists(sites_csv):
        d = pd.read_csv(sites_csv)
    else:
        d = pdz.scan_sites(R, lengths=(18, 19, 20), max_degen=3, max_fold=8)
        d.to_csv(sites_csv, index=False)

    rows = []
    for kmin in [0.8, 0.75, 0.7, 0.6, 0.5, 0.4, 0.3]:
        s = d[d.key_min >= kmin]
        for n in (1, 2):
            r = pdz.optimize_pools(s[s.orient == "F"], s[s.orient == "R"], n, AMP, MIN_GAP, (), V_WEIGHTS,
                                   special_range={"V4": (AMP[0], V4_MAX)})
            cov = [v for v, c in r["vc"].items() if c > .999]
            rows.append({"key min >=": kmin, "pools": n, "n V": len(cov), "complete V": " ".join(cov),
                         "amplicons": " ".join(f"{a['lo']}-{a['hi']}{'AB'[p]}" for a, p in r["sol"])})
    pd.DataFrame(rows).to_csv(os.path.join(out_dir, "frontier.csv"), index=False)

    out = {}
    for name, p in PLANS.items():
        s = d[(d.key_min >= p["kmin"]) & (d.key_mean >= p["kmean"])]
        res = pdz.optimize_pools(s[s.orient == "F"], s[s.orient == "R"], p["n"], AMP, MIN_GAP, ("V4",), V_WEIGHTS,
                                 special_range={"V4": (AMP[0], V4_MAX)})
        prim, amps = [], []
        for k, (a, pl) in enumerate(res["sol"], 1):
            tuned = {}
            for o in "FR":
                r = a[o]._asdict()
                t = pdz.tune_tm(R, r, TARGET_TM) or dict(start=r["start"], end=r["end"], primer=r["primer"])
                ev = pdz.evaluate_primer(R, t["primer"], o, t["start"])
                b = ev[ev.domain == "Bacteria"]
                per = b.groupby("phylum").hit.mean().reindex(pdz.KEY_PHYLA)
                lo, hi = pdz.tm_range(t["primer"])
                prim.append(dict(amplicon=f"A{k}", pool="AB"[pl], name=f"A{k}-{o}", orient=o,
                                 site=f"{t['start']}-{t['end']}", primer=t["primer"], length=len(t["primer"]),
                                 fold=len(pdz.expand(t["primer"])), GC=round(pdz.gc_content(t["primer"]), 2),
                                 Tm=f"{lo:.1f}-{hi:.1f}", bacteria=round(b.hit.mean(), 3),
                                 key_mean=round(per.mean(), 3), key_min=round(per.min(), 3), worst=per.idxmin(),
                                 archaea=round(ev[ev.domain == "Archaea"].hit.mean(), 3),
                                 **{kk: round(v, 2) for kk, v in per.items()}))
                tuned[o] = t
            amps.append(dict(amplicon=f"A{k}", pool="AB"[pl], F_site=f"{tuned['F']['start']}-{tuned['F']['end']}",
                             R_site=f"{tuned['R']['start']}-{tuned['R']['end']}", start=tuned["F"]["start"],
                             end=tuned["R"]["end"], length=tuned["R"]["end"] - tuned["F"]["start"] + 1,
                             regions="+".join(v for v, c in a["vc"].items() if c > .999)))
        P = pd.DataFrame(prim)
        dims = []
        for pl in sorted(set(P.pool)):
            sub = P[P.pool == pl]
            for (n1, s1), (n2, s2) in itertools.combinations_with_replacement(zip(sub.name, sub.primer), 2):
                k = max(pdz.three_prime_dimer(s1, s2, 5), pdz.three_prime_dimer(s2, s1, 5))
                if k >= 5:
                    dims.append(dict(pool=pl, pair=f"{n1} × {n2}", nt=k))
        out[name] = dict(amplicons=amps, primers=P.to_dict("records"), dimers=dims,
                         vc={k: float(v) for k, v in res["vc"].items()})
    json.dump(out, open(os.path.join(out_dir, "designs_final.json"), "w"), indent=1, ensure_ascii=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("gg_dir"); ap.add_argument("out_dir"); ap.add_argument("--skip-scan", action="store_true")
    a = ap.parse_args()
    main(a.gg_dir, a.out_dir, a.skip_scan)
