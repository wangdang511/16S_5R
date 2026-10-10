"""udp_designer 公共函数：配置、panel 读取、引物展开、热力学计算、等级、序列规则、每周期配色。
所有步骤脚本都 `from udp_common import *`。热力学参数全部来自 config.json，不写死。"""
import os, sys, re, json, csv, pickle, itertools, math, random
import numpy as np

DEFAULT_CFG = {
    "thermo": {"temp_c": 58.0, "mv_conc": 50, "dv_conc": 2, "dntp_conc": 0.2, "dna_conc": 250},
    "tier": {"A_H": -7.0, "A_E": -4.5, "B_H": -8.5, "B_E": -5.5,    # H=全局异源二聚体最差ΔG，E=3′端锚定最差ΔG（kcal/mol）
             "mode": "delta",          # delta：不带尾巴就超标的引物对按“加尾巴后额外恶化”计分，其余用绝对阈值（推荐）；
                                       # auto：基线超标时整体放宽阈值（旧做法，一对极端引物会让阈值失去意义）；absolute：不做任何处理
             "delta_grace": 1.0,       # delta 模式容差：基线已超标的引物对，加尾巴后再恶化不超过这么多（kcal/mol）不计；任何尾巴挂在完美双链末端都会带来约 0.7–0.8 的通用稳定
             "margin_A_H": 1.3, "margin_A_E": 0.7, "gap_B_H": 1.5, "gap_B_E": 1.0},
    "tail": {"len": 10, "gc_min": 4, "gc_max": 6, "min_edit": 4, "cross_rc_min": 3, "palin_min": 4},
    "per_cycle": {"lo": 0.4, "hi": 0.6, "base_min": 0.15},
    "universe": {"n_each": 480, "order_a": "random", "order_b": "random", "seed": 1},
    "pool": {"n_each": 144, "iters": 800000, "seed": 1},
    "pairs": {"n_pairs": 96, "n_core": 48, "iters": 800000, "seed": 1},
    "hot": {"margin_H": 2.5, "margin_E": 2.0, "lint_min_k": 5, "extra_file": None},
    "procs": os.cpu_count() or 4,
}
ADAPTERS = ["AGATCGGAAGAGC", "CTGTCTCTTATACACATCT", "AGATGTGTATAAGAGACAG", "AATGATACGGCGACCACCGAGATCTACAC",
            "CAAGCAGAAGACGGCATACGAGAT", "TCGTCGGCAGCGTC", "GTCTCGTGGGCTCGG"]       # TruSeq / Nextera / ME / P5 / P7 片段

def deep_update(a, b):
    for k, v in b.items():
        if isinstance(v, dict) and isinstance(a.get(k), dict): deep_update(a[k], v)
        else: a[k] = v
    return a
def load_cfg(path=None):
    cfg = json.loads(json.dumps(DEFAULT_CFG))
    if path: deep_update(cfg, json.load(open(path)))
    return cfg

COMP = str.maketrans("ACGTacgt", "TGCAtgca")
def revcomp(s): return s[::-1].translate(COMP)
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
         "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
def expand(seq):
    return ["".join(p) for p in itertools.product(*[IUPAC[c] for c in seq.upper()])]

# ---------------- panel ----------------
def load_panel(path):
    """panel.csv 列：name, site, orient(F/R，缺省取 site 最后一个字符), seq(5′→3′，允许 IUPAC)。返回 oligo 列表。"""
    out = []
    for r in csv.DictReader(open(path, encoding="utf-8-sig")):
        o = (r.get("orient") or r["site"][-1]).strip().upper()
        assert o in "FR", f"orient 必须是 F 或 R：{r}"
        out.append(dict(name=r["name"].strip(), site=r["site"].strip(), orient=o, seq=r["seq"].strip().upper()))
    names = [x["name"] for x in out]; assert len(set(names)) == len(names), "引物名重复"
    return out
def build_ex(panel):
    """展开后的全部引物：EX[i] = (name, orient, seq)。"""
    return [(o["name"], o["orient"], e) for o in panel for e in expand(o["seq"])]

# ---------------- 热力学 ----------------
_T = {}; _OFF = {}
def init_thermo(cfg):
    _T.update(cfg["thermo"])
    p = os.path.join(cfg.get("_work", ""), "offsets.pkl")
    if cfg.get("_work") and cfg["tier"].get("mode", "delta") == "delta" and os.path.exists(p) and not _OFF: _OFF.update(pickle.load(open(p, "rb")))
def adj(h, e, x, y):
    """delta 模式：(x, y) 是不带尾巴的两条引物；若它们不带尾巴时就已比 A 阈值差，把这部分基线超标量加回去，
    只按尾巴带来的额外恶化计分。其他引物对原样返回。"""
    o = _OFF.get((x, y))
    return (h - o[0], e - o[1]) if o else (h, e)
def _kw(): return dict(mv_conc=_T["mv_conc"], dv_conc=_T["dv_conc"], dntp_conc=_T["dntp_conc"], dna_conc=_T["dna_conc"])
def H(a, b):
    import primer3
    return primer3.calc_heterodimer(a, b, temp_c=_T["temp_c"], **_kw()).dg / 1000
def Edir(a, b):
    """a 的 3′端锚定在 b 上的最差稳定性"""
    import primer3
    return primer3.calc_end_stability(a, b, temp_c=_T["temp_c"], **_kw()).dg / 1000
def E(a, b): return min(Edir(a, b), Edir(b, a))
def HP(a):
    import primer3
    return primer3.calc_hairpin(a, temp_c=_T["temp_c"], **_kw()).dg / 1000

# ---------------- 等级 ----------------
def tier_of(h, e, th):
    return "A" if (h >= th["A_H"] and e >= th["A_E"]) else "B" if (h >= th["B_H"] and e >= th["B_E"]) else "C"
def tier_cost(h, e, th):
    t = {"A": 0.0, "B": 3.0, "C": 30.0}[tier_of(h, e, th)]
    return t + 0.2 * max(0.0, th["A_H"] - h) + 0.2 * max(0.0, th["A_E"] - e)
def load_thresholds(work, cfg):
    p = os.path.join(work, "thresholds.json")
    th = json.load(open(p)) if os.path.exists(p) else dict(cfg["tier"])
    return {k: th[k] for k in ("A_H", "A_E", "B_H", "B_E")}

# ---------------- 尾巴规则 ----------------
def kmers6():
    ks = set()
    for a in ADAPTERS:
        for t in (a, revcomp(a)): ks |= {t[i:i + 6] for i in range(len(t) - 5)}
    return ks
def single_ok(s, cfg, K6, bad=frozenset()):
    t = cfg["tail"]; g = sum(c in "GC" for c in s)
    if not t["gc_min"] <= g <= t["gc_max"]: return False
    if re.search(r"(.)\1\1", s) or re.search(r"(..)\1\1", s) or "GGGG" in s or s[:2] == "GG": return False
    if any(s[i:i + 6] in K6 for i in range(len(s) - 5)): return False
    return s not in bad and revcomp(s) not in bad

# ---------------- 每周期配色 ----------------
BASE = {"A": 0, "C": 1, "G": 2, "T": 3}
def to_mat(seqs): return np.array([[BASE[c] for c in s] for s in seqs])
def cycle_pen_counts(cnt, N, pc):
    """cnt: (L,4) 每个周期 A/C/G/T 计数；返回越界平方和（GC、A+C、A+T 占比 lo–hi，各碱基 ≥ base_min）"""
    gc = cnt[:, 1] + cnt[:, 2]; red = cnt[:, 0] + cnt[:, 1]; green = cnt[:, 0] + cnt[:, 3]; v = 0.0
    for x in (gc, red, green): v += (np.maximum(0, pc["lo"] * N - x) ** 2 + np.maximum(0, x - pc["hi"] * N) ** 2).sum()
    return v + (np.maximum(0, pc["base_min"] * N - cnt) ** 2).sum()
def cycle_ok(seqs, pc):
    N = len(seqs); M = to_mat(seqs); L = M.shape[1]; cnt = np.zeros((L, 4), int)
    for p in range(L): cnt[p] = np.bincount(M[:, p], minlength=4)
    return cycle_pen_counts(cnt, N, pc) == 0

# ---------------- 其他 ----------------
def lev_cdist(A, B, workers=2):
    from rapidfuzz.process import cdist
    from rapidfuzz.distance import Levenshtein
    return cdist(A, B, scorer=Levenshtein.distance, dtype=np.uint8, workers=workers)
def read_tsv(path): return [l.rstrip("\n").split("\t") for l in open(path, encoding="utf-8") if l.strip() and not l.startswith("#")]
def common_args(desc):
    import argparse
    ap = argparse.ArgumentParser(description=desc)
    ap.add_argument("--work", required=True, help="工作目录（所有中间结果都放这里，可断点续算）")
    ap.add_argument("--cfg", default=None, help="config.json（只写要覆盖的项）")
    ap.add_argument("--panel", default=None, help="panel.csv")
    return ap
def setup(args):
    os.makedirs(args.work, exist_ok=True); cfg = load_cfg(args.cfg); cfg["_work"] = args.work; init_thermo(cfg); return cfg
def wp(work, name): return os.path.join(work, name)
def pool_init(cfg): init_thermo(cfg)
