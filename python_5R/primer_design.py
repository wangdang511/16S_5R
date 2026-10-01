"""
primer_design.py -- 16S 多重扩增引物的保守性分析与平铺（tiling）方案优化
=======================================================================

数据：Greengenes 13_8 的 85% OTU 代表序列 PyNAST 比对（4,797 条，4,536 细菌 + 261 古菌），
来自 PyPI 包 ``qiime-default-reference``。85% 聚类使每个“序列”大致代表一个科/属级谱系，
统计保守性时不会被大量重复的近缘序列（如肠杆菌）主导。

坐标：*E. coli* 16S 编号。比对中以 Salmonella（GG 4384058，与 E. coli 锚点距离仅差 0–2 nt）为参照，
再用 11 个通用引物锚点做分段线性换算，误差约 ±1 nt。

主要函数
  load_reference()            载入比对并建立 比对列 <-> E. coli 坐标
  conservation_profile()      每个 E. coli 位置的保守度
  scan_sites()                对每个候选引物位点（正向 / 反向）求“允许 ≤k 个简并碱基”时的最佳引物与覆盖率
  locate_primer()             把任意引物（含简并碱基）定位到 E. coli 坐标
  optimize_tiling()           在长度、覆盖率、不重叠约束下选扩增子组合（单管 / 两管）
"""
from __future__ import annotations

import gzip
import itertools
import os
import re
from dataclasses import dataclass

import numpy as np
import pandas as pd

# E. coli 可变区（Chakravorty et al. 2007, J Microbiol Methods）
V_REGIONS = {"V1": (69, 99), "V2": (137, 242), "V3": (433, 497), "V4": (576, 682), "V5": (822, 879),
             "V6": (986, 1043), "V7": (1117, 1173), "V8": (1243, 1294), "V9": (1435, 1465)}

# 锚点: (E. coli 起点, 正链序列)；在参照序列上找到后做分段线性换算
ANCHORS = [(8, "AGAGTTTGATC[ACM]TGGCTCAG"), (104, "GGCGGACGGGTGAGTAA"), (338, "ACTCCTACGGGAGGCAG"),
           (515, "GTGCCAGCAGCCGCGGTAA"), (685, "GTGTAGCGGTGAAATGCG"), (787, "ATTAGATACCCTGGTAGTCC"),
           (907, "AAACTCAAATGAATTGACGG"), (967, "CAACGCGAAGAACCTTACC"), (1100, "CAACGAGCGCAACCC"),
           (1175, "GGAGGAAGGTGGGGATGAC"), (1392, "GTACACACCGCCCGT"), (1492, "AAGTCGTAACAAGGTAACC")]

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "K": "GT", "M": "AC", "S": "CG",
         "W": "AT", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
CODE = {frozenset(v): k for k, v in IUPAC.items()}
COMP = str.maketrans("ACGTRYKMSWBDHVN", "TGCAYRMKSWVHDBN")


def revcomp(s):
    return s.translate(COMP)[::-1]


def _read_fasta_gz(path):
    out, h = {}, None
    with gzip.open(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                h = line[1:].split()[0]
                out[h] = []
            elif h:
                out[h].append(line)
    return {k: "".join(v) for k, v in out.items()}


@dataclass
class Reference:
    aln: np.ndarray            # (N, 1542) uint8：每条序列在 E. coli 位置 1..1542 上的碱基；ord('-') = 缺失，0 = 序列未覆盖（末端截短）
    ids: list
    domain: np.ndarray         # 'Bacteria' / 'Archaea'
    phylum: np.ndarray
    col_of_pos: dict           # E. coli 位置 -> 矩阵列（恒等映射 p -> p-1）
    ref_seq: str               # E. coli 坐标下的参照序列（Salmonella 碱基，长度 1542）
    taxonomy: dict


def _ecoli_reference(seqs97, ref_id):
    """把参照序列（Salmonella）按锚点分段线性换算到 E. coli 坐标，返回长度 1542 的字符串"""
    ref = seqs97[ref_id]
    pts = []
    for ecoli_pos, pat in ANCHORS:
        m = re.search(pat, ref)
        if m:
            pts.append((m.start(), ecoli_pos - 1))
    pts = np.array(pts)
    idx = np.arange(len(ref))
    e0 = np.interp(idx, pts[:, 0], pts[:, 1])
    e0[idx < pts[0, 0]] = pts[0, 1] - (pts[0, 0] - idx[idx < pts[0, 0]])
    e0[idx > pts[-1, 0]] = pts[-1, 1] + (idx[idx > pts[-1, 0]] - pts[-1, 0])
    pos_base = {}
    for i, e in enumerate(np.round(e0).astype(int)):
        if 0 <= e < 1542:
            pos_base.setdefault(e, ref[i])
    return "".join(pos_base.get(p, "N") for p in range(1542))


def _project_one(args):
    ref, seq = args
    from Bio import Align
    al = Align.PairwiseAligner(mode="global", match_score=2, mismatch_score=-3,
                               open_gap_score=-6, extend_gap_score=-2)
    al.end_insertion_score = 0
    al.end_deletion_score = 0
    a = al.align(ref, seq)[0]
    row = np.zeros(len(ref), np.uint8)
    blocks_t, blocks_q = a.aligned
    if len(blocks_t) == 0:
        return row
    first, last = blocks_t[0][0], blocks_t[-1][1]
    row[first:last] = ord("-")                       # 覆盖范围内默认“缺失”
    for (ts, te), (qs, qe) in zip(blocks_t, blocks_q):
        row[ts:te] = np.frombuffer(seq[qs:qe].encode(), np.uint8)
    return row


def load_reference(gg_dir, ref_id="4384058", n_jobs=4, cache=None):
    """
    gg_dir = .../qiime_default_reference/gg_13_8_otus
    把 85% OTU 代表序列逐条与 E. coli 坐标参照做两两全局比对（末端缺口不罚分），
    得到 (N, 1542) 的位点矩阵。对引物位点这类保守区，两两比对足够可靠。
    """
    if cache and os.path.exists(cache):
        import pickle
        return pickle.load(open(cache, "rb"))
    aligned = _read_fasta_gz(os.path.join(gg_dir, "rep_set_aligned", "85_otus.pynast.fasta.gz"))
    seqs97 = _read_fasta_gz(os.path.join(gg_dir, "rep_set", "97_otus.fasta.gz"))
    tax = {}
    with gzip.open(os.path.join(gg_dir, "taxonomy", "97_otu_taxonomy.txt.gz"), "rt") as fh:
        for line in fh:
            k, v = line.rstrip("\n").split("\t")
            tax[k] = [x.strip()[3:] for x in v.split(";")]
    ids = list(aligned)
    ref = _ecoli_reference(seqs97, ref_id)
    raw = [re.sub("[^ACGTN]", "N", aligned[i].replace("-", "").replace(".", "").upper()) for i in ids]
    from multiprocessing import Pool
    with Pool(n_jobs) as pool:
        rows = pool.map(_project_one, [(ref, s) for s in raw], chunksize=50)
    A = np.vstack(rows)
    domain = np.array([tax[i][0] for i in ids])
    phylum = np.array([tax[i][1] or "unclassified" for i in ids])
    R = Reference(A, ids, domain, phylum, {p: p - 1 for p in range(1, 1543)}, ref, tax)
    if cache:
        import pickle
        pickle.dump(R, open(cache, "wb"))
    return R


def site_matrix(R: Reference, start, length):
    """E. coli [start, start+length) 对应的比对列 -> (N, length) 碱基矩阵；位置缺失时返回 None"""
    cols = [R.col_of_pos.get(p) for p in range(start, start + length)]
    if any(c is None for c in cols):
        return None
    M = R.aln[:, cols]
    # 两列之间若有其他序列的插入（参照中为缺口），这里忽略（引物位点内插入极少）
    return M


def conservation_profile(R: Reference, domain="Bacteria"):
    rows = []
    sel = R.domain == domain
    for p in range(1, 1543):
        c = R.col_of_pos.get(p)
        if c is None:
            continue
        col = R.aln[sel, c]
        col = col[col != 0]
        if len(col) == 0:
            continue
        acgt = col[np.isin(col, np.frombuffer(b"ACGT", np.uint8))]
        if len(acgt):
            vals, cnt = np.unique(acgt, return_counts=True)
            top = cnt.max() / len(col)
            f = cnt / cnt.sum()
            ent = -(f * np.log2(f)).sum()
        else:
            top, ent = 0, 2
        rows.append({"pos": p, "top_freq": top, "entropy": ent, "gap": (col == ord("-")).mean(), "n": len(col)})
    return pd.DataFrame(rows).set_index("pos")


# --------------------------------------------------------------------------
# 引物位点评估
# --------------------------------------------------------------------------
def _match(M, primer_top, three_prime_side, max_mm=1, clamp=3):
    """
    M: (N, L) 正链碱基；primer_top: 引物在正链方向的序列（反向引物传其反向互补），可含简并碱基
    匹配定义（常用的 in-silico PCR 标准）：总错配 ≤ max_mm，且 3′ 端 clamp 个碱基内无错配
    three_prime_side: 'right'（正向引物）或 'left'（反向引物）
    """
    ok_base = np.zeros(M.shape, bool)
    for j, ch in enumerate(primer_top):
        ok_base[:, j] = np.isin(M[:, j], np.frombuffer(IUPAC[ch].encode(), np.uint8))
    mm = (~ok_base).sum(1)
    cl = ~ok_base[:, -clamp:] if three_prime_side == "right" else ~ok_base[:, :clamp]
    return (mm <= max_mm) & (cl.sum(1) == 0)


def best_degenerate(M, three_prime_side, max_degen=2, max_fold=4, weights=None):
    """
    贪心设计：先取每列多数碱基，再依次在“能挽救最多序列”的位置加入简并碱基，
    总简并度 ≤ max_fold，简并位置数 ≤ max_degen。返回 (引物正链序列, 覆盖率)
    """
    L = M.shape[1]
    data = (M != 0).all(1) & (M != ord("-")).all(1)
    Md = M[data]
    w = np.ones(len(Md)) if weights is None else weights[data]
    if len(Md) < 20:
        return None, np.nan
    sets = []
    for j in range(L):
        v, c = np.unique(Md[:, j], return_counts=True)
        v = [chr(x) for x in v]
        sets.append({v[int(np.argmax(c))]} if set(v) & set("ACGT") else {"A"})

    def seq_of(st):
        return "".join(CODE[frozenset(s)] for s in st)

    def cov(st):
        return (w * _match(Md, seq_of(st), three_prime_side)).sum() / w.sum()

    cur = cov(sets)
    for _ in range(max_degen):
        best = (cur, None)
        fold = np.prod([len(s) for s in sets])
        for j in range(L):
            for b in "ACGT":
                if b in sets[j] or fold * (len(sets[j]) + 1) / len(sets[j]) > max_fold:
                    continue
                trial = [set(s) for s in sets]
                trial[j].add(b)
                c = cov(trial)
                if c > best[0] + 0.005:
                    best = (c, trial)
        if best[1] is None:
            break
        cur, sets = best
    return seq_of(sets), cur


KEY_PHYLA = ["Firmicutes", "Proteobacteria", "Bacteroidetes", "Actinobacteria", "Fusobacteria", "Verrucomicrobia",
             "Tenericutes", "Spirochaetes", "Cyanobacteria", "Chloroflexi", "Acidobacteria", "Planctomycetes"]


def _weights(R: Reference, domain="Bacteria", key_phyla=KEY_PHYLA):
    """设计权重：主要门各占相同总权重；其他细菌门合计再占 1 份（避免完全忽略）"""
    sel = R.domain == domain
    ph = R.phylum[sel]
    w = np.zeros(sel.sum())
    groups = list(key_phyla) + ["__other__"]
    lab = np.where(np.isin(ph, key_phyla), ph, "__other__")
    for g in groups:
        m = lab == g
        if m.any():
            w[m] = 1.0 / m.sum()
    return sel, ph, w


def _scan_worker(args):
    R, positions, lengths, max_degen, max_fold, domain, key_phyla = args
    sel, ph, w = _weights(R, domain, key_phyla)
    rows = []
    for start in positions:
        for L in lengths:
            M = site_matrix(R, start, L)
            if M is None:
                continue
            M = M[sel]
            for orient, side in (("F", "right"), ("R", "left")):
                top, _ = best_degenerate(M, side, max_degen, max_fold, weights=w)
                if top is None:
                    continue
                data = (M != 0).all(1) & (M != ord("-")).all(1)
                hit = _match(M[data], top, side)
                phd = ph[data]
                per = {k: hit[phd == k].mean() for k in key_phyla if (phd == k).sum() >= 10}
                # 全部细菌门等权
                u, inv, cnt = np.unique(phd, return_inverse=True, return_counts=True)
                wb = 1.0 / cnt[inv]
                primer = top if orient == "F" else revcomp(top)
                rows.append({"start": start, "end": start + L - 1, "len": L, "orient": orient, "primer": primer,
                             "fold": int(np.prod([len(IUPAC[ch]) for ch in primer])),
                             "coverage": hit.mean(),
                             "coverage_balanced": (wb * hit).sum() / wb.sum(),
                             "key_mean": float(np.mean(list(per.values()))),
                             "key_min": float(min(per.values())),
                             "key_worst": min(per, key=per.get),
                             **{f"cov_{k}": v for k, v in per.items()}})
    return rows


def scan_sites(R: Reference, lengths=(18, 19, 20), max_degen=3, max_fold=8, domain="Bacteria",
               positions=None, key_phyla=KEY_PHYLA, n_jobs=4):
    """
    对每个起点和长度，求正向/反向两种用法的最佳简并引物（贪心加入简并碱基，目标为主要门等权覆盖率），并报告：
      coverage           所有细菌序列的命中比例
      coverage_balanced  所有细菌门等权
      key_mean / key_min 主要门（KEY_PHYLA）的平均 / 最低覆盖率；key_worst 为最差的门
    """
    positions = list(positions if positions is not None else range(1, 1524))
    chunks = [positions[i::n_jobs] for i in range(n_jobs)]
    from multiprocessing import Pool
    with Pool(n_jobs) as pool:
        parts = pool.map(_scan_worker, [(R, c, lengths, max_degen, max_fold, domain, key_phyla) for c in chunks])
    df = pd.DataFrame([r for p in parts for r in p])
    return df.sort_values(["start", "len", "orient"]).reset_index(drop=True)


def evaluate_primer(R: Reference, primer, orient, start, by="phylum", top=None):
    """给定引物（5'->3'）与 E. coli 起点，计算细菌/古菌/各门的覆盖率"""
    top_seq = primer if orient == "F" else revcomp(primer)
    side = "right" if orient == "F" else "left"
    M = site_matrix(R, start, len(top_seq))
    data = (M != 0).all(1) & (M != ord("-")).all(1)
    hit = np.zeros(len(M), bool)
    hit[data] = _match(M[data], top_seq, side)
    df = pd.DataFrame({"domain": R.domain, "phylum": R.phylum, "hit": hit, "data": data})
    df = df[df.data]
    return df


def locate_primer(ref_seq_ecoli_fn, primer, orient, max_mm=3):
    """在 E. coli 坐标的参照序列上定位引物（含简并），返回 (start, end, mismatches)"""
    top = primer if orient == "F" else revcomp(primer)
    best = (None, None, 99)
    s = ref_seq_ecoli_fn
    for i in range(len(s) - len(top) + 1):
        mm = sum(b not in IUPAC[c] for b, c in zip(s[i:i + len(top)], top))
        if mm < best[2]:
            best = (i, i + len(top) - 1, mm)
    return best


# --------------------------------------------------------------------------
# 平铺方案优化
# --------------------------------------------------------------------------
def v_coverage(amplicons, inner_only=True):
    """每个可变区被扩增子（去掉引物后的内部）覆盖的比例"""
    res = {}
    for v, (a, b) in V_REGIONS.items():
        covered = np.zeros(b - a + 1, bool)
        for f_s, f_e, r_s, r_e in amplicons:
            lo, hi = (f_e + 1, r_s - 1) if inner_only else (f_s, r_e)
            for p in range(max(a, lo), min(b, hi) + 1):
                covered[p - a] = True
        res[v] = covered.mean()
    return res


def candidate_amplicons(F_sites, R_sites, amp_range=(180, 250), score_col="key_min", special_range=None):
    """
    所有满足长度约束的 (F, R) 组合；长度含引物。
    special_range: {可变区: (min, max)}，完整覆盖该可变区的扩增子使用单独的长度范围（例如 V4 放宽到 292 bp）
    """
    special_range = special_range or {}
    out = []
    for f in F_sites.itertuples():
        for r in R_sites.itertuples():
            L = r.end - f.start + 1
            if r.start <= f.end:
                continue
            lo_, hi_ = amp_range
            for v, (a_, b_) in special_range.items():
                x, y = V_REGIONS[v]
                if f.end < x and r.start > y:
                    lo_, hi_ = min(lo_, a_), max(hi_, b_)
            if lo_ <= L <= hi_ and (L <= amp_range[1] or any(
                    f.end < V_REGIONS[v][0] and r.start > V_REGIONS[v][1] for v in special_range)):
                out.append({"F": f, "R": r, "len": L, "lo": f.start, "hi": r.end,
                            "inner": (f.end + 1, r.start - 1),
                            "minc": min(getattr(f, score_col), getattr(r, score_col))})
    return out


def _amp_score(a, covered, weights):
    """扩增子内部（去引物）对“尚未覆盖”的可变区位置的加权覆盖"""
    s = 0.0
    lo, hi = a["inner"]
    for v, (x, y) in V_REGIONS.items():
        n = y - x + 1
        seg = covered[v][max(x, lo) - x: min(y, hi) - x + 1] if min(y, hi) >= max(x, lo) else []
        if len(seg):
            s += weights[v] * (len(seg) - seg.sum()) / n
    return s


def _dp_one_pool(amps, covered, weights, min_gap, must, tie=0.02):
    """
    带权区间调度（同一管内扩增子区间互不相交、间隔 ≥ min_gap）。
    状态带一个比特：是否已包含完整覆盖 must 中可变区的扩增子。
    分数 = Σ 新增加权覆盖 + tie × 引物最低覆盖率（用于并列时偏好更通用的引物）
    """
    amps = sorted(amps, key=lambda a: a["hi"])
    n = len(amps)
    his = [a["hi"] for a in amps]
    import bisect
    prev = [bisect.bisect_right(his, a["lo"] - min_gap - 1) - 1 for a in amps]
    gain = [_amp_score(a, covered, weights) for a in amps]
    keep = [g > 1e-9 for g in gain]                 # 不带来新覆盖的扩增子不入选
    amps = [a for a, k in zip(amps, keep) if k]
    gain = [g for g, k in zip(gain, keep) if k]
    n = len(amps)
    his = [a["hi"] for a in amps]
    prev = [bisect.bisect_right(his, a["lo"] - min_gap - 1) - 1 for a in amps]
    val = [g + tie * a["minc"] for g, a in zip(gain, amps)]

    def has_must(a):
        lo, hi = a["inner"]
        return all(lo <= V_REGIONS[m][0] and V_REGIONS[m][1] <= hi for m in must)

    mflag = [has_must(a) for a in amps]
    NEG = -1e18
    # best[i][k]: 前 i 个扩增子中的最优，k=1 表示已含 must
    best = [[0.0, NEG]]
    choice = [[None, None]]
    for i in range(n):
        b0, b1 = best[i]
        c0, c1 = (choice[i][0], choice[i][1])
        p = prev[i] + 1
        t0 = best[p][0] + val[i]
        t1 = best[p][1] + val[i]
        new0, new1 = b0, b1
        ch0, ch1 = c0, c1
        if not mflag[i] and t0 > new0:
            new0, ch0 = t0, ("take", i, 0)
        k_from0 = 1 if mflag[i] else 0
        if mflag[i] and t0 > new1:
            new1, ch1 = t0, ("take", i, 0)
        if t1 > new1:
            new1, ch1 = t1, ("take", i, 1)
        if ch0 is c0:
            ch0 = ("skip", i, 0)
        if ch1 is c1:
            ch1 = ("skip", i, 1)
        best.append([new0, new1])
        choice.append([ch0, ch1])
    k = 1 if must else 0
    if not must:
        k = 0 if best[n][0] >= best[n][1] else 1
    if best[n][k] <= NEG / 2:
        return []
    sol, i = [], n
    while i > 0:
        act, idx, kk = choice[i][k]
        if act == "skip":
            i -= 1
            continue
        sol.append(amps[idx])
        k = kk
        i = prev[idx] + 1
    return sol[::-1]


def optimize_tiling(F_sites, R_sites, amp_range=(180, 250), min_gap=0, must=("V4",), weights=None,
                    n_pools=1, rounds=4, score_col="key_min", special_range=None):
    """
    F_sites / R_sites: 已按覆盖率阈值筛选的候选位点（scan_sites 的输出）。

    同一管内的约束：扩增子区间 [F.start, R.end] 两两不相交，间隔 ≥ min_gap。
      * 若下一个扩增子的正向引物落在上一个扩增子内部（区间重叠），F(n+1) 与 R(n) 会扩出短产物，
        短产物扩增效率最高，会大量消耗这两条引物 —— 这就是要避免的“重复区域引物消耗”。
      * 区间不相交时，夹在中间的 R(n) 与 F(n+1) 3′ 端背向，不形成产物；只剩外侧 F(n)+R(n+1) 的长产物（>360 bp，实测极少）。
    单管：带权区间调度的动态规划，给出全局最优。
    两管：交替优化（固定一管，对另一管在“尚未覆盖的位置”上做 DP），收敛到局部最优；管与管之间允许重叠。
    """
    weights = weights or {v: 1.0 for v in V_REGIONS}
    amps = candidate_amplicons(F_sites, R_sites, amp_range, score_col, special_range)

    def cover(sol_list):
        covered = {v: np.zeros(y - x + 1, bool) for v, (x, y) in V_REGIONS.items()}
        for a in sol_list:
            lo, hi = a["inner"]
            for v, (x, y) in V_REGIONS.items():
                s, e = max(x, lo), min(y, hi)
                if e >= s:
                    covered[v][s - x:e - x + 1] = True
        return covered

    pools = [[] for _ in range(n_pools)]
    for _ in range(rounds if n_pools > 1 else 1):
        for k in range(n_pools):
            others = [a for j, p in enumerate(pools) if j != k for a in p]
            req = must if (k == 0) else ()
            pools[k] = _dp_one_pool(amps, cover(others), weights, min_gap, req)
    allsol = [(a, k) for k, p in enumerate(pools) for a in p]
    allsol.sort(key=lambda x: x[0]["lo"])
    vc = {v: c.mean() for v, c in cover([a for a, _ in allsol]).items()}
    for a, _ in allsol:
        a["vc"] = {v: c.mean() for v, c in cover([a]).items()}
    return {"sol": allsol, "vc": vc, "score": sum(weights[v] * c for v, c in vc.items())}


def summarize(sol):
    rows = []
    for k, (a, pl) in enumerate(sol, 1):
        covered = [v for v, c in a["vc"].items() if c > 0]
        rows.append({"amplicon": f"A{k}", "pool": "AB"[pl], "F site": f"{a['F'].start}-{a['F'].end}", "F primer": a["F"].primer,
                     "F key min": a["F"].key_min, "R site": f"{a['R'].start}-{a['R'].end}",
                     "R primer": a["R"].primer, "R key min": a["R"].key_min, "length": a["len"],
                     "variable regions": "+".join(f"{v}({a['vc'][v]:.0%})" if a['vc'][v] < .999 else v for v in covered)})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# 多重体系的引物相容性检查
# --------------------------------------------------------------------------
def expand(seq):
    return ["".join(p) for p in itertools.product(*[IUPAC[c] for c in seq])]


def tm_range(seq, Na=50, Mg=2.0, dNTPs=0.2, dnac1=250):
    """所有简并展开版本的最近邻 Tm（Biopython Tm_NN, SantaLucia 2004 参数，含 Mg2+ 校正）"""
    from Bio.SeqUtils import MeltingTemp as mt
    v = [mt.Tm_NN(s, nn_table=mt.DNA_NN4, Na=Na, Mg=Mg, dNTPs=dNTPs, dnac1=dnac1, dnac2=0) for s in expand(seq)]
    return min(v), max(v)


def gc_content(seq):
    vs = expand(seq)
    return np.mean([(s.count("G") + s.count("C")) / len(s) for s in vs])


def three_prime_dimer(p1, p2, min_len=4):
    """
    两条引物 3′ 端是否能互相配对（可被聚合酶延伸的二聚体）：
    p1 的 3′ 端 k 个碱基与 p2 任意位置的反向互补完全配对，返回最长 k（任一简并版本组合）
    """
    best = 0
    for a in expand(p1)[:16]:
        for b in expand(p2)[:16]:
            rb = revcomp(b)
            for k in range(min_len, min(len(a), len(b)) + 1):
                if a[-k:] in rb:
                    best = max(best, k)
    return best


def pool_report(primers):
    """primers: list of (name, seq)；返回每条引物的 Tm/GC，以及两两 3′ 端二聚体"""
    rows = []
    for n, s in primers:
        lo, hi = tm_range(s)
        rows.append({"primer": n, "seq": s, "len": len(s), "fold": len(expand(s)),
                     "GC": round(gc_content(s), 2), "Tm min": round(lo, 1), "Tm max": round(hi, 1)})
    dimers = []
    for (n1, s1), (n2, s2) in itertools.combinations_with_replacement(primers, 2):
        k = max(three_prime_dimer(s1, s2), three_prime_dimer(s2, s1))
        if k >= 4:
            dimers.append({"pair": f"{n1} × {n2}", "3' complementary nt": k})
    return pd.DataFrame(rows), pd.DataFrame(dimers)


def tune_tm(R: Reference, site_row, target_tm, min_len=17, max_len=24, max_degen=3, max_fold=8,
            max_loss=0.03, domain="Bacteria", key_phyla=KEY_PHYLA):
    """
    固定 3′ 端，在 5′ 端缩短或延长引物（min_len–max_len）以接近目标 Tm；
    主要门最低覆盖率与平均覆盖率的下降都不超过 max_loss。
    返回 dict(start, end, primer, key_mean, key_min, coverage, Tm_lo, Tm_hi)
    """
    sel, ph, w = _weights(R, domain, key_phyla)
    best = None
    for L in range(min_len, max_len + 1):
        if site_row["orient"] == "F":
            s, e, side = site_row["end"] - L + 1, site_row["end"], "right"
        else:
            s, e, side = site_row["start"], site_row["start"] + L - 1, "left"
        M = site_matrix(R, s, L)
        if M is None:
            continue
        M = M[sel]
        top, _ = best_degenerate(M, side, max_degen, max_fold, weights=w)
        if top is None:
            continue
        data = (M != 0).all(1) & (M != ord("-")).all(1)
        hit = _match(M[data], top, side)
        phd = ph[data]
        per = [hit[phd == k].mean() for k in key_phyla if (phd == k).sum() >= 10]
        kmean, kmin = float(np.mean(per)), float(min(per))
        if kmean < site_row["key_mean"] - max_loss or kmin < site_row["key_min"] - max_loss:
            continue
        primer = top if site_row["orient"] == "F" else revcomp(top)
        lo, hi = tm_range(primer)
        d = abs((lo + hi) / 2 - target_tm)
        if best is None or d < best[0] - 0.3:
            best = (d, dict(start=s, end=e, primer=primer, key_mean=kmean, key_min=kmin, coverage=float(hit.mean()),
                            Tm_lo=round(lo, 1), Tm_hi=round(hi, 1)))
    return best[1] if best else None


def optimize_pools(F_sites, R_sites, n_pools=2, amp_range=(180, 250), min_gap=20, must=("V4",), weights=None,
                   score_col="key_min", special_range=None, per_key=4):
    """
    精确求解 n_pools（1–3）管的最优平铺：
      * 同一管内扩增子区间互不相交且间隔 ≥ min_gap（避免短的跨区产物与引物消耗）
      * 不同管之间可以重叠
    按起点排序后做动态规划，状态 = 每一管最后一个扩增子。由于同管扩增子互不相交，
    新扩增子只可能与“其他管的最后一个扩增子”重叠，因此新增覆盖可以精确计算。
    目标：Σ_V 权重 × 覆盖比例（按碱基计），并列时偏好引物最低覆盖率更高的方案。
    """
    weights = weights or {v: 1.0 for v in V_REGIONS}
    amps = candidate_amplicons(F_sites, R_sites, amp_range, score_col, special_range)
    wpos = np.zeros(1600)
    for v, (x, y) in V_REGIONS.items():
        wpos[x:y + 1] = weights[v] / (y - x + 1)
    # 剪枝：相同覆盖“签名”的扩增子只保留引物最好的几个
    for a in amps:
        lo, hi = a["inner"]
        a["w"] = wpos[lo:hi + 1].sum()
        a["sig"] = tuple(v for v, (x, y) in V_REGIONS.items() if lo <= x and y <= hi) + \
            tuple(f"{v}~" for v, (x, y) in V_REGIONS.items() if not (lo <= x and y <= hi) and hi >= x and lo <= y)
    amps = [a for a in amps if a["w"] > 0]
    by = {}
    for a in amps:
        by.setdefault(a["sig"], []).append(a)
    pool = []
    for lst in by.values():
        lst.sort(key=lambda a: (-a["minc"], a["len"]))
        pool.extend(lst[:per_key])
    pool.sort(key=lambda a: (a["lo"], a["hi"]))
    n = len(pool)

    def overlap_w(a, b):
        if b is None:
            return 0.0
        lo, hi = max(a["inner"][0], b["inner"][0]), min(a["inner"][1], b["inner"][1])
        return wpos[lo:hi + 1].sum() if hi >= lo else 0.0

    def has_must(a):
        return all(a["inner"][0] <= V_REGIONS[m][0] and V_REGIONS[m][1] <= a["inner"][1] for m in must)

    # 状态: (tuple(各管最后扩增子下标, -1 表示空), must_flag) -> (score, minc, path)
    start = (tuple([-1] * n_pools), False)
    states = {start: (0.0, 1.0, ())}
    for i, a in enumerate(pool):
        new_states = dict(states)
        for (lasts, flag), (sc, mc, path) in states.items():
            for t in range(n_pools):
                lt = lasts[t]
                if lt >= 0 and not (pool[lt]["hi"] + min_gap < a["lo"]):
                    continue
                # 与其他管最后一个扩增子的重叠（同管已不相交）
                ov = 0.0
                others = [pool[lasts[u]] for u in range(n_pools) if u != t and lasts[u] >= 0]
                # 多于一个其他管时，取并集的近似：逐个扣除（3 管时可能略低估增益）
                cov_mask_lo, cov_mask_hi = a["inner"]
                seg = wpos[cov_mask_lo:cov_mask_hi + 1].copy()
                for b in others:
                    lo, hi = max(a["inner"][0], b["inner"][0]), min(a["inner"][1], b["inner"][1])
                    if hi >= lo:
                        seg[lo - cov_mask_lo:hi - cov_mask_lo + 1] = 0
                gain = seg.sum()
                if gain <= 1e-9:
                    continue
                nl = list(lasts)
                nl[t] = i
                # 管的编号无序：规范化为排序后的元组，减少重复状态
                key = (tuple(sorted(nl)), flag or has_must(a))
                val = (sc + gain, min(mc, a["minc"]), path + ((i, t),))
                old = new_states.get(key)
                if old is None or (round(val[0], 4), val[1]) > (round(old[0], 4), old[1]):
                    new_states[key] = val
        states = new_states
    cands = [(v, k) for k, v in states.items() if (k[1] or not must)]
    if not cands:
        return {"sol": [], "vc": {}, "score": 0}
    (score, minc, path), _ = max(cands, key=lambda x: (round(x[0][0], 4), x[0][1]))
    chosen = [pool[i] for i, _ in path]
    # 重新分配管号（区间图着色：按起点贪心）
    ends = [-10 ** 9] * n_pools
    sol = []
    for a in sorted(chosen, key=lambda a: a["lo"]):
        for t in range(n_pools):
            if ends[t] + min_gap < a["lo"]:
                ends[t] = a["hi"]
                sol.append((a, t))
                break
    covered = {v: np.zeros(y - x + 1, bool) for v, (x, y) in V_REGIONS.items()}
    for a, _ in sol:
        lo, hi = a["inner"]
        a["vc"] = {}
        for v, (x, y) in V_REGIONS.items():
            s, e = max(x, lo), min(y, hi)
            m = np.zeros(y - x + 1, bool)
            if e >= s:
                m[s - x:e - x + 1] = True
            covered[v] |= m
            a["vc"][v] = m.mean()
    vc = {v: c.mean() for v, c in covered.items()}
    return {"sol": sol, "vc": vc, "score": sum(weights[v] * c for v, c in vc.items())}
