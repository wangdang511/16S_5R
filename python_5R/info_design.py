"""
info_design.py -- 以“每条序列平均测到的信息量”为目标的单管多扩增子设计

信息量 = 被测序碱基位置的保守度熵之和（bit，熵来自细菌序列的位点碱基分布）。
一个扩增子的期望信息量 = 其信息量 × 序列被该扩增子扩出的概率（正反引物都匹配）；
方案的期望信息量 = 各扩增子之和，所以可以在“同管内扩增子区间互不相交”的约束下用动态规划精确求最优。
"""
import bisect
import numpy as np
import pandas as pd
import primer_design as pdz


def seq_weights(phylum, key=pdz.KEY_PHYLA):
    """主要门等权，其余细菌门合计 1 份"""
    lab = np.where(np.isin(phylum, key), phylum, "__other__")
    w = np.zeros(len(phylum))
    for g in np.unique(lab):
        m = lab == g
        w[m] = 1.0 / m.sum()
    return w / w.sum()


def sequenced_positions(fs, fl, rs, rl, read_len):
    """正向位点 [fs, fs+fl-1]，反向位点 [rs, rs+rl-1]；read_len 含引物；返回被测序的内部位置"""
    fe, re_ = fs + fl - 1, rs + rl - 1
    fwd = range(fe + 1, min(fs + read_len - 1, rs - 1) + 1)
    rev = range(max(re_ - read_len + 1, fe + 1), rs)
    return sorted(set(fwd) | set(rev))


def site_hits(ref, seq_top, start, side, rows):
    M = pdz.site_matrix(ref, start, len(seq_top))[rows]
    ok = (M != 0).all(1) & (M != ord("-")).all(1)
    h = np.zeros(len(M), bool)
    h[ok] = pdz._match(M[ok], seq_top, side)
    return h


def build_candidates(sites, refs, entropy, read_len=126, amp_range=(180, 250), min_pair_cov=0.5, min_gap_end=0):
    """
    sites: DataFrame(start, end, len, orient, primer)；refs: [(Reference, rows, weights)] 设计集（可多个，权重按各库样本数归一后合并）
    返回候选扩增子列表 dict(F, R, lo, hi, len, bits, cov, value, seq_cov)
    """
    # 每个 (start, orient) 只保留引物质量最好的长度版本
    s = sites.sort_values("key_min", ascending=False).drop_duplicates(["start", "orient"])
    F = s[s.orient == "F"].sort_values("start").reset_index(drop=True)
    R = s[s.orient == "R"].sort_values("start").reset_index(drop=True)

    def hit_matrix(df, orient):
        out = []
        for ref, rows, _ in refs:
            H = np.zeros((len(df), len(rows)), bool)
            for i, r in enumerate(df.itertuples()):
                top = r.primer if orient == "F" else pdz.revcomp(r.primer)
                H[i] = site_hits(ref, top, r.start, "right" if orient == "F" else "left", rows)
            out.append(H)
        return np.hstack(out)

    HF, HR = hit_matrix(F, "F"), hit_matrix(R, "R")
    w = np.concatenate([wt for _, _, wt in refs])
    w = w / w.sum()
    r_end = (R.end).values
    cands = []
    for i, f in enumerate(F.itertuples()):
        lo_end = f.start + amp_range[0] - 1
        hi_end = f.start + amp_range[1] - 1
        j0 = bisect.bisect_left(list(r_end), lo_end)
        j1 = bisect.bisect_right(list(r_end), hi_end)
        if j1 <= j0:
            continue
        hf = HF[i] * w
        cov = HR[j0:j1].astype(np.float64) @ hf
        for k, c in enumerate(cov):
            if c < min_pair_cov:
                continue
            r = R.iloc[j0 + k]
            if r.start <= f.end:
                continue
            pos = sequenced_positions(f.start, f.len, r.start, r.len, read_len)
            bits = float(entropy.reindex(pos).fillna(0).sum())
            cands.append(dict(F=f._asdict(), R=r.to_dict(), lo=int(f.start), hi=int(r.end), len=int(r.end - f.start + 1),
                              bits=bits, cov=float(c), value=bits * float(c), nbases=len(pos), fi=i, rj=j0 + k))
    return cands, (F, R, HF, HR, w)


def best_tiling(cands, K, min_gap=20):
    """动态规划：至多 K 个扩增子，区间互不相交（间隔 ≥ min_gap），最大化 Σ 期望信息量"""
    cs = sorted(cands, key=lambda a: a["hi"])
    his = [a["hi"] for a in cs]
    prev = [bisect.bisect_right(his, a["lo"] - min_gap - 1) - 1 for a in cs]
    n = len(cs)
    best = [[0.0] * (K + 1) for _ in range(n + 1)]
    take = [[False] * (K + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        a = cs[i - 1]
        p = prev[i - 1] + 1
        for k in range(K + 1):
            best[i][k] = best[i - 1][k]
            if k >= 1:
                v = best[p][k - 1] + a["value"]
                if v > best[i][k] + 1e-12:
                    best[i][k] = v
                    take[i][k] = True
    sol, i, k = [], n, K
    while i > 0 and k >= 0:
        if take[i][k]:
            sol.append(cs[i - 1]); k -= 1; i = prev[i - 1] + 1
        else:
            i -= 1
    return sol[::-1], best[n][K]


def best_tiling_forced(cands, K, forced, min_gap=20):
    """必须包含 forced（一个候选扩增子）的最优平铺：左右两侧各自做 DP，枚举 k 的分配"""
    left = [a for a in cands if a["hi"] + min_gap < forced["lo"]]
    right = [a for a in cands if a["lo"] > forced["hi"] + min_gap]
    best = (-1, None)
    for k1 in range(K):
        k2 = K - 1 - k1
        sl, vl = best_tiling(left, k1) if (k1 > 0 and left) else ([], 0.0)
        sr, vr = best_tiling(right, k2) if (k2 > 0 and right) else ([], 0.0)
        v = vl + vr + forced["value"]
        if v > best[0]:
            best = (v, sl + [forced] + sr)
    return best[1], best[0]


def v4_fraction(a, read_len, v4=(576, 682)):
    pos = set(sequenced_positions(a["F"]["start"], a["F"]["len"], a["R"]["start"], a["R"]["len"], read_len))
    return sum(p in pos for p in range(v4[0], v4[1] + 1)) / (v4[1] - v4[0] + 1)
