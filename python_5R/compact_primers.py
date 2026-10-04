"""
compact_primers.py -- 给定扩增子的引物位点，在覆盖率目标下找最少的寡核苷酸，并均衡 Tm、避开二聚体

覆盖率 = 12 个主要细菌门等权的平均覆盖率（最多 1 个错配、3′ 端 3 个碱基必须匹配）。
每个位点：3′ 端位置允许 ±shift nt 微调，长度 17–22 nt（5′ 端伸缩），用贪心集合覆盖设计 k=1,2,... 条简并引物，
取满足目标的最小 k。再在各位点的候选配置中做坐标下降，使 Tm 接近目标、严重二聚体最少、寡核苷酸最少。
"""
import itertools
import numpy as np
import primer3
import primer_design as pdz
import iterate_5R as it

PARAMS = dict(mv_conc=50, dv_conc=2, dntp_conc=0.2, dna_conc=250)
_tm_cache, _dimer_cache, _hp_cache = {}, {}, {}


def tm(seq):
    if seq not in _tm_cache:
        _tm_cache[seq] = primer3.calc_tm(seq, **PARAMS)
    return _tm_cache[seq]


def oligo_tm(p):
    v = [tm(s) for s in pdz.expand(p)]
    return float(np.mean(v)), float(min(v)), float(max(v))


def dimer(a, b, temp=60.0):
    key = (a, b) if a <= b else (b, a)
    if key not in _dimer_cache:
        h = primer3.calc_heterodimer(a, b, temp_c=temp, **PARAMS).dg / 1000
        e1 = primer3.bindings.calc_end_stability(a, b, temp_c=temp, **PARAMS).dg / 1000
        e2 = primer3.bindings.calc_end_stability(b, a, temp_c=temp, **PARAMS).dg / 1000
        _dimer_cache[key] = (h, min(e1, e2))
    return _dimer_cache[key]


def hairpin(a, temp=60.0):
    if a not in _hp_cache:
        _hp_cache[a] = primer3.calc_hairpin(a, temp_c=temp, **PARAMS).dg / 1000
    return _hp_cache[a]


def coverage(M, ph, tops, side, key=pdz.KEY_PHYLA, minn=10):
    h = it.hit_any(M, tops, side)
    per = [h[ph == k].mean() for k in key if (ph == k).sum() >= minn]
    return float(np.mean(per)), float(h.mean()), float(np.min(per))


def coverage_blocks(M, ph, tops, side, blocks):
    """每个参考库（块）各自的主要门平均覆盖率；返回 (最小值, 各块的值)"""
    h = it.hit_any(M, tops, side)
    vals = []
    for sl in blocks:
        per = [h[sl][ph[sl] == k].mean() for k in pdz.KEY_PHYLA if (ph[sl] == k).sum() >= 10]
        vals.append(float(np.mean(per)))
    return min(vals), vals


def design_min_set(M, ph, side, target, blocks=None, kmax=6, folds=((4, 2), (8, 3))):
    """
    返回满足 target 的最小 k 的引物集合（顶链方向）；同 k 下取总简并数小者。
    blocks: 各参考库在 M 中的切片；每个库的覆盖率都要 >= target（防止混合后被大库掩盖）。
    贪心集合是嵌套的：每种简并上限只设计一次到 kmax 条，再取前缀判断。
    """
    blocks = blocks or [slice(0, len(M))]
    best = None
    for fold, deg in folds:
        tops = it.design_set(M, ph, side, kmax, fold, deg)
        for k in range(1, len(tops) + 1):
            worst, vals = coverage_blocks(M, ph, tops[:k], side, blocks)
            if worst >= target:
                tot = sum(len(pdz.expand(x)) for x in tops[:k])
                if best is None or (k, tot) < (len(best[0]), best[1]):
                    c = coverage(M, ph, tops[:k], side)
                    best = (tops[:k], tot, (worst,) + tuple(vals) + (c[2],))
                break
    return best


def site_geometry(orient, p3, L):
    """返回 (start, L)：F 的 3′ 端在 start+L-1；R 的 3′ 端在 start"""
    return (p3 - L + 1, L) if orient == "F" else (p3, L)


def enumerate_configs(site, design_refs, target, shifts=range(-3, 4), lengths=range(17, 23)):
    """site = (name, orient, st, L)。返回候选配置列表（每个都满足覆盖率目标）"""
    name, orient, st, L0 = site
    p3_0 = st + L0 - 1 if orient == "F" else st
    side = "right" if orient == "F" else "left"
    out = []
    for sh in shifts:
        for L in lengths:
            s0, L_ = site_geometry(orient, p3_0 + sh, L)
            Md = [it.site_data(r, s0, L_) for r in design_refs]
            M = np.vstack([x[0] for x in Md]); ph = np.concatenate([x[1] for x in Md])
            if len(M) < 200:
                continue
            cuts = np.cumsum([0] + [len(x[0]) for x in Md])
            blocks = [slice(cuts[i], cuts[i + 1]) for i in range(len(Md))]
            res = design_min_set(M, ph, side, target, blocks)
            if res is None:
                continue
            tops, tot, c = res
            prim = [t if orient == "F" else pdz.revcomp(t) for t in tops]
            tms = [oligo_tm(p)[0] for p in prim]
            out.append(dict(name=name, orient=orient, start=s0, L=L_, shift=sh, tops=tops, prim=prim, n=len(prim), fold=tot,
                            cov_design=c, tm_mean=float(np.mean(tms)), tm_dev=max(abs(t - 0) for t in tms)))
    return out


def pool_issues(prim_by_site, dg_pair=-9.0, dg_end=-6.0, dg_hp=-3.0):
    """所有寡核苷酸（含简并展开）两两检查；返回严重二聚体列表和发夹列表。键为寡核苷酸（含简并码）"""
    oligos = [(site, p) for site, ps in prim_by_site.items() for p in ps]
    exp = {(s, p): pdz.expand(p) for s, p in oligos}
    severe, hairp = [], []
    for (i, (s1, p1)), (j, (s2, p2)) in itertools.combinations_with_replacement(list(enumerate(oligos)), 2):
        w_h, w_e = 0.0, 0.0
        for a in exp[(s1, p1)]:
            for b in exp[(s2, p2)]:
                h, e = dimer(a, b)
                w_h, w_e = min(w_h, h), min(w_e, e)
        if w_h <= dg_pair or w_e <= dg_end:
            severe.append((s1, p1, s2, p2, round(w_h, 1), round(w_e, 1)))
    for s, p in oligos:
        w = min(hairpin(a) for a in exp[(s, p)])
        if w <= dg_hp:
            hairp.append((s, p, round(w, 1)))
    return severe, hairp
