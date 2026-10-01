"""
iterate_5R.py -- 在现有 5R 的 10 个引物位点上做简并碱基迭代，以及在同一位点放多条引物（贪心集合覆盖）

位置和长度保持不变（3′ 端位置不动、长度不变），所以扩增子定义、SMURF 的引物切除长度、k-mer 数据库的截取方式都不用改。
设计集 = Greengenes 85% OTU（4,797）+ SILVA 随机 12,000 条；验证 = SILVA 其余 48,000 条 + Greengenes 12,000 条留出序列。
匹配规则：最多 1 个错配，且 3′ 端 3 个碱基必须匹配。
"""
import itertools, json, pickle, numpy as np, pandas as pd, warnings
import primer_design as pdz
warnings.filterwarnings("ignore")

FIVE = [("R1-F", "TGGCGAACGGGTGAGTAA", "F", 103), ("R1-R", "CCGTGTCTCAGTCCCARTG", "R", 314),
        ("R2-F", "ACTCCTACGGGAGGCAGC", "F", 338), ("R2-R", "GTATTACCGCGGCTGCTG", "R", 519),
        ("R3-F", "GTGTAGCGGTGRAATGCG", "F", 685), ("R3-R", "CCCGTCAATTCMTTTGAGTT", "R", 908),
        ("R4-F", "GGAGCATGTGGWTTAATTCGA", "F", 944), ("R4-R", "CGTTGCGGGACTTAACCC", "R", 1087),
        ("R5-F", "GGAGGAAGGTGGGGATGAC", "F", 1175), ("R5-R", "AAGGCCCGGGAACGTATT", "R", 1374)]


def site_data(ref, start, L, domain="Bacteria"):
    M = pdz.site_matrix(ref, start, L)
    sel = (ref.domain == domain)
    ok = (M != 0).all(1) & (M != ord("-")).all(1) & sel
    return M[ok], ref.phylum[ok], np.flatnonzero(ok)


def hit_any(M, primers_top, side):
    h = np.zeros(len(M), bool)
    for p in primers_top:
        h |= pdz._match(M, p, side)
    return h


def design_set(M, ph, side, k, max_fold, max_degen, key=pdz.KEY_PHYLA):
    """贪心：每一轮针对尚未被覆盖的序列设计一条简并引物（主要门等权）"""
    lab = np.where(np.isin(ph, key), ph, "__other__")
    w = np.zeros(len(M))
    for g in np.unique(lab):
        m = lab == g
        w[m] = 1.0 / m.sum()
    covered = np.zeros(len(M), bool)
    out = []
    for _ in range(k):
        wi = w * (~covered)
        if wi.sum() <= 0:
            break
        top, _ = pdz.best_degenerate(M, side, max_degen, max_fold, weights=wi)
        if top is None:
            break
        h = pdz._match(M, top, side)
        if (h & ~covered).sum() == 0:
            break
        out.append(top)
        covered |= h
    return out


def metrics(M, ph, primers_top, side, key=pdz.KEY_PHYLA, minn=30):
    h = hit_any(M, primers_top, side)
    per = {k: h[ph == k].mean() for k in key if (ph == k).sum() >= minn}
    v = np.array(list(per.values()))
    return dict(all_bact=float(h.mean()), key_mean=float(v.mean()), key_min=float(v.min()), worst=min(per, key=per.get))


def run(R, S, H, out_prefix):
    rng = np.random.default_rng(3)
    idx12 = np.sort(rng.choice(len(S.ids), 12000, replace=False))      # 与 SILVA 扫描用的子集相同
    rest = np.setdiff1d(np.arange(len(S.ids)), idx12)
    S12 = pdz.Reference(S.aln[idx12], [S.ids[i] for i in idx12], S.domain[idx12], S.phylum[idx12], S.col_of_pos, S.ref_seq, {})
    Srest = pdz.Reference(S.aln[rest], [S.ids[i] for i in rest], S.domain[rest], S.phylum[rest], S.col_of_pos, S.ref_seq, {})
    rows, sets = [], {}
    for name, seq, o, st in FIVE:
        L = len(seq)
        side = "right" if o == "F" else "left"
        Md = [site_data(r, st, L) for r in (R, S12)]
        M = np.vstack([m[0] for m in Md]); ph = np.concatenate([m[1] for m in Md])
        tests = {"SILVA(48k)": site_data(Srest, st, L), "GG留出(12k)": site_data(H, st, L)}
        top0 = [seq if o == "F" else pdz.revcomp(seq)]
        # 原引物（含简并碱基）本身
        configs = {"原引物": top0}
        configs["1条 ≤4简并"] = design_set(M, ph, side, 1, 4, 2)
        configs["1条 ≤16简并"] = design_set(M, ph, side, 1, 16, 4)
        configs["2条 各≤4"] = design_set(M, ph, side, 2, 4, 2)
        configs["3条 各≤4"] = design_set(M, ph, side, 3, 4, 2)
        for cn, tops in configs.items():
            prim = [t if o == "F" else pdz.revcomp(t) for t in tops]
            fold = sum(len(pdz.expand(p)) for p in prim)
            r = dict(primer=name, site=f"{st}-{st + L - 1}", config=cn, n_oligos=len(prim), total_fold=fold, seqs=" | ".join(prim))
            r.update({f"design_{k}": v for k, v in metrics(M, ph, tops, side).items()})
            for tn, (Mt, pht, _) in tests.items():
                r.update({f"{tn}_{k}": v for k, v in metrics(Mt, pht, tops, side).items()})
            rows.append(r)
            sets[(name, cn)] = prim
    df = pd.DataFrame(rows)
    df.to_csv(out_prefix + "_primers.csv", index=False)
    pickle.dump(sets, open(out_prefix + "_sets.pkl", "wb"))
    return df, sets, Srest, S12


if __name__ == "__main__":
    import sys
    SC = sys.argv[1]
    R = pickle.load(open(SC + "R.pkl", "rb")); S = pickle.load(open(SC + "S.pkl", "rb")); H = pickle.load(open(SC + "H.pkl", "rb"))
    df, sets, Srest, S12 = run(R, S, H, SC + "iter5R")
    pd.set_option("display.width", 250)
    cols = ["primer", "config", "n_oligos", "total_fold", "SILVA(48k)_key_mean", "SILVA(48k)_key_min", "GG留出(12k)_key_mean", "GG留出(12k)_key_min"]
    print(df[cols].round(3).to_string(index=False))
    print(df.groupby("config")[cols[4:]].mean().round(3))
