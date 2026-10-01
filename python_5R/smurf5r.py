"""
smurf5r.py -- 5R / SMURF 流程的 Python 逐函数复现
=====================================================

本模块与 ``mFiles/`` 下的 MATLAB 实现一一对应（函数名后的注释标出了对应的 .m 文件），
目的是让流程的每一步都可以在 Jupyter Notebook 中单独运行、检查中间结果。

设计原则：
* 默认行为尽量与 MATLAB 版本**逐位一致**（包括几个已知的小 bug，见 ``faithful`` 参数）。
* 计算量最大的“reads 与数据库 k-mer 比对”一步，MATLAB 使用暴力 Hamming 距离；
  这里用“鸽巢原理 + 分段哈希索引”实现**结果完全相同**但快几个数量级的检索。

依赖: numpy, scipy, pandas
"""
from __future__ import annotations

import glob
import os
import re
import time
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

EPS = np.finfo(float).eps

# --------------------------------------------------------------------------
# 0. 配置   (get_configs.m)
# --------------------------------------------------------------------------
PRIMERS = [
    # (forward, reverse)                              名称
    ("TGGCGAACGGGTGAGTAA", "CCGTGTCTCAGTCCCARTG"),    # R1: 143-16S_DS-F1b / 145-16S_DS-R1a  (V2)
    ("ACTCCTACGGGAGGCAGC", "GTATTACCGCGGCTGCTG"),     # R2: 130-16S_DS-F2  / 131-16S_DS-R2   (V3)
    ("GTGTAGCGGTGRAATGCG", "CCCGTCAATTCMTTTGAGTT"),   # R3: 132-16S_DS-F3  / 133-16S_DS-R3   (V4-V5)
    ("GGAGCATGTGGWTTAATTCGA", "CGTTGCGGGACTTAACCC"),  # R4: 134-16S_DS-F4  / 135-16S_DS-R4   (V6)
    ("GGAGGAAGGTGGGGATGAC", "AAGGCCCGGGAACGTATT"),    # R5: 149-16S_DS-F5b / 153-16S_DS-R5c  (V8)
]


@dataclass
class Config:
    """所有参数集中在一起；默认值与 get_configs.m 相同。"""
    kmer_len: int = 126          # 用于重建的读长（每端截取的长度）
    db_kmer_len: int = 160       # 数据库中预存的 k-mer 长度 (RL160)
    allowed_mm: int = 2          # 数据库构建时引物允许的错配数 (2mm)
    # --- 质量过滤 (PrepConfig) ---
    qual_th: int = 30
    prc_high_qual: float = 0.75
    low10_th: int = 3
    max_num_Ns: int = 0
    max_err_inprimer: int = 2
    # --- 算法 (AlgoConfig) ---
    min_read_freq: float = 1e-4
    min_read_count: int = 2
    pe: float = 0.005            # 每碱基测序错误率
    nMM_cut: int = 2             # read 与 k-mer 允许的最大错配
    tol: float = 5e-7
    num_iter: int = 10000
    do_filter: int = 1
    filter_included_bacteria: int = 1
    # True: 复现 MATLAB 中 R2 质量过滤误用 q1 的 bug（见 read_fastq_save_unireads.m）
    faithful: bool = True
    primers: list = field(default_factory=lambda: list(PRIMERS))

    def __post_init__(self):
        self.kmer_len = min(self.kmer_len, self.db_kmer_len)

    @property
    def primers_len(self):
        return np.array([[len(f), len(r)] for f, r in self.primers])

    def db_file_prefix(self, db_name="GreenGenes_201305", max_non_acgt=3):
        return (f"{db_name}_unique_up_to_{max_non_acgt}_ambiguous_16S"
                f"_ffpe5regions_{self.allowed_mm}mm_RL{self.db_kmer_len}")


# --------------------------------------------------------------------------
# 统计对象   (ReadsStats.m)
# --------------------------------------------------------------------------
class ReadsStats:
    def __init__(self, nR):
        self.nR = nR
        self.stats = {}            # 全局统计: name -> (unique, counts)
        self.region_stats = {}     # 区域统计: name -> {rr: (unique, counts)}

    def add(self, desc, n_unique, n_counts, rr=0):
        if rr == 0:
            self.stats[desc] = (n_unique, n_counts)
        else:
            self.region_stats.setdefault(desc, {})[rr] = (n_unique, n_counts)

    def to_series(self):
        out = {k: v[1] for k, v in self.stats.items()}
        for k, d in self.region_stats.items():
            for rr in range(1, self.nR + 1):
                out[f"{k}{rr}"] = d.get(rr, (0, 0))[1]
        return pd.Series(out)


# --------------------------------------------------------------------------
# 工具函数
# --------------------------------------------------------------------------
_IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
          "R": "AG", "Y": "CT", "K": "GT", "M": "AC", "S": "CG", "W": "AT",
          "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
_COMP = bytes.maketrans(b"ACGT", b"TGCA")


def expand_degenerate(seq: str, max_non_acgt=3):
    """把含简并碱基的引物展开为所有具体序列   (unambiguit_one_seq.m)"""
    seqs = [""]
    n_amb = sum(c not in "ACGT" for c in seq)
    if n_amb > max_non_acgt:
        return []
    for c in seq:
        seqs = [s + b for s in seqs for b in _IUPAC[c]]
    return seqs


def revcomp(seq: str) -> str:
    return seq.encode().translate(_COMP)[::-1].decode()


def to_u8(strings, width=None):
    """字符串列表 -> (N, L) uint8 矩阵"""
    if len(strings) == 0:
        return np.zeros((0, width or 0), np.uint8)
    width = width or len(strings[0])
    return np.frombuffer("".join(strings).encode(), np.uint8).reshape(-1, width).copy()


def u8_to_str(mat):
    return [r.tobytes().decode() for r in mat]


def unique_rows_with_counts(S, freq):
    """对矩阵行去重并累计频数 (MATLAB: sortrows + unique 'first'/'last' + cumsum)"""
    if S.shape[0] == 0:
        return S, freq
    v = np.ascontiguousarray(S).view(np.dtype((np.void, S.shape[1])))[:, 0]
    uniq, inv = np.unique(v, return_inverse=True)
    counts = np.bincount(inv.ravel(), weights=freq).astype(np.int64)
    U = np.frombuffer(uniq.tobytes(), np.uint8).reshape(len(uniq), S.shape[1])
    return U.copy(), counts


# --------------------------------------------------------------------------
# 1. 样本拆分 / 样本名    (split_files2directories.m, extract_sample_name.m)
# --------------------------------------------------------------------------
def extract_sample_name(file_name: str) -> str:
    """'RDB1_TTGGTGCA_L001_R1_001.fastq' -> 'RDB1_TTGGTGCA'"""
    m = re.match(r"^(.*?)_L0\d\d_R[123]_", file_name)
    return m.group(1) if m else ""


def find_sample_pairs(fastq_dir):
    """返回 {sample_name: [(R1, R2), ...]}；支持 .fastq / .fastq.gz"""
    files = sorted(glob.glob(os.path.join(fastq_dir, "**", "*.fastq*"), recursive=True))
    samples = {}
    for f in files:
        if f.endswith(".zip"):
            continue
        b = os.path.basename(f)
        if "_R1_" not in b:
            continue
        r2 = os.path.join(os.path.dirname(f), b.replace("_R1_", "_R2_"))
        if not os.path.exists(r2):
            raise FileNotFoundError(f"Couldn't match R1-R2 pair for {f}")
        samples.setdefault(extract_sample_name(b), []).append((f, r2))
    return dict(sorted(samples.items()))


# --------------------------------------------------------------------------
# 2. 质量过滤 + 去重     (read_fastq_save_unireads.m)
# --------------------------------------------------------------------------
def _read_fastq(path):
    import gzip
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        lines = fh.read().split("\n")
    seqs = lines[1::4]
    quals = lines[3::4]
    n = min(len(seqs), len(quals))
    return seqs[:n], quals[:n]


def _trim_block(seqs, quals, L):
    """长度 >= L 的 read 截到 L；过短 read 用 'N'/0 填充（之后会被过滤）"""
    n = len(seqs)
    long_enough = np.fromiter((len(s) >= L for s in seqs), bool, n)
    S = np.full((n, L), ord("N"), np.uint8)
    Q = np.zeros((n, L), np.int16)
    idx = np.flatnonzero(long_enough)
    if len(idx):
        S[idx] = to_u8([seqs[i][:L] for i in idx], L)
        Q[idx] = to_u8([quals[i][:L] for i in idx], L).astype(np.int16) - 33
    return S, Q, long_enough


def quality_filter_pairs(r1_path, r2_path, cfg: Config, stats: ReadsStats | None = None):
    """
    规则（逐 read）：
      * 两端读长均 >= kmer_len，截取前 kmer_len 个碱基
      * 无 N (max_num_Ns = 0)
      * Q>30 的碱基比例 >= 75%，且 Q<10 的碱基少于 3 个
    合格的 read 对拼接为 [R1 | R2]（R2 尚未反向互补），再去重计数。
    """
    L = cfg.kmer_len
    s1, q1 = _read_fastq(r1_path)
    s2, q2 = _read_fastq(r2_path)
    S1, Q1, long1 = _trim_block(s1, q1, L)
    S2, Q2, long2 = _trim_block(s2, q2, L)

    acgt = np.isin(S1, np.frombuffer(b"ACGT", np.uint8))
    nonambig_1 = (~acgt).sum(1) <= cfg.max_num_Ns
    nonambig_2 = (~np.isin(S2, np.frombuffer(b"ACGT", np.uint8))).sum(1) <= cfg.max_num_Ns

    good_1 = ((Q1 > cfg.qual_th).mean(1) > cfg.prc_high_qual - EPS) & ((Q1 < 10).sum(1) < cfg.low10_th)
    low10_src = Q1 if cfg.faithful else Q2       # MATLAB 原代码此处误用了 q1
    good_2 = ((Q2 > cfg.qual_th).mean(1) > cfg.prc_high_qual - EPS) & ((low10_src < 10).sum(1) < cfg.low10_th)

    keep = good_1 & good_2 & nonambig_1 & nonambig_2
    S = np.hstack([S1[keep], S2[keep]])
    Suni, freq = unique_rows_with_counts(S, np.ones(S.shape[0]))

    if stats is not None:
        stats.add("Number of loaded reads", np.nan, len(s1))
        stats.add("Number of long reads", np.nan, int((long1 & long2).sum()))
        stats.add("Number of good reads", len(freq), int(keep.sum()))
    return Suni, freq


# --------------------------------------------------------------------------
# 3. 按引物把 read 分配到 5 个区域     (read_unireads_save_split_to_regions.m)
# --------------------------------------------------------------------------
def split_to_regions(Suni, freq, cfg: Config, stats: ReadsStats | None = None):
    """
    * 把 R2 部分反向互补，使整条 [R1 | rc(R2)] 与 16S 正链同向
    * 区域 r 的判定：R1 开头与正向引物、序列末尾与反向引物的反向互补，
      各自错配数 <= max_err_inprimer(2)
    * 去掉两端引物后得到长度 2L - len(Fp) - len(Rp) 的“配对 read”
    """
    L = cfg.kmer_len
    S = Suni.copy()
    tail = S[:, L:]
    comp = tail.copy()
    for a, b in zip(b"ACGT", b"TGCA"):
        comp[tail == a] = b
    S[:, L:] = comp[:, ::-1]

    regions = []
    mapped_any = np.zeros(len(freq), bool)
    for rr, (fp, rp) in enumerate(cfg.primers, start=1):
        fps = to_u8(expand_degenerate(fp))
        rps = to_u8([revcomp(x) for x in expand_degenerate(rp)])
        flen, rlen = fps.shape[1], rps.shape[1]
        head, end = S[:, :flen], S[:, -rlen:]
        notN_h, notN_e = head != ord("N"), end != ord("N")
        d_fwd = np.min([((head != p) & notN_h).sum(1) for p in fps], axis=0)
        d_rvs = np.min([((end != p) & notN_e).sum(1) for p in rps], axis=0)
        use = (d_fwd <= cfg.max_err_inprimer) & (d_rvs <= cfg.max_err_inprimer)
        mapped_any |= use
        reads = S[use, flen:S.shape[1] - rlen]
        f = freq[use]
        if stats is not None:
            stats.add("Number of reads mapped to region", len(f), int(f.sum()), rr)
        regions.append({"reads": reads, "freq": f})
    pct_unique = 100 * sum(len(r["freq"]) for r in regions) / max(len(freq), 1)
    pct_counts = 100 * freq[mapped_any].sum() / max(freq.sum(), 1)
    print(f"Mapped to primers {pct_unique:.0f}% of unique reads, {pct_counts:.0f}% of read counts")
    return regions


# --------------------------------------------------------------------------
# 4. 载入 k-mer 数据库     (load_bact_DB.m)
# --------------------------------------------------------------------------
def load_region_db(db_dir, cfg: Config, rr: int, db_name="GreenGenes_201305"):
    """
    数据库文件（每个区域一个 .mat）中:
      values           (K x 320) 每个唯一扩增子的 [前160nt | 后160nt]（含引物，正链方向）
      indInValue       (nB x 1)  每条 16S 序列对应 values 的行号；0 表示该区域不被扩增
      is_perfect_match (nB x 1)  该序列与本区域引物是否 0 错配
    这里把 k-mer 截成与实验 read 相同的形状：去掉引物、每端只保留 kmer_len。
    """
    path = os.path.join(db_dir, f"{cfg.db_file_prefix(db_name)}_region{rr}.mat")
    d = sio.loadmat(path, variable_names=["values", "indInValue", "is_perfect_match"])
    v = d["values"]
    W = 2 * cfg.db_kmer_len
    values = np.frombuffer(v.astype(f"S{W}").tobytes(), np.uint8).reshape(len(v), W)
    ind = d["indInValue"].ravel().astype(np.int64)            # 1-based, 0 = 不扩增
    is_pm = d["is_perfect_match"].ravel().astype(bool)
    L = cfg.kmer_len
    pf, pr = cfg.primers_len[rr - 1]
    kmers = np.hstack([values[:, pf:L], values[:, W - L:W - pr]])
    U, inv = np.unique(kmers, axis=0, return_inverse=True)
    inv = inv.ravel()
    amplified = ind > 0
    ind_in_kmers = np.zeros(len(ind), np.int64)
    ind_in_kmers[amplified] = inv[ind[amplified] - 1] + 1
    # 只保留至少被一条序列使用的 k-mer，并重新编号（与 MATLAB 第二次 unique 一致）
    used, new_ind = np.unique(ind_in_kmers, return_inverse=True)
    if used[0] == 0:
        U = U[used[1:] - 1]
        new_ind = new_ind.ravel()          # 0 = 不扩增, 1..K
    else:                                  # MATLAB 在此情况下会丢掉第一个 k-mer（边界 bug）
        U = U[used - 1]
        new_ind = new_ind.ravel() + 1
    is_pm = is_pm & (new_ind > 0)
    return {"kmers": U, "indInSeqs": new_ind, "is_perfect_match": is_pm}


def load_bact_db(db_dir, cfg: Config, db_name="GreenGenes_201305", regions=(1, 2, 3, 4, 5)):
    dbs = []
    for rr in regions:
        t = time.time()
        dbs.append(load_region_db(db_dir, cfg, rr, db_name))
        print(f"Region {rr}: {dbs[-1]['kmers'].shape[0]:,} unique k-mers "
              f"(len {dbs[-1]['kmers'].shape[1]}), {(dbs[-1]['indInSeqs'] > 0).sum():,} amplified seqs "
              f"[{time.time() - t:.1f}s]")
    return dbs


def load_headers_taxonomy(db_dir, db_name="GreenGenes_201305"):
    h = sio.loadmat(os.path.join(db_dir, f"{db_name}_unique_up_to_3_ambiguous_16S_headers.mat"))["Header_uni"]
    headers = np.array([int(x.ravel()[0]) for x in h.ravel()])
    t = sio.loadmat(os.path.join(db_dir, "taxonomy_db.mat"))
    ranks = [str(x.ravel()[0]) for x in t["ranks_to_extract"].ravel()]
    tax = t["taxa_name_calls"]
    taxa = pd.DataFrame([[str(c.ravel()[0]) if c.size else "" for c in row] for row in tax], columns=ranks)
    return headers, taxa


# --------------------------------------------------------------------------
# 5. 构建 read x 细菌 的似然矩阵 A     (build_A_matrices.m)
# --------------------------------------------------------------------------
class KmerIndex:
    """
    鸽巢原理：长度为 n 的两条序列若 Hamming 距离 <= m，把序列切成 m+1 段，
    至少有一段完全相同。于是只需在 m+1 个“段 -> k-mer 列表”哈希表中查找候选，
    再对候选精确计算 Hamming 距离 —— 结果与 MATLAB 的暴力比对完全一致。
    """

    def __init__(self, kmers: np.ndarray, max_mm: int):
        self.kmers = kmers
        n = kmers.shape[1]
        self.bounds = np.linspace(0, n, max_mm + 2).astype(int)
        self.tables = []
        for a, b in zip(self.bounds[:-1], self.bounds[1:]):
            seg = np.ascontiguousarray(kmers[:, a:b]).view(np.dtype((np.void, b - a)))[:, 0]
            order = np.argsort(seg, kind="stable")
            seg_sorted = seg[order]
            uniq, start = np.unique(seg_sorted, return_index=True)
            stop = np.append(start[1:], len(seg_sorted))
            self.tables.append({u.tobytes(): order[s:e] for u, s, e in zip(uniq, start, stop)})

    def query(self, read: np.ndarray, max_mm: int):
        cand = [t.get(read[a:b].tobytes()) for t, a, b in
                zip(self.tables, self.bounds[:-1], self.bounds[1:])]
        cand = [c for c in cand if c is not None]
        if not cand:
            return np.zeros(0, np.int64), np.zeros(0, np.int64)
        cand = np.unique(np.concatenate(cand))
        d = (self.kmers[cand] != read).sum(1)
        ok = d <= max_mm
        return cand[ok], d[ok]


def build_A_matrices(dbs, regions, cfg: Config, stats: ReadsStats | None = None):
    """
    A_r[y, j] = Pr(read y | 细菌 j) = (1-pe)^(L-d) * (pe/3)^d ,  d = Hamming(read y, 细菌 j 在区域 r 的 k-mer) <= 2
    """
    nB = len(dbs[0]["indInSeqs"])
    dat0 = {"A": [], "F": [], "reads": []}
    for rr, (db, reg) in enumerate(zip(dbs, regions), start=1):
        t = time.time()
        # 5a. 用于本区域：先在区域内再去重（去掉引物后可能合并）
        reads, F = unique_rows_with_counts(reg["reads"], reg["freq"].astype(float))
        # 5b. 低丰度过滤: count < max(1e-4 * total, 2) 的 read 被丢弃
        thr = max(cfg.min_read_freq * F.sum(), cfg.min_read_count)
        keep = F >= thr
        reads, F = reads[keep], F[keep]
        if stats is not None:
            stats.add("After low abundance filter", len(F), int(F.sum()), rr)

        # 5c. k-mer -> 细菌 的倒排（每个细菌在每个区域只有一个 k-mer）
        idx = db["indInSeqs"]
        amp = np.flatnonzero(idx > 0)
        order = amp[np.argsort(idx[amp], kind="stable")]
        k_of = idx[order] - 1
        K = db["kmers"].shape[0]
        starts = np.searchsorted(k_of, np.arange(K))
        stops = np.searchsorted(k_of, np.arange(K), side="right")

        L = db["kmers"].shape[1]
        index = KmerIndex(db["kmers"], cfg.nMM_cut)
        rows, cols, vals = [], [], []
        for y in range(len(F)):
            ks, d = index.query(reads[y], cfg.nMM_cut)
            if len(ks) == 0:
                continue
            p = (1 - cfg.pe) ** (L - d) * (cfg.pe / 3) ** d
            # M 矩阵按列归一化（PE 下每列只有一个 1，归一化后为 1/(1+eps)）
            p = p * (1.0 / (1.0 + EPS))
            for k, pk in zip(ks, p):
                bact = order[starts[k]:stops[k]]
                rows.append(np.full(len(bact), y))
                cols.append(bact)
                vals.append(np.full(len(bact), pk))
        if rows:
            A = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                              shape=(len(F), nB))
        else:
            A = sp.csr_matrix((len(F), nB))
        dat0["A"].append(A)
        dat0["F"].append(F)
        dat0["reads"].append(reads)
        print(f"Region {rr}: {len(F):,} unique reads (>= {thr:.1f}), A nnz = {A.nnz:,} [{time.time() - t:.1f}s]")
    dat0["is_perfect_match"] = np.column_stack([db["is_perfect_match"] for db in dbs])
    dat0["indInSeqs"] = np.column_stack([db["indInSeqs"] for db in dbs])
    return dat0


# --------------------------------------------------------------------------
# 6. 过滤候选细菌 + 合并不可区分者 + EM 求解     (solve_iterative_noisy.m, ml_em_iterative.m)
# --------------------------------------------------------------------------
def ml_em_iterative(A, y, cfg: Config, x0=None, verbose=True):
    """
    最大似然 EM（与 Richardson-Lucy / MLEM 同形）:
        theta = A x ;  x <- x * A^T (y / theta)
    收敛判据: sum_j |1 - factor_j| * x_j < tol
    （MATLAB 版本每运行 60 秒会剪除 x<1e-10 的列，这依赖机器速度，本实现不做此剪枝）
    """
    A = sp.csc_matrix(A)
    AT = A.T.tocsr()
    x = AT @ y if x0 is None else x0.copy()
    x = x / x.sum()
    t0 = time.time()
    for it in range(1, cfg.num_iter + 1):
        theta = A @ x
        factor = AT @ (y / (theta + EPS))
        err = np.abs(1 - factor) @ x
        x = x * factor
        if err < cfg.tol:
            break
    if verbose:
        print(f"EM: {it} iterations, L1 error {err:.2e}, {time.time() - t0:.1f}s")
    x[x < 1e-10] = 0
    return x / x.sum(), it


def _column_keys(A):
    A = sp.csc_matrix(A)
    A.sort_indices()
    keys = []
    for j in range(A.shape[1]):
        s, e = A.indptr[j], A.indptr[j + 1]
        keys.append((A.indices[s:e].tobytes(), A.data[s:e].tobytes()))
    return keys


def solve_iterative_noisy(dat0, cfg: Config, stats: ReadsStats | None = None, return_debug=False):
    nR = len(dat0["A"])
    nB_all = dat0["A"][0].shape[1]
    sumA = np.zeros((nR, nB_all))
    n_pm = np.zeros((nR, nB_all))
    A_mat, F_vec = [], []
    for i in range(nR):
        A, F = dat0["A"][i], dat0["F"][i]
        rL = dat0["reads"][i].shape[1]
        total_pr_thresh = (1 - cfg.pe) ** (rL - cfg.nMM_cut) * (cfg.pe / 3) ** cfg.nMM_cut
        perfect_pr = (1 - cfg.pe) ** rL
        rs = np.asarray(A.sum(1)).ravel()
        keep = rs > 0
        if stats is not None:
            stats.add("Number of reads mathced to DB", int(keep.sum()), int(F[keep].sum()), i + 1)
            stats.add("Count of maximal unmatched read", 1, int(max([0, *F[~keep]])), i + 1)
            stats.add("Count of all unmatched read", int((~keep).sum()), int(F[~keep].sum()), i + 1)
        A_mat.append(A[keep])
        F_vec.append(F[keep])
        sumA[i] = np.asarray(A_mat[-1].sum(0)).ravel()
        n_pm[i] = np.asarray((A_mat[-1] > perfect_pr - 0.1 * total_pr_thresh - EPS).sum(0)).ravel()

    comb = dat0["indInSeqs"].T > 0                      # nR x nB : 区域是否被扩增
    norm_factor = comb.sum(0)                           # 每个细菌被扩增的区域数

    # (1) 候选过滤：凡是引物 0 错配（理论上必被扩增）的区域，必须至少有一条完全匹配的 read
    missing = np.zeros((nR, nB_all), bool)
    if cfg.do_filter == 1:
        pm = dat0["is_perfect_match"].T
        missing = (pm & (n_pm < 1)) | (~pm & (n_pm > 0) & (n_pm < 1))
    keep_col = np.flatnonzero((missing.sum(0) == 0) & (sumA.sum(0) > 0) & (norm_factor > 0))
    n_after_filter = len(keep_col)

    # (2) 各区域 read 计数统一按总数归一化 -> y
    total = sum(F.sum() for F in F_vec)
    y = np.concatenate([F / total for F in F_vec])
    A_L2 = sp.vstack([A[:, keep_col] for A in A_mat]).tocsc()
    A_L2 = A_L2 @ sp.diags(1.0 / norm_factor[keep_col])

    # (3) 列完全相同的细菌 = 在这 5 个区域上不可区分 -> 合并为一个 "group"
    keys = _column_keys(A_L2)
    first = {}
    groups = []
    rep = []
    for j, k in enumerate(keys):
        if k not in first:
            first[k] = len(rep)
            rep.append(j)
            groups.append([])
        groups[first[k]].append(keep_col[j])
    # MATLAB 的 unique(A','rows') 会按列内容排序；这里按代表列在 keep_col 中出现的顺序排列（只影响输出顺序）
    rep = np.array(rep)
    A_L2 = A_L2[:, rep]
    keep_col = keep_col[rep]
    n_groups = len(rep)

    # (4) 去除“被包含”的细菌：若 i 的非零 read 集合是 j 的真子集，则删除 i
    if cfg.filter_included_bacteria:
        TF = (A_L2 > 0).astype(np.float64).tocsc()
        nnz = np.asarray(TF.sum(0)).ravel()
        inter = (TF.T @ TF).tocsr()
        remove = np.zeros(len(nnz), bool)
        for i in range(len(nnz)):
            s, e = inter.indptr[i], inter.indptr[i + 1]
            jj, cnt = inter.indices[s:e], inter.data[s:e]
            if np.any((cnt == nnz[i]) & (nnz[jj] > nnz[i])):
                remove[i] = True
        A_L2 = A_L2[:, ~remove]
        keep_col = keep_col[~remove]
        groups = [g for g, r in zip(groups, remove) if not r]
        print(f"Removed {remove.sum()} included bacteria out of {len(remove)}")

    # (5) EM
    A_em = A_L2
    x, n_iter = ml_em_iterative(A_L2, y, cfg)
    use = np.flatnonzero(x > 1e-10)
    x = x[use]
    keep_col = keep_col[use]
    groups = [groups[u] for u in use]
    A_L2 = A_L2[:, use]

    # (6) 读数比例 -> 细胞比例：除以被扩增区域数
    freq = x / norm_factor[keep_col]
    freq = freq / freq.sum()
    meta = [{"db_ind": np.array(g), "comb_vec": comb[:, k]} for g, k in zip(groups, keep_col)]
    debug = dict(n_candidates=n_after_filter, n_groups=n_groups, n_final=len(x), n_iter=n_iter,
                 A_em=A_em, A_L2=A_L2, y=y)
    if return_debug:
        return freq, meta, keep_col, debug
    return freq, meta, keep_col


def assign_reads(dat0, freq, keep_col):
    """按后验 Pr(j | read) 把 read 计数分给各细菌   (reconstruction_func.m 末尾)"""
    counts = np.zeros((len(keep_col), len(dat0["A"])))
    for r, (A, F) in enumerate(zip(dat0["A"], dat0["F"])):
        P = A[:, keep_col].multiply(freq[None, :]).tocsr()
        rs = np.asarray(P.sum(1)).ravel() + EPS
        post = sp.diags(F / rs) @ P
        counts[:, r] = np.round(np.asarray(post.sum(0)).ravel())
    return counts.sum(1)


# --------------------------------------------------------------------------
# 7. 分类注释 + 汇总成表     (save_reconstruction_new_nogroups.m, scott_format_newer_func.m)
# --------------------------------------------------------------------------
RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]


def annotate_groups(freq, assigned_reads, meta, taxa: pd.DataFrame, level="species"):
    """
    每个 group 内可能包含多条参考序列（分类可能不同）。
    MATLAB 输出采用“硬决策”：整个 group 的频率赋给 group 内占比最高的分类
    （并列时取字典序最小者），reads 数按该分类的占比折算。
    """
    tl = RANKS.index(level) + 1
    rows = []
    for f, r, m in zip(freq, assigned_reads, meta):
        sub = taxa.iloc[np.sort(m["db_ind"]), :tl]
        vc = sub.value_counts(sort=False).sort_index()        # 字典序
        frac = vc.values / vc.values.sum()
        best = int(np.argmax(frac))
        rows.append(list(vc.index[best]) + [f, round(r * frac[best]), len(m["db_ind"]), frac[best]])
    df = pd.DataFrame(rows, columns=RANKS[:tl] + ["freq", "reads", "group_size", "majority_fraction"])
    for c in RANKS[:tl]:
        df[c] = df[c].str.replace("Assigned", "Unknown")
    return df


def merge_samples(per_sample: dict, level="species"):
    tl = RANKS.index(level) + 1
    cols = RANKS[:tl]
    tables, totals = [], {}
    for name, df in per_sample.items():
        tables.append(df.groupby(cols, sort=True)["freq"].sum().rename(name))
        totals[name] = int(df["reads"].sum())
    out = pd.concat(tables, axis=1).fillna(0.0).sort_index()
    return out, pd.Series(totals)


def write_matlab_style_table(table: pd.DataFrame, totals: pd.Series, path):
    """与 saveCellFile.m 相同的 tab 分隔格式"""
    nlev = table.index.nlevels
    with open(path, "w") as fh:
        fh.write("Total # of reads" + "\t" * nlev + "\t".join(str(totals[c]) for c in table.columns) + "\n")
        fh.write("\t".join(list(table.index.names) + list(table.columns)) + "\n")
        for idx, row in table.iterrows():
            fh.write("\t".join(list(idx) + [f"{v:f}" for v in row.values]) + "\n")


# --------------------------------------------------------------------------
# 一键运行  (main_5R.m / main_multiple_regions.m)
# --------------------------------------------------------------------------
def run_sample(r1, r2, dbs, cfg: Config, taxa=None):
    stats = ReadsStats(len(cfg.primers))
    Suni, freq = quality_filter_pairs(r1, r2, cfg, stats)
    regions = split_to_regions(Suni, freq, cfg, stats)
    dat0 = build_A_matrices(dbs, regions, cfg, stats)
    bfreq, meta, keep_col = solve_iterative_noisy(dat0, cfg, stats)
    reads = assign_reads(dat0, bfreq, keep_col)
    res = {"freq": bfreq, "meta": meta, "keep_col": keep_col, "assigned_reads": reads, "stats": stats}
    if taxa is not None:
        res["table"] = annotate_groups(bfreq, reads, meta, taxa)
    return res


# --------------------------------------------------------------------------
# 8. 扩展工具：引物覆盖度评估 & 从任意 16S FASTA 构建 5R k-mer 数据库
#    （原仓库没有提供数据库构建脚本；以下实现用于更换/更新数据库，如 SILVA / GTDB / Greengenes2）
# --------------------------------------------------------------------------
def primer_coverage(db_dir, taxa: pd.DataFrame, cfg: Config, rank="phylum", db_name="GreenGenes_201305", top=15):
    """统计每个区域在各分类群中的 (a) 可扩增比例(<=2 错配) 与 (b) 引物 0 错配比例"""
    rows = {}
    for rr in range(1, len(cfg.primers) + 1):
        path = os.path.join(db_dir, f"{cfg.db_file_prefix(db_name)}_region{rr}.mat")
        d = sio.loadmat(path, variable_names=["indInValue", "is_perfect_match"])
        rows[f"R{rr}_amp"] = d["indInValue"].ravel() > 0
        rows[f"R{rr}_0mm"] = d["is_perfect_match"].ravel().astype(bool) & rows[f"R{rr}_amp"]
    df = pd.DataFrame(rows)
    amp = df[[c for c in df if c.endswith("_amp")]]
    df["n_regions"] = amp.sum(1)
    df[rank] = taxa[rank].values
    top_taxa = df[rank].value_counts().index[:top]
    g = df[df[rank].isin(top_taxa)].groupby(rank)
    out = g.mean().loc[top_taxa]
    out.insert(0, "n_seqs", g.size().loc[top_taxa])
    overall = df.drop(columns=[rank]).mean()
    overall["n_seqs"] = len(df)
    out.loc["ALL"] = overall
    return out


def _best_hit(seq_u8, primer_variants_u8, max_mm, lo=0, hi=None):
    """在 seq[lo:hi] 中滑动比对（任一简并展开版本），返回 (位置, 错配数)；找不到返回 (-1, inf)"""
    hi = len(seq_u8) if hi is None else hi
    plen = primer_variants_u8.shape[1]
    if hi - lo < plen:
        return -1, np.inf
    win = np.lib.stride_tricks.sliding_window_view(seq_u8[lo:hi], plen)       # (n_pos, plen)
    mm = np.min([(win != p).sum(1) for p in primer_variants_u8], axis=0)
    pos = int(np.argmin(mm))
    return (lo + pos, int(mm[pos])) if mm[pos] <= max_mm else (-1, np.inf)


def build_region_db_from_fasta(seqs: dict, cfg: Config, rr: int, max_amplicon=600):
    """
    对一组全长 16S 序列做 in-silico PCR，生成与原 .mat 同结构的区域数据库：
      values (K x 2*db_kmer_len) / indInValue (1-based, 0=不扩增) / is_perfect_match
    seqs: {header: sequence}
    """
    fp, rp = cfg.primers[rr - 1]
    F = to_u8(expand_degenerate(fp))
    R = to_u8([revcomp(x) for x in expand_degenerate(rp)])
    k = cfg.db_kmer_len
    vals, ind, pm = {}, [], []
    for h, s in seqs.items():
        s = s.upper().replace("U", "T")
        u = np.frombuffer(s.encode(), np.uint8)
        fpos, fmm = _best_hit(u, F, cfg.allowed_mm)
        if fpos < 0:
            ind.append(0); pm.append(False); continue
        rpos, rmm = _best_hit(u, R, cfg.allowed_mm, fpos + F.shape[1], min(len(u), fpos + max_amplicon))
        if rpos < 0:
            ind.append(0); pm.append(False); continue
        amp = s[fpos:rpos + R.shape[1]]
        if len(amp) < k:
            ind.append(0); pm.append(False); continue
        v = amp[:k] + amp[-k:]
        ind.append(vals.setdefault(v, len(vals) + 1))
        pm.append(fmm == 0 and rmm == 0)
    values = np.array(list(vals.keys()), dtype=f"<U{2 * k}")
    return {"values": values, "indInValue": np.array(ind, float)[:, None],
            "is_perfect_match": np.array(pm)[:, None]}
