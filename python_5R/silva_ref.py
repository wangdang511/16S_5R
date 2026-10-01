"""
silva_ref.py -- 把 SILVA 128 (SEPP 参考包里的 99% OTU 全长比对) 投影到 E. coli 坐标，并用 Greengenes 训练的
8-mer 朴素贝叶斯分类器推断门/域，得到与 primer_design.Reference 相同结构的对象。
"""
import re, pickle, numpy as np
import primer_design as pdz


def sample_alignment(fa_path, n, seed=0):
    """蓄水池采样 n 条序列；返回 (ids, (n, ncol) uint8 比对矩阵)"""
    rng = np.random.default_rng(seed)
    keep_ids, rows = [], []
    seen = 0
    hdr, buf = None, []

    def flush():
        nonlocal seen
        if hdr is None:
            return
        seq = "".join(buf)
        if len(rows) < n:
            keep_ids.append(hdr); rows.append(np.frombuffer(seq.encode(), np.uint8))
        else:
            j = rng.integers(0, seen + 1)
            if j < n:
                keep_ids[j] = hdr; rows[j] = np.frombuffer(seq.encode(), np.uint8)
        seen += 1

    with open(fa_path) as fh:
        for line in fh:
            if line[0] == ">":
                flush(); hdr, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
        flush()
    return keep_ids, np.vstack(rows), seen


def ecoli_columns(A, ref_seq):
    """
    在抽样中找一条与 E. coli 参照高度相似的序列，用两两比对把 E. coli 位置 1..1542 映射到比对列。
    缺失位置用相邻位置插值。
    """
    from Bio import Align
    al = Align.PairwiseAligner(mode="global", match_score=2, mismatch_score=-3, open_gap_score=-6, extend_gap_score=-2)
    al.end_insertion_score = 0; al.end_deletion_score = 0
    best = None
    anchors = [b"AAACTCAAATGAATTGACGG", b"ACTCCTACGGGAGGCAG", b"GTGCCAGCAGCCGCGGTAA", b"GGAGGAAGGTGGGGATGAC"]
    for i in range(A.shape[0]):
        raw = A[i][(A[i] != ord("-")) & (A[i] != ord("."))].tobytes()
        if len(raw) < 1450 or not all(a in raw for a in anchors):
            continue
        s = raw.decode().upper().replace("U", "T")
        sc = al.score(ref_seq.replace("N", "A"), s) / len(s)
        if best is None or sc > best[0]:
            best = (sc, i, s)
        if best[0] > 1.9:
            break
    sc, i, s = best
    cols = np.flatnonzero((A[i] != ord("-")) & (A[i] != ord(".")))        # 序列第 k 个碱基所在的比对列
    aln = al.align(ref_seq.replace("N", "A"), s)[0]
    pos2col = np.full(len(ref_seq), -1)
    for (rs, re_), (qs, qe) in zip(*aln.aligned):
        pos2col[rs:re_] = cols[qs:qe]
    miss = np.flatnonzero(pos2col < 0)
    ok = np.flatnonzero(pos2col >= 0)
    pos2col[miss] = np.round(np.interp(miss, ok, pos2col[ok])).astype(int)
    return pos2col, best[0], i


def classify_phylum(seqs_ungapped, gg_dir, k=8, max_train=None):
    """Greengenes 97% 代表序列训练 8-mer 朴素贝叶斯 -> 预测 (domain, phylum)"""
    import gzip, scipy.sparse as sp
    from sklearn.naive_bayes import MultinomialNB
    tr = pdz._read_fasta_gz(f"{gg_dir}/rep_set/97_otus.fasta.gz")
    tax = {}
    with gzip.open(f"{gg_dir}/taxonomy/97_otu_taxonomy.txt.gz", "rt") as fh:
        for l in fh:
            a, b = l.rstrip("\n").split("\t"); tax[a] = [x.strip()[3:] for x in b.split(";")]
    ids = [i for i in tr if i in tax and tax[i][1]]
    if max_train:
        ids = ids[:max_train]
    lut = np.full(256, -1, np.int64)
    for j, c in enumerate("ACGT"):
        lut[ord(c)] = j

    def kmer_counts(seq):
        a = lut[np.frombuffer(seq.encode(), np.uint8)]
        n = len(a) - k + 1
        if n <= 0:
            return np.zeros(0, np.int64)
        idx = np.zeros(n, np.int64)
        ok = np.ones(n, bool)
        for t in range(k):
            x = a[t:t + n]; ok &= x >= 0; idx = idx * 4 + np.where(x >= 0, x, 0)
        return np.unique(idx[ok])

    def matrix(seqs):
        rows, cols = [], []
        for r, s in enumerate(seqs):
            c = kmer_counts(s); rows.append(np.full(len(c), r)); cols.append(c)
        r = np.concatenate(rows); c = np.concatenate(cols)
        return sp.csr_matrix((np.ones(len(r), np.float32), (r, c)), shape=(len(seqs), 4 ** k))

    Xtr = matrix([tr[i].upper() for i in ids])
    label = np.array([tax[i][0] + "|" + tax[i][1] for i in ids])
    nb = MultinomialNB(alpha=0.1).fit(Xtr, label)
    Xte = matrix(seqs_ungapped)
    lp = nb.predict_log_proba(Xte)
    pred = nb.classes_[lp.argmax(1)]
    # 置信度: 最大类的后验概率
    post = np.exp(lp.max(1))
    return np.array([p.split("|")[0] for p in pred]), np.array([p.split("|")[1] for p in pred]), post, (nb, matrix, ids, label)


def build_reference(fa_path, gg_dir, ref: "pdz.Reference", n=60000, seed=0, out_pkl=None):
    ids, A, total = sample_alignment(fa_path, n, seed)
    pos2col, score, gi = ecoli_columns(A, ref.ref_seq)
    M = A[:, pos2col].copy()
    M[(M == ord(".")) | (M == 0)] = ord("-")
    M = np.char.upper(M.view("S1")).view(np.uint8) if False else M
    low = (M >= ord("a")) & (M <= ord("z"))
    M[low] -= 32
    M[M == ord("U")] = ord("T")
    M[~np.isin(M, np.frombuffer(b"ACGT-", np.uint8))] = ord("N")
    # 序列两端的缺口 = 没有数据
    first = np.array([np.flatnonzero(A[r] != ord("-"))[0] if (A[r] != ord("-")).any() else A.shape[1] for r in range(len(A))])
    last = np.array([np.flatnonzero(A[r] != ord("-"))[-1] if (A[r] != ord("-")).any() else -1 for r in range(len(A))])
    for r in range(len(A)):
        M[r, pos2col < first[r]] = 0
        M[r, pos2col > last[r]] = 0
    ungapped = [A[r][A[r] != ord("-")].tobytes().decode().upper().replace("U", "T") for r in range(len(A))]
    dom, phy, post, _ = classify_phylum(ungapped, gg_dir)
    S = pdz.Reference(M, ids, dom, phy, {p: p - 1 for p in range(1, 1543)}, ref.ref_seq, {})
    S.meta = dict(total=total, guide=gi, guide_score=score, post=post)
    if out_pkl:
        pickle.dump(S, open(out_pkl, "wb"))
    return S
