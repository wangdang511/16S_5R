"""
tiling_search.py -- 以“属水平最近邻准确率”为目标，在单管约束下搜索扩增子平铺

关键性质：扩增子长度 ≤ 2×读长−引物 时，内部碱基全部被测到，所以一个扩增子的匹配数矩阵 =
C[内部终点] − C[内部起点−1]（C 为逐位点累计的“两两匹配数”矩阵），任意平铺的距离矩阵是若干个差的和，
评估一次只需几毫秒，才能做模拟退火。
"""
import numpy as np, pandas as pd

CODE = {ord('A'): 0, ord('C'): 1, ord('G'): 2, ord('T'): 3, ord('-'): 4, 0: 5, ord('N'): 4}


def encode(aln):
    x = np.full(aln.shape, 4, np.uint8)
    for k, v in CODE.items():
        x[aln == k] = v
    return x


def cumulative_matches(x):
    """C[k][i,j] = 序列 i、j 在位置 1..k 上相同的个数（int16），形状 (L+1, n, n)"""
    n, L = x.shape
    C = np.zeros((L + 1, n, n), np.int16)
    for p in range(L):
        C[p + 1] = C[p] + (x[:, p, None] == x[None, :, p])
    return C


def nn_accuracy(mism, labels, valid=None):
    """mism: (n,n) 错配数；留一最近邻；并列按比例计分"""
    n = len(labels)
    d = mism.astype(np.int32).copy()
    np.fill_diagonal(d, 10 ** 6)
    m = d.min(1, keepdims=True)
    tied = d == m
    same = labels[None, :] == labels[:, None]
    credit = (tied & same).sum(1) / tied.sum(1)
    if valid is not None:
        credit = credit[valid]
    return float(credit.mean())


def accuracy_direct(x, intervals, labels, valid=None, chunk=2500):
    """用 one-hot 矩阵乘法直接算一个平铺的属准确率（不用累计矩阵，适合大验证集）"""
    pos = np.concatenate([np.arange(a - 1, b) for a, b in intervals]) if intervals else np.array([], int)
    n = len(x)
    S = np.zeros((n, n), np.float32)
    for s in range(0, len(pos), 400):
        P = pos[s:s + 400]
        X = np.zeros((n, len(P) * 6), np.float32)
        for k, p in enumerate(P):
            X[np.arange(n), k * 6 + x[:, p].astype(np.int64)] = 1
        S += X @ X.T
    mism = len(pos) - S
    return nn_accuracy(mism, labels, valid)


class NN:
    """预计算同属矩阵，加速重复评估"""
    def __init__(self, labels, valid=None):
        self.same = labels[None, :] == labels[:, None]
        np.fill_diagonal(self.same, False)
        self.n = len(labels)
        self.valid = valid if valid is not None else np.ones(len(labels), bool)
        self.diag = np.eye(len(labels), dtype=bool)

    def acc(self, mism):
        d = mism.astype(np.int16)
        d[self.diag] = 30000
        m = d.min(1)
        tied = d == m[:, None]
        credit = (tied & self.same).sum(1) / tied.sum(1)
        return float(credit[self.valid].mean())
