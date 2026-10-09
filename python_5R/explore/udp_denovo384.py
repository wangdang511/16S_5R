"""自建 UDP：生成 384 对（384 个 i5 + 384 个 i7，分成 8 个 48 组块，每块的 i5、i7 各自配色平衡）。
单条规则：10 nt，GC 4–6，无同聚物≥3，无二核苷酸重复≥3 次，非回文（与自身反向互补编辑距离≥4），无 GGGG，前两位不都是 G，不含接头 6-mer，
          不与 Illumina UDP0001–0384 的任一序列（含反向互补）相同。
集合规则：768 条序列两两 Levenshtein ≥4，且 x 与 rc(y) 的编辑距离 ≥4；
块规则：每块 48 个 i5 和 48 个 i7 各自每个周期 GC、A+C、A+T 占比 40–60% 且 A/C/G/T 各 ≥15%（8 块的任意并集也满足）。
用法：python3 udp_denovo384.py <订购清单 xlsx> <Illumina tsv> <输出目录> <DMIN> <RCMIN> <Hmin> <Emin> <seed>"""
import sys, itertools, random, pickle, json, math, re
import numpy as np, openpyxl, primer3
from multiprocessing import Pool
from rapidfuzz.process import cdist
from rapidfuzz.distance import Levenshtein
import primer_design as pdz
xlsx, ill_tsv, outdir, DMIN, RCMIN, HMIN, EMIN, seed = sys.argv[1:9]; DMIN, RCMIN, HMIN, EMIN, seed = int(DMIN), int(RCMIN), float(HMIN), float(EMIN), int(seed)
rc = lambda s: s[::-1].translate(str.maketrans('ACGT', 'TGCA'))
ADAPT = ['AGATCGGAAGAGC', 'CTGTCTCTTATACACATCT', 'AGATGTGTATAAGAGACAG', 'AATGATACGGCGACCACCGAGATCTACAC', 'CAAGCAGAAGACGGCATACGAGAT', 'TCGTCGGCAGCGTC', 'GTCTCGTGGGCTCGG']
K6 = set()
for a in ADAPT:
    for t in (a, rc(a)): K6 |= {t[i:i + 6] for i in range(len(t) - 5)}
ILL = set()
for l in open(ill_tsv):
    r = l.split(); ILL |= {r[1], r[3]}
ILLR = ILL | {rc(s) for s in ILL}
def single_ok(s):
    g = sum(c in 'GC' for c in s)
    if not 4 <= g <= 6: return False
    if re.search(r'(.)\1\1', s) or re.search(r'(..)\1\1', s) or 'GGGG' in s or (s[0] == 'G' and s[1] == 'G'): return False
    if any(s[i:i + 6] in K6 for i in range(5)): return False
    return s not in ILLR
if __name__ == '__main__':
    rnd = random.Random(seed)
    allc = [''.join(p) for p in itertools.product('ACGT', repeat=10)]
    allc = [s for s in allc if single_ok(s)]
    rcs = [rc(s) for s in allc]
    pal = np.array([Levenshtein.distance(a, b) for a, b in zip(allc, rcs)]); allc = [s for s, d in zip(allc, pal) if d >= 4]
    print('单条规则通过', len(allc), flush=True)
    import os
    if os.environ.get('ORDER') != 'lex': rnd.shuffle(allc)
    S = []; Src = []
    BS = 3000
    for b0 in range(0, len(allc), BS):
        B = allc[b0:b0 + BS]
        if S:
            C1 = cdist(B, S, scorer=Levenshtein.distance, dtype=np.uint8, workers=4).min(1) < DMIN
            C2 = cdist(B, Src, scorer=Levenshtein.distance, dtype=np.uint8, workers=4).min(1) < RCMIN
            B = [x for x, c1, c2 in zip(B, C1, C2) if not (c1 or c2)]
        acc = []
        for x in B:
            if all(Levenshtein.distance(x, y) >= DMIN and Levenshtein.distance(x, rc(y)) >= RCMIN for y in acc): acc.append(x)
        S += acc; Src += [rc(y) for y in acc]
        if (b0 // BS) % 20 == 0: print(b0, len(S), flush=True)
    print('独立集大小', len(S), flush=True)
    pickle.dump(S, open(outdir + '/indep.pkl', 'wb'))
