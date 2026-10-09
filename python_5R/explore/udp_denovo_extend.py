"""在已有 768 条自建序列（384 对）的基础上，找可以追加的序列：与已有序列及彼此都满足 Levenshtein ≥DMIN、与反向互补 ≥RCMIN。
用法：python3 udp_denovo_extend.py <final json> <Illumina tsv> <DMIN> <RCMIN> <SYM> <输出 pkl>"""
import sys, itertools, pickle, re, json
import numpy as np
from rapidfuzz.process import cdist
from rapidfuzz.distance import Levenshtein
sys.argv_ = sys.argv
fin, ill_tsv, DMIN, RCMIN, SYM, out = sys.argv[1:7]; DMIN, RCMIN = int(DMIN), int(RCMIN); TYPE = sys.argv[7] if len(sys.argv) > 7 else 'all'
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
d = json.load(open(fin)); both = set(d['i5'] + d['i7']); old = d['i5'] + d['i7'] if TYPE == 'all' else d[TYPE]
allc = [''.join(p) for p in itertools.product('ACGT', repeat=10)]
allc = [s for s in allc if single_ok(s) and s not in both and rc(s) not in both]
pal = np.array([Levenshtein.distance(a, rc(a)) for a in allc]); allc = [s for s, x in zip(allc, pal) if x >= 4]
oldr = [rc(s) for s in old]; keep = []
for b0 in range(0, len(allc), 20000):
    B = allc[b0:b0 + 20000]
    c1 = cdist(B, old, scorer=Levenshtein.distance, dtype=np.uint8, workers=2).min(1) < DMIN
    c2 = cdist(B, oldr, scorer=Levenshtein.distance, dtype=np.uint8, workers=2).min(1) < RCMIN
    keep += [x for x, a, b in zip(B, c1, c2) if not (a or b)]
print('与已有 768 条兼容的候选', len(keep), flush=True)
tb = str.maketrans(SYM, 'ACGT'); keep.sort(key=lambda x: x.translate(tb))
S = []; Src = []
for b0 in range(0, len(keep), 3000):
    B = keep[b0:b0 + 3000]
    if S:
        c1 = cdist(B, S, scorer=Levenshtein.distance, dtype=np.uint8, workers=2).min(1) < DMIN
        c2 = cdist(B, Src, scorer=Levenshtein.distance, dtype=np.uint8, workers=2).min(1) < RCMIN
        B = [x for x, a, b in zip(B, c1, c2) if not (a or b)]
    acc = []
    for x in B:
        if all(Levenshtein.distance(x, y) >= DMIN and Levenshtein.distance(x, rc(y)) >= RCMIN for y in acc): acc.append(x)
    S += acc; Src += [rc(y) for y in acc]
print('可追加的独立集大小', len(S), '（需要 192）', flush=True)
pickle.dump(S, open(out, 'wb'))
