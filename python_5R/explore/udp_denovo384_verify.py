import sys, json, itertools, collections, re
import numpy as np
from rapidfuzz.process import cdist
from rapidfuzz.distance import Levenshtein
sel_json, ill_tsv, out = sys.argv[1:4]
d = json.load(open(sel_json)); sl = d['slots']; assert len(sl) == 768
rc = lambda s: s[::-1].translate(str.maketrans('ACGT', 'TGCA'))
i5 = [sl[g * 48 + k] for g in range(8) for k in range(48)]; i7 = [sl[(8 + g) * 48 + k] for g in range(8) for k in range(48)]
allseq = i5 + i7; res = {}
D = cdist(allseq, allseq, scorer=Levenshtein.distance, dtype=np.uint8, workers=4); np.fill_diagonal(D, 99)
Dr = cdist(allseq, [rc(s) for s in allseq], scorer=Levenshtein.distance, dtype=np.uint8, workers=4)
res['768 条互不相同'] = len(set(allseq)) == 768
res['全部序列两两最小编辑距离'] = int(D.min()); res['i5 之间最小'] = int(D[:384, :384].min()); res['i7 之间最小'] = int(D[384:, 384:].min()); res['i5 与 i7 之间最小'] = int(D[:384, 384:].min())
Dr2 = Dr.copy(); res['与反向互补最小编辑距离（含自身）'] = int(Dr2.min()); np.fill_diagonal(Dr2, 99); res['x 与 rc(y)（y≠x）最小'] = int(Dr2.min())
gcs = [sum(c in 'GC' for c in s) for s in allseq]; res['单个 GC 范围'] = (min(gcs), max(gcs))
res['同聚物≥3'] = sum(bool(re.search(r'(.)\1\1', s)) for s in allseq); res['含 GGGG'] = sum('GGGG' in s for s in allseq); res['前两位 GG'] = sum(s[:2] == 'GG' for s in allseq)
ill = set()
for l in open(ill_tsv):
    r = l.split(); ill |= {r[1], r[3]}
ham = lambda a, b: sum(x != y for x, y in zip(a, b))
res['与 Illumina 384 对：相同序列数'] = sum(s in ill or rc(s) in ill for s in allseq)
res['与 Illumina 的最小汉明距离'] = min(ham(s, t) for s in allseq for t in ill)
def cyc(seqs, N):
    ok = True; worst = []
    for p in range(10):
        c = collections.Counter(s[p] for s in seqs); gc = c['G'] + c['C']; red = c['A'] + c['C']; green = c['A'] + c['T']
        for v in (gc, red, green): ok &= (0.4 * N <= v <= 0.6 * N)
        ok &= all(c[b] >= 0.15 * N for b in 'ACGT')
        worst.append((min(c.values()), gc, red, green))
    return ok, worst
groups = {}
for name, seqs, N in [(f'i5 48组 {g + 1}', i5[g * 48:(g + 1) * 48], 48) for g in range(8)] + [(f'i7 48组 {g + 1}', i7[g * 48:(g + 1) * 48], 48) for g in range(8)] + \
        [(f'i5 96块 {b + 1}', i5[b * 96:(b + 1) * 96], 96) for b in range(4)] + [(f'i7 96块 {b + 1}', i7[b * 96:(b + 1) * 96], 96) for b in range(4)] + [('i5 全部384', i5, 384), ('i7 全部384', i7, 384)]:
    ok, w = cyc(seqs, N); groups[name] = ok
res['每周期达标的组'] = f'{sum(groups.values())}/{len(groups)}'; res['未达标'] = [k for k, v in groups.items() if not v]
print(json.dumps(res, ensure_ascii=False, indent=1))
json.dump(dict(res=res, i5=i5, i7=i7, H=d['H'], E=d['E']), open(out, 'w'), ensure_ascii=False)
