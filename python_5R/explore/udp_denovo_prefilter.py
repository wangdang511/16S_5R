"""从头生成 10-nt UDP 候选并做预筛。
规则：GC 4–6/10、无同聚物≥3、无二核苷酸重复≥3 次、非近回文（与自身反向互补的编辑距离≥4）、与 Illumina UDP0001–0384 所有序列（及其反向互补）编辑距离≥3（不继承）。
预筛指标（58 °C，primer3）：尾巴 × 池里 198 个展开序列（不带尾巴）的最差异源二聚体 ΔG（H）、池里引物 3′端 对尾巴的最差 3′锚定 ΔG（E）。
用法：python3 udp_denovo_prefilter.py <订购清单 xlsx> <Illumina udp tsv> <候选数> <输出 pkl> <seed>"""
import sys, random, pickle, re
import numpy as np, openpyxl, primer3
from multiprocessing import Pool
import primer_design as pdz
xlsx, ill_tsv, NC, out, seed = sys.argv[1:6]; NC = int(NC); random.seed(int(seed))
KW = dict(mv_conc=50, dv_conc=2, dntp_conc=0.2, dna_conc=250); T = 58.0
ws = openpyxl.load_workbook(xlsx)['订购清单']
BODY = sorted({e for r in ws.iter_rows(min_row=4, max_row=36, values_only=True) for e in pdz.expand(r[4])})
rc = lambda s: s[::-1].translate(str.maketrans('ACGT', 'TGCA'))
def lev(a, b):
    p = list(range(len(b) + 1))
    for i, x in enumerate(a, 1):
        c = [i]
        for j, y in enumerate(b, 1): c.append(min(p[j] + 1, c[-1] + 1, p[j - 1] + (x != y)))
        p = c
    return p[-1]
ILL = set()
for l in open(ill_tsv):
    r = l.split(); ILL |= {r[1], r[3]}
ILLR = ILL | {rc(s) for s in ILL}
gc = lambda s: sum(c in 'GC' for c in s)
def ok(s):
    if not 4 <= gc(s) <= 6: return False
    if re.search(r'(.)\1\1', s): return False
    if re.search(r'(..)\1\1', s): return False
    if lev(s, rc(s)) < 4: return False
    return all(lev(s, x) >= 3 for x in ILLR)
def proxy(s):
    H = min(primer3.calc_heterodimer(s, b, temp_c=T, **KW).dg for b in BODY) / 1000
    E = min(primer3.calc_end_stability(b, s, temp_c=T, **KW).dg for b in BODY) / 1000
    return s, H, E
if __name__ == '__main__':
    cand = set()
    while len(cand) < NC:
        s = ''.join(random.choice('ACGT') for _ in range(10))
        if s not in cand and ok(s): cand.add(s)
    with Pool(4) as p: res = p.map(proxy, sorted(cand), chunksize=20)
    pickle.dump(res, open(out, 'wb')); print('完成', len(res))
