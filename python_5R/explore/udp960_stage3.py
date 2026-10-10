"""第三阶段：候选池里每个 i5 × 每个 i7 的“正向×反向”热点相互作用（两个尾巴都加上）。
只算热点引物对里 正向×反向 的展开序列对；热点按重要性排序，一旦已经比 B 级还差（全局 < −8.5 或 3′端 < −5.5）就提前终止（该组合记为 C，并标记 early）。
输出：H[i,j]、E[i,j]、early[i,j]（i 为 i5 序号，j 为 i7 序号）。可断点续算。
用法：python3 udp960_stage3.py <订购清单 xlsx> <hotpairs json> <pool json> <tails tsv> <输出 pkl>"""
import sys, os, json, pickle
import numpy as np, openpyxl, primer3
from multiprocessing import Pool
import primer_design as pdz
xlsx, hotj, poolj, tsv, out = sys.argv[1:6]
KW = dict(mv_conc=50, dv_conc=2, dntp_conc=0.2, dna_conc=250); T = 58.0
ws = openpyxl.load_workbook(xlsx)['订购清单']
EXP = {r[1]: pdz.expand(r[4]) for r in ws.iter_rows(min_row=4, max_row=36, values_only=True)}
ORI = {r[1]: r[2][-1] for r in ws.iter_rows(min_row=4, max_row=36, values_only=True)}
HOT = [tuple(x) for x in json.load(open(hotj))['hot']]
FR = [((a, b) if ORI[a] == 'F' else (b, a)) for a, b in HOT if ORI[a] != ORI[b]]
pool = json.load(open(poolj)); seqs = {l.split()[0]: l.split()[1] for l in open(tsv)}
T5 = [seqs[n] for n in pool['i5']]; T7 = [seqs[n] for n in pool['i7']]
def H(a, b): return primer3.calc_heterodimer(a, b, temp_c=T, **KW).dg / 1000
def E(a, b): return min(primer3.calc_end_stability(a, b, temp_c=T, **KW).dg, primer3.calc_end_stability(b, a, temp_c=T, **KW).dg) / 1000
def row(i):
    t5 = T5[i]; Fl = {f: [t5 + e for e in EXP[f]] for f, r in FR}; hs = np.zeros(len(T7)); es = np.zeros(len(T7)); early = np.zeros(len(T7), bool)
    for j, t7 in enumerate(T7):
        h = e = 0.0
        for f, r in FR:
            Rl = [t7 + x for x in EXP[r]]
            for a in Fl[f]:
                for b in Rl:
                    h = min(h, H(a, b)); e = min(e, E(a, b))
            if h < -8.5 or e < -5.5: early[j] = True; break
        hs[j], es[j] = h, e
    return i, hs, es, early
if __name__ == '__main__':
    st = pickle.load(open(out, 'rb')) if os.path.exists(out) else dict(rows={})
    todo = [i for i in range(len(T5)) if i not in st['rows']]; print('待算行', len(todo), flush=True)
    with Pool(4) as p:
        for k, (i, hs, es, early) in enumerate(p.imap_unordered(row, todo), 1):
            st['rows'][i] = (hs, es, early)
            pickle.dump(st, open(out, "wb")); print("完成行", len(st["rows"]), flush=True)
    st['i5'] = pool['i5']; st['i7'] = pool['i7']; pickle.dump(st, open(out, 'wb')); print('全部完成', flush=True)
