"""第一阶段：对 960 条尾巴（U001–U960）各自、分别作为 i5（加在正向引物）和 i7（加在反向引物）评估。
只算“热点引物对”（见 udp960_hotpairs.py，38 个名称对，覆盖 99.6% 的整管分级）：
  F 角色：正向×正向（两个都带尾巴）FF、带尾巴正向×不带尾巴反向 FRu、带尾巴正向的发夹；
  R 角色：反向×反向 RR、带尾巴反向×不带尾巴正向 RFu、带尾巴反向的发夹。
每个指标 = (全局异源二聚体最差 ΔG, 3′端锚定最差 ΔG)；58 °C，50 mM Na⁺ / 2 mM Mg²⁺ / 0.2 mM dNTP / 250 nM。可断点续算。
用法：python3 udp960_stage1.py <订购清单 xlsx> <tails tsv> <hotpairs json> <输出 pkl>"""
import sys, os, json, pickle
import numpy as np, openpyxl, primer3
from multiprocessing import Pool
import primer_design as pdz
xlsx, tsv, hotj, out = sys.argv[1:5]
KW = dict(mv_conc=50, dv_conc=2, dntp_conc=0.2, dna_conc=250); T = 58.0
ws = openpyxl.load_workbook(xlsx)['订购清单']
EXP = {r[1]: pdz.expand(r[4]) for r in ws.iter_rows(min_row=4, max_row=36, values_only=True)}
ORI = {r[1]: r[2][-1] for r in ws.iter_rows(min_row=4, max_row=36, values_only=True)}
HOT = [tuple(x) for x in json.load(open(hotj))['hot']]
TAILS = [l.split()[:2] for l in open(tsv)]
def H(a, b): return primer3.calc_heterodimer(a, b, temp_c=T, **KW).dg / 1000
def E(a, b): return min(primer3.calc_end_stability(a, b, temp_c=T, **KW).dg, primer3.calc_end_stability(b, a, temp_c=T, **KW).dg) / 1000
def pm(A, B, same):
    h = e = 0.0
    for i, a in enumerate(A):
        for j, b in enumerate(B):
            if same and j < i: continue
            h = min(h, H(a, b)); e = min(e, E(a, b))
    return h, e
def do(tail):
    name, t = tail; res = dict(name=name, seq=t)
    for role, cl in (('F', 'FF'), ('R', 'RR')):
        res[cl] = (0.0, 0.0); res['FRu' if role == 'F' else 'RFu'] = (0.0, 0.0)
    for a, b in HOT:
        oa, ob = ORI[a], ORI[b]
        if oa == ob:
            cl = 'FF' if oa == 'F' else 'RR'; A = [t + e for e in EXP[a]]; B = [t + e for e in EXP[b]]
            h, e = pm(A, B, a == b); res[cl] = (min(res[cl][0], h), min(res[cl][1], e))
        else:
            f, r = (a, b) if oa == 'F' else (b, a)
            h, e = pm([t + x for x in EXP[f]], EXP[r], False); res['FRu'] = (min(res['FRu'][0], h), min(res['FRu'][1], e))
            h, e = pm([t + x for x in EXP[r]], EXP[f], False); res['RFu'] = (min(res['RFu'][0], h), min(res['RFu'][1], e))
    res['hpF'] = min(primer3.calc_hairpin(t + e, temp_c=T, **KW).dg for n in EXP if ORI[n] == 'F' for e in EXP[n]) / 1000
    res['hpR'] = min(primer3.calc_hairpin(t + e, temp_c=T, **KW).dg for n in EXP if ORI[n] == 'R' for e in EXP[n]) / 1000
    return res
if __name__ == '__main__':
    done = pickle.load(open(out, 'rb')) if os.path.exists(out) else {}
    todo = [t for t in TAILS if t[0] not in done]; print('待算', len(todo), flush=True)
    with Pool(4) as p:
        for k, r in enumerate(p.imap_unordered(do, todo, chunksize=2), 1):
            done[r['name']] = r
            if k % 20 == 0: pickle.dump(done, open(out, 'wb')); print('完成', len(done), flush=True)
    pickle.dump(done, open(out, 'wb')); print('全部完成', len(done), flush=True)
