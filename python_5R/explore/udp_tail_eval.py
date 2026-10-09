"""把 Illumina UDP 10 nt 内联尾巴（i5 加在正向引物、i7 加在反向引物 5′端）加到 33 条池的 198 个展开序列上，
逐个管（UDP0001..UDP0096，同一编号的 i5/i7 在同一管）评估二聚体、发夹和 Olivar 坏度，并与无尾巴的基线比较。
用法：PYTHONPATH=python_5R python3 python_5R/explore/udp_tail_eval.py <订购清单 xlsx> <udp tsv> <输出 pkl> [ANNEAL_C]"""
import sys, pickle, json, itertools, collections
import numpy as np, openpyxl, primer3
from multiprocessing import Pool
import primer_design as pdz
sys.path.insert(0, '/home/user/wangdang511/olivar_primer/src/olivar')
import design as od
xlsx, tsv, out = sys.argv[1:4]; T = float(sys.argv[4]) if len(sys.argv) > 4 else 58.0
KW = dict(mv_conc=50, dv_conc=2, dntp_conc=0.2, dna_conc=250)
ws = openpyxl.load_workbook(xlsx)['订购清单']
OL = [(r[1], r[2], r[4], r[7], r[11]) for r in ws.iter_rows(min_row=4, max_row=36, values_only=True)]
V = {}
for site in dict.fromkeys(o[1] for o in OL):
    g = [o for o in OL if o[1] == site]; sup = [o for o in g if o[4] == '补']; main = [o for o in g if o[4] != '补']
    for o in main: V[o[0]] = 2.5 * (0.9 if sup else 1) * o[3] / sum(x[3] for x in main)
    for o in sup: V[o[0]] = 2.5 * 0.1 * o[3] / sum(x[3] for x in sup)
EX = [(n, s[-1], e, V[n] / len(pdz.expand(q))) for n, s, q, ne, grp in OL for e in pdz.expand(q)]
UDP = [l.split() for l in open(tsv)]
N = len(EX)

def mats(seqs):
    H = np.zeros((N, N), np.float32); E = np.zeros((N, N), np.float32)
    for i in range(N):
        for j in range(i, N):
            H[i, j] = H[j, i] = primer3.calc_heterodimer(seqs[i], seqs[j], temp_c=T, **KW).dg / 1000
            E[i, j] = primer3.calc_end_stability(seqs[i], seqs[j], temp_c=T, **KW).dg / 1000
            E[j, i] = primer3.calc_end_stability(seqs[j], seqs[i], temp_c=T, **KW).dg / 1000 if j != i else E[i, j]
    return H, E

def hairs(seqs):
    r = [primer3.calc_hairpin(s, temp_c=T, **KW) for s in seqs]
    return np.array([x.dg / 1000 for x in r], np.float32), np.array([x.tm for x in r], np.float32)

def badness(seqs):
    F = [s.lower() for s, x in zip(seqs, EX) if x[1] == 'F']; R = [s.lower() for s, x in zip(seqs, EX) if x[1] == 'R']
    cf = [x[3] for x in EX if x[1] == 'F']; cr = [x[3] for x in EX if x[1] == 'R']; m = np.mean(cf + cr)
    return od.PrimerSetBadnessFast(F, R, fP_conc=[c / m for c in cf], rP_conc=[c / m for c in cr])[0]

def tube(k):
    name, i7, _, i5, _ = UDP[k]
    seqs = [(i5 if x[1] == 'F' else i7) + x[2] for x in EX]
    H, E = mats(seqs); hg, ht = hairs(seqs)
    tt = (primer3.calc_heterodimer(i5, i7, temp_c=T, **KW).dg / 1000, primer3.calc_end_stability(i5, i7, temp_c=T, **KW).dg / 1000,
          primer3.calc_end_stability(i7, i5, temp_c=T, **KW).dg / 1000)
    return k, dict(name=name, i5=i5, i7=i7, H=H, E=E, hg=hg, ht=ht, bad=badness(seqs), tail_tail=tt,
                   tail_selfhp=(primer3.calc_hairpin(i5, temp_c=T, **KW).dg / 1000, primer3.calc_hairpin(i7, temp_c=T, **KW).dg / 1000))

if __name__ == '__main__':
    base_seqs = [x[2] for x in EX]
    Hb, Eb = mats(base_seqs); hgb, htb = hairs(base_seqs)
    res = dict(EX=EX, base=dict(H=Hb, E=Eb, hg=hgb, ht=htb, bad=badness(base_seqs)), T=T, tubes={})
    print('基线完成', res['base']['bad'], flush=True)
    with Pool(4) as p:
        for k, d in p.imap_unordered(tube, range(len(UDP))):
            res['tubes'][k] = d; print('完成', k + 1, d['name'], flush=True)
            if len(res["tubes"]) % 4 == 0: pickle.dump(res, open(out, 'wb'))
    pickle.dump(res, open(out, 'wb')); print('全部完成', flush=True)
