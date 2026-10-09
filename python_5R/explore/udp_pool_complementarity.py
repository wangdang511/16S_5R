"""混合建库后（连接接头、PCR 3–10 个循环）不同样本 UDP 之间是否会互补。
链的末端：正链 5′端=i5，3′端=rc(i7)；负链 5′端=i7，3′端=rc(i5)。
 - 头尾配对（会形成串联/锅柄结构）：i5_a×i7_b、i7_a×i5_b，含 a=b（锅柄，同一分子两端自己配对）。需要 i5_a 与 i7_b 序列相同或部分相同才会互补。
 - 同位置异源双链（i5_a×i5_b、i7_a×i7_b）：不同样本扩增子退火时本来就会形成，UDP 区有错配，错配越多越不稳定。
指标：最长精确互补段（nt）、10-mer 对 10-mer 的 primer3 异源二聚体 ΔG（60 °C 退火、72 °C 延伸）。
用法：python3 udp_pool_complementarity.py <选择 json> <48|96> <udp tsv> <输出 json> [选择 json 的 key，默认 res]"""
import sys, json, itertools, collections
import primer3
sel_json, N, tsv, out = sys.argv[1:5]
d = json.load(open(sel_json)); names = d['res'][N]['sel']
U = {l.split()[0]: l.split() for l in open(tsv)}
S = [(n, U[n][3], U[n][1]) for n in names]           # name, i5, i7
rc = lambda s: s[::-1].translate(str.maketrans('ACGT', 'TGCA'))
KW = dict(mv_conc=50, dv_conc=2, dntp_conc=0.2, dna_conc=250)
def lcs_comp(s, t):   # s 与 t 之间最长精确互补段 = s 与 rc(t) 的最长公共子串
    r = rc(t); best = 0
    for i in range(len(s)):
        for j in range(len(r)):
            k = 0
            while i + k < len(s) and j + k < len(r) and s[i + k] == r[j + k]: k += 1
            best = max(best, k)
    return best
def dg(s, t, T): return primer3.calc_heterodimer(s, t, temp_c=T, **KW).dg / 1000
res = {}
def scan(label, pairs):
    rows = []
    for (na, s), (nb, t) in pairs:
        rows.append((lcs_comp(s, t), dg(s, t, 60), dg(s, t, 72), na, nb, s, t))
    l = [r[0] for r in rows]; worst60 = min(rows, key=lambda r: r[1]); lmax = max(rows, key=lambda r: (r[0], -r[1]))
    res[label] = dict(n=len(rows), lcs_hist=dict(sorted(collections.Counter(l).items())), worst_dg60=worst60[1], worst_dg60_pair=worst60[3:],
                      worst_dg72=min(r[2] for r in rows), max_lcs=lmax[0], max_lcs_pair=lmax[3:], n_dg60_le6=sum(r[1] <= -6 for r in rows),
                      n_dg60_le4=sum(r[1] <= -4 for r in rows), n_dg60_le2=sum(r[1] <= -2 for r in rows))
    print(label, res[label]['n'], '最长互补段', res[label]['max_lcs'], res[label]['max_lcs_pair'], '60°C 最差 ΔG %.1f' % res[label]['worst_dg60'], 'ΔG≤-4: %d ≤-2: %d' % (res[label]['n_dg60_le4'], res[label]['n_dg60_le2']), flush=True)
# 头尾：top_a 5′ 端 i5_a（5′→3′）与 top_b 3′ 端 rc(i7_b)
scan('头尾 i5_a×i7_b（含 a=b 锅柄）', [((a[0] + ':i5', a[1]), (b[0] + ':i7', rc(b[2]))) for a in S for b in S])
scan('锅柄（仅 a=b）i5_a×i7_a', [((a[0] + ':i5', a[1]), (a[0] + ':i7', rc(a[2]))) for a in S])
# 同位置异源双链：top_a 5′ i5_a × bottom_b 3′ rc(i5_b)，a≠b
scan('异源双链 i5_a×i5_b', [((a[0] + ':i5', a[1]), (b[0] + ':i5', rc(b[1]))) for a, b in itertools.permutations(S, 2)])
scan('异源双链 i7_a×i7_b', [((a[0] + ':i7', a[2]), (b[0] + ':i7', rc(b[2]))) for a, b in itertools.permutations(S, 2)])
json.dump(res, open(out, 'w'), ensure_ascii=False, indent=1, default=str)
