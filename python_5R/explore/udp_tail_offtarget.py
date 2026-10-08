"""UDP 尾巴对脱靶的影响：对基因特异部分已有的脱靶命中（≤5 错配），取其 5′ 侧 10 nt 的基因组侧翼，
与每个 UDP 尾巴（正向引物 i5、反向引物 i7）比较，得到带尾巴的全长错配数。
用法：PYTHONPATH=python_5R python3 python_5R/explore/udp_tail_offtarget.py <human|pig> <基因组路径> <hits 合并 pkl> <udp tsv> <输出 pkl>"""
import sys, pickle, collections
import numpy as np
import primer_design as pdz
from offtarget import read_2bit, read_fa, enc
sys.path.insert(0, '/home/user/16S_5R/python_5R/explore')
which, gpath, hitpkl, tsv, out = sys.argv[1:6]
hits = pickle.load(open(hitpkl, 'rb'))[which]         # (chrom,pos,strand,label,mm,n)
UDP = [l.split() for l in open(tsv)]
byc = collections.defaultdict(list)
for h in hits: byc[h[0]].append(h)
if which == 'human':
    gen = ((n, c) for n, c in read_2bit(gpath))
else:
    import glob
    from offtarget_fasta import read_fasta_records
    def g():
        for f in sorted(glob.glob(gpath + '/*.fa')):
            yield from read_fasta_records(f)
    gen = g()
rows = []
for name, code in gen:
    if name not in byc: continue
    for c, pos, st, lab, mm, n in byc[name]:
        ori = 'F' if '-F.' in lab else 'R'
        if st == '+': a, b = pos - 10, pos
        else: a, b = pos + n, pos + n + 10
        if a < 0 or b > len(code): continue
        rows.append((ori, st, mm, code[a:b].copy(), lab, name, pos))
    print(name, len(rows), flush=True)
nh = len(rows); TM = np.zeros((nh, len(UDP)), np.int8)
for k, (nm, i7, _, i5, _) in enumerate(UDP):
    for ori, tail in (('F', i5), ('R', i7)):
        idx = [i for i, r in enumerate(rows) if r[0] == ori]
        for st in '+-':
            ii = [i for i in idx if rows[i][1] == st]
            if not ii: continue
            e = enc(tail if st == '+' else pdz.revcomp(tail))
            TM[ii, k] = np.array([(rows[i][3] != e).sum() for i in ii])
pickle.dump(dict(rows=[(r[0], r[1], r[2], r[4], r[5], r[6]) for r in rows], TM=TM), open(out, 'wb')); print('完成', nh)
