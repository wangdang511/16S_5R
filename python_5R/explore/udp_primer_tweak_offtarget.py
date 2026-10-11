"""改动引物后的脱靶：基因特异部分已有命中（≤5 错配）里，逐位比较，重新数改动后 ≤2 错配的条目数。
条目 = (染色体, 位置, 链, 展开序列)，与订购清单“≤2错配位点”同一口径（原始 F.2 人 17、猪 12；R.1 人 28、猪 11）。
用法：PYTHONPATH=python_5R python3 udp_primer_tweak_offtarget.py <human|pig> <基因组路径> <hits pkl> <输出 json>"""
import sys, pickle, json, glob, collections
import numpy as np
import primer_design as pdz
sys.path.insert(0, '/home/user/16S_5R/python_5R/explore')
from offtarget import read_2bit, enc
which, gpath, hitpkl, out = sys.argv[1:5]
SEQ = {'16S-A6-F.2': 'GGCKACACACGTGMTACAA', '16S-A2-R.1': 'GCTGGCACGKARTTAGCCR'}
EXP = {k: pdz.expand(v) for k, v in SEQ.items()}
# (名, 变体名, 5′ 删几个, 3′ 删几个)
VAR = [('16S-A6-F.2', '原始', 0, 0), ('16S-A6-F.2', '5′修剪1', 1, 0), ('16S-A6-F.2', '5′修剪2', 2, 0),
       ('16S-A2-R.1', '原始', 0, 0), ('16S-A2-R.1', '3′修剪1', 0, 1), ('16S-A2-R.1', '3′修剪2', 0, 2)]
hits = [h for h in pickle.load(open(hitpkl, 'rb'))[which] if h[3].split('#')[0] in SEQ and h[4] <= 4]
byc = collections.defaultdict(list)
for h in hits: byc[h[0]].append(h)
if which == 'human': gen = read_2bit(gpath)
else:
    from offtarget_fasta import read_fasta_records
    def g():
        for f in sorted(glob.glob(gpath + '/*.fa')): yield from read_fasta_records(f)
    gen = g()
cnt = collections.defaultdict(set); chk = collections.Counter()
for name, code in gen:
    for c, pos, st, lab, mm, n in byc.get(name, []):
        base, k = lab.split('#'); p = EXP[base][int(k)]
        W = code[pos:pos + n]
        if len(W) < n: continue
        T = enc(p) if st == '+' else enc(pdz.revcomp(p)); d = (W != T)
        chk[int(d.sum() == mm)] += 1
        for b, vn, d5, d3 in VAR:
            if b != base: continue
            # d 是顶链方向；'+' 时 5′ 在左，'-' 时 5′ 在右
            l, r = (d5, d3) if st == '+' else (d3, d5)
            if int(d[l:n - r].sum()) <= 2:
                q = p[d5:len(p) - d3]; cnt[(b, vn)].add((c, pos + (l if st == '+' else l), st, q))
    print(name, dict(chk), flush=True)
json.dump({f'{b}|{v}': len(s) for (b, v), s in cnt.items()} | {'校验(1=重算与原错配一致)': dict(chk)}, open(out, 'w'), ensure_ascii=False, indent=1)
print(open(out).read())
