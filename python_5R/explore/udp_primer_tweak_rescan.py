"""3′ 端修剪的变体需要重新扫描基因组：原搜索要求 3′ 端 8 nt 完全匹配，3′ 末位错配的位点不在 hits 里。
对变体自己的展开序列直接做全基因组搜索（≤2 错配，3′ 端 8 nt 种子），按条目 (染色体, 位置, 链, 展开序列) 计数。
用法：PYTHONPATH=python_5R python3 udp_primer_tweak_rescan.py <human|pig> <基因组路径> <输出 json>"""
import sys, glob, json, collections
import primer_design as pdz
sys.path.insert(0, '/home/user/16S_5R/python_5R/explore')
from offtarget import read_2bit, search
which, gpath, out = sys.argv[1:4]
import os
V = json.loads(os.environ['VARS']) if os.environ.get('VARS') else {'16S-A2-R.1|原始': 'GCTGGCACGKARTTAGCCR', '16S-A2-R.1|3′修剪1': 'GCTGGCACGKARTTAGCC', '16S-A2-R.1|3′修剪2': 'GCTGGCACGKARTTAGC'}
prims = {f'{k}#{i}': e for k, s in V.items() for i, e in enumerate(pdz.expand(s))}
if which == 'human': gen = read_2bit(gpath)
else:
    from offtarget_fasta import read_fasta_records
    def g():
        for f in sorted(glob.glob(gpath + '/*.fa')): yield from read_fasta_records(f)
    gen = g()
seen = collections.defaultdict(set)
for name, code in gen:
    for c, pos, st, lab, mm, n in search(code, name, prims, 2):
        seen[lab.rsplit('#', 1)[0]].add((c, pos, st, prims[lab]))
    print(name, {k: len(v) for k, v in seen.items()}, flush=True)
json.dump({k: len(v) for k, v in seen.items()}, open(out, 'w'), ensure_ascii=False, indent=1)
