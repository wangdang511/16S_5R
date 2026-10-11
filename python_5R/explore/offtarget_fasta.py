"""与 offtarget.py 同一套规则（3′端 8 nt 完全匹配、整体错配 ≤5），但输入是按染色体拆分的 fasta 文件，适用于没有 2bit 的基因组（如猪）。
用法：PYTHONPATH=python_5R python3 python_5R/explore/offtarget_fasta.py <fasta 目录或单个 fasta> <线粒体染色体名，如 chrM> <输出目录> [集合名]
输出：offtarget_hits_<集合名>.pkl，与 offtarget.py 的格式相同（可直接用 summarize_offtarget.py 汇总）。"""
import sys, os, json, glob, pickle
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import primer_design as pdz
from offtarget import search, pairs, ENC

M4 = np.full(256, 4, np.uint8)
for c, v in zip("ACGTacgt", [0, 1, 2, 3, 0, 1, 2, 3]): M4[ord(c)] = v

def read_fasta_records(path):
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line[0] == ">":
                if name is not None: yield name, np.concatenate(buf) if buf else np.zeros(0, np.uint8)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(M4[np.frombuffer(line.strip().encode(), np.uint8)])
        if name is not None: yield name, np.concatenate(buf) if buf else np.zeros(0, np.uint8)

if __name__ == "__main__":
    src, mito, outdir = sys.argv[1:4]; SET = sys.argv[4] if len(sys.argv) > 4 else None; MAXMM = 5
    sets = json.load(open(outdir + "/offtarget_primers.json"))
    if SET: sets = {SET: sets[SET]}
    exp = {}
    for sname, ol in sets.items():
        d = {}
        for o in ol:
            for k, p in enumerate(pdz.expand(o["seq"])): d[f'{o["name"]}#{k}'] = p
        exp[sname] = d
    files = sorted(glob.glob(os.path.join(src, "*.fa*"))) if os.path.isdir(src) else [src]
    res = {}; allhits = {s: [] for s in exp}
    for f in files:
        for name, code in read_fasta_records(f):
            for sname, d in exp.items():
                h = search(code, name, d, MAXMM); allhits[sname] += h
                if name == mito: res[(sname, "mito")] = h
            print(name, len(code), {s: len(v) for s, v in allhits.items()}, flush=True)
    pickle.dump((res, allhits), open(outdir + "/offtarget_hits_%s.pkl" % (SET or "all"), "wb"))
