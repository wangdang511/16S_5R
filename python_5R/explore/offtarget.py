"""引物对人基因组（hg19，含 chrM）和 rCRS 线粒体的脱靶检查。
用法：PYTHONPATH=python_5R python3 explore/offtarget.py <hg19 2bit> <rCRS fasta> <输出目录>
规则：引物 3′ 端 8 nt 必须完全匹配，整体错配 ≤MAXMM；正向/反向两个方向都找；
再把相对方向、3′ 端相对的命中配成潜在产物（≤2000 bp）。"""
import sys, json, struct, itertools, collections
import numpy as np
import primer_design as pdz

def read_2bit(path):
    f = open(path, "rb"); sig = f.read(4)
    le = struct.unpack("<I", sig)[0] == 0x1A412743; E = "<" if le else ">"
    ver, n, _ = struct.unpack(E + "3I", f.read(12)); idx = []
    for _ in range(n):
        ln = f.read(1)[0]; name = f.read(ln).decode(); off = struct.unpack(E + "I", f.read(4))[0]; idx.append((name, off))
    for name, off in idx:
        f.seek(off); size, nb = struct.unpack(E + "2I", f.read(8))
        ns = np.frombuffer(f.read(4 * nb), dtype=E + "u4"); nz = np.frombuffer(f.read(4 * nb), dtype=E + "u4")
        mb = struct.unpack(E + "I", f.read(4))[0]; f.read(8 * mb + 4)
        packed = np.frombuffer(f.read((size + 3) // 4), dtype=np.uint8)
        a = np.empty(len(packed) * 4, np.uint8)
        a[0::4] = packed >> 6; a[1::4] = (packed >> 4) & 3; a[2::4] = (packed >> 2) & 3; a[3::4] = packed & 3
        a = a[:size]; code = np.array([3, 1, 0, 2], np.uint8)[a]   # 2bit: T=0 C=1 A=2 G=3 -> A=0 C=1 G=2 T=3
        for s, z in zip(ns, nz): code[s:s + z] = 4
        yield name, code

def read_fa(path):
    seq = "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()
    m = {"A": 0, "C": 1, "G": 2, "T": 3}; return np.array([m.get(c, 4) for c in seq], np.uint8)

ENC = {"A": 0, "C": 1, "G": 2, "T": 3}
def enc(s): return np.array([ENC[c] for c in s], np.uint8)

def search(code, name, prims, MAXMM, SEED=8):
    """prims: dict label -> (seq 5'->3'); 每条是展开后的非简并序列。返回命中列表。"""
    L = len(code); out = []
    # 种子表：3′ 端 SEED nt。+ 方向：窗口末 SEED nt = 引物 3′ 端；− 方向：窗口起始 SEED nt = rc(引物)[:SEED]
    for strand in "+-":
        table = collections.defaultdict(list)
        for lab, p in prims.items():
            q = p if strand == "+" else pdz.revcomp(p)
            seed = q[-SEED:] if strand == "+" else q[:SEED]
            table[int("".join(str(ENC[c]) for c in seed), 5)].append((lab, q, enc(q)))
        keys = np.array(sorted(table), np.int64); lut = np.zeros(5 ** SEED, np.int32) - 1
        lut[keys] = np.arange(len(keys))
        step = 20_000_000
        for s0 in range(0, L, step):
            e0 = min(L, s0 + step + 40); c = code[s0:e0].astype(np.int64)
            if len(c) < SEED: continue
            h = np.zeros(len(c) - SEED + 1, np.int64)
            for i in range(SEED): h = h * 5 + c[i:len(c) - SEED + 1 + i]
            hit = np.flatnonzero(lut[h] >= 0); hit = hit[hit < step]
            if len(hit) == 0: continue
            for pos in hit:
                grp = table[int(keys[lut[h[pos]]])]
                for lab, q, eq in grp:
                    n = len(q)
                    w_start = s0 + pos - (n - SEED) if strand == "+" else s0 + pos
                    if w_start < 0 or w_start + n > L: continue
                    w = code[w_start:w_start + n]
                    mm = int((w != eq).sum())
                    if mm <= MAXMM: out.append((name, int(w_start), strand, lab, mm, n))
    return out

def pairs(hits, maxlen=2000, maxmm=4):
    """+ 方向命中（3′ 端在窗口右端）与其下游的 − 方向命中（3′ 端在窗口左端）配成产物。"""
    res = []; by = collections.defaultdict(lambda: ([], []))
    for name, pos, st, lab, mm, n in hits:
        if mm > maxmm: continue
        by[name][0 if st == "+" else 1].append((pos, n, lab, mm))
    for name, (P, M) in by.items():
        M = sorted(M); P = sorted(P); ms = [m[0] for m in M]
        import bisect
        for pos, n, lab, mm in P:
            end = pos + n
            i = bisect.bisect_left(ms, end)
            while i < len(M) and M[i][0] + M[i][1] - pos <= maxlen:
                res.append((name, pos, M[i][0] + M[i][1], M[i][0] + M[i][1] - pos, lab, mm, M[i][2], M[i][3])); i += 1
    return res

if __name__ == "__main__":
    twobit, rcrs, outdir = sys.argv[1:4]; MAXMM = 5; ONLY = sys.argv[4] if len(sys.argv) > 4 else None
    sets = json.load(open(outdir + "/offtarget_primers.json"))
    if ONLY: sets = {ONLY: sets[ONLY]}
    exp = {}
    for sname, ol in sets.items():
        d = {}
        for o in ol:
            for k, p in enumerate(pdz.expand(o["seq"])): d[f'{o["name"]}#{k}'] = p
        exp[sname] = d
    res = {}
    for sname, d in exp.items():
        hits = search(read_fa(rcrs), "rCRS", d, MAXMM)
        res[(sname, "rCRS")] = hits
        print(sname, "rCRS hits ≤%d mm:" % MAXMM, len(hits), flush=True)
    allhits = {s: [] for s in exp}
    for name, code in read_2bit(twobit):
        if "_" in name and not name.endswith("_hap1"): pass
        for sname, d in exp.items():
            h = search(code, name, d, MAXMM); allhits[sname] += h
        print(name, len(code), {s: len(v) for s, v in allhits.items()}, flush=True)
    import pickle; pickle.dump((res, allhits), open(outdir + "/offtarget_hits_%s.pkl" % (ONLY or "all"), "wb"))
