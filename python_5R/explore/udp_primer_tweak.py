"""对池里一条引物做最小改动（5′端修剪 / 5′端首位替换 / 3′端修剪），评估：
 1) 位点覆盖率（SILVA、Greengenes，主要门平均，≤1 错配且 3′ 端 3 nt 无错配，与原设计同一口径）与 Tm
 2) 在最终 96 管里与伙伴引物的尾巴相互作用（全局二聚体、3′端稳定性）
 3) 把改动后的引物换进整管后，96 管重新分级
注：覆盖率部分依赖会话 scratchpad 里的 SILVA / Greengenes 参考（final_eval.py 的预处理），路径需自行替换。
用法：python3 udp_primer_tweak.py <pairs json> <final eval pkl> <out json> [full]
  不带 full 只做 1)+2)；带 full 对通过筛选的改动做 3)。"""
import sys, os, json, pickle, itertools
import numpy as np, openpyxl, primer3
from multiprocessing import Pool
sys.path.insert(0, '/home/user/16S_5R/python_5R'); sys.path.insert(0, '/home/user/16S_5R/python_5R/explore')
SC = '/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
pj, evp, outp = sys.argv[1:4]; FULL = len(sys.argv) > 4 and sys.argv[4] == 'full'
KW = dict(mv_conc=50, dv_conc=2, dntp_conc=0.2, dna_conc=250); T = 58.0
import warnings; warnings.filterwarnings('ignore')
here = os.getcwd(); os.chdir(SC)
exec(open(SC + 'final_eval.py').read().split('def evaluate(key,sites)')[0])   # 提供 VAL、site_hit、pdz、cp（oligo_tm / pool_issues）
os.chdir(here)

ws = openpyxl.load_workbook('/home/user/16S_5R/docs/primer_design/primer_order_6amp_pool33_v3.xlsx')['订购清单']
OL = {}
for r in range(4, 37):
    nm, site, seq, pos = (ws.cell(r, c).value for c in (2, 3, 5, 6)); a, b = [int(x) for x in str(pos).replace('–', '-').split('-')]
    OL[nm] = dict(site=site, seq=seq, start=a, L=b - a + 1, orient=site[-1])

def spec(o, seq, start, L): return dict(orient=o, start=start, L=L, prim=[seq], ext=[0])
def site_cov(site, repl=None):
    specs = []
    for nm, o in OL.items():
        if o['site'] != site: continue
        if repl and nm == repl[0]: specs.append(spec(o['orient'], *repl[1:]))
        else: specs.append(spec(o['orient'], o['seq'], o['start'], o['L']))
    out = {}
    for vn, (ref, rows_) in VAL.items():
        phs = ref.phylum[rows_]; h = np.zeros(len(rows_), bool)
        for s in specs: h |= site_hit(ref, rows_, s)
        out[vn] = round(float(np.mean([h[phs == g].mean() for g in pdz.KEY_PHYLA if (phs == g).sum() >= 30])), 3)
    return out

def variants(nm):
    o = OL[nm]; s = o['seq']; v = {'原始': (s, o['start'], o['L'])}
    if o['orient'] == 'F':           # 5′ 端在左
        for k in (1, 2, 3): v[f'5′修剪{k}'] = (s[k:], o['start'] + k, o['L'] - k)
        for b in 'ACT': v[f'5′首位G→{b}'] = (b + s[1:], o['start'], o['L'])
        for b in 'ACT': v[f'5′第二位→{b}'] = (s[0] + b + s[2:], o['start'], o['L']) if s[1] != b else None
    else:                            # 3′ 端在左
        for k in (1, 2): v[f'3′修剪{k}'] = (s[:-k], o['start'] + k, o['L'] - k)
        for k in (1, 2, 3): v[f'5′修剪{k}'] = (s[k:], o['start'], o['L'] - k)
    return {k: x for k, x in v.items() if x}

pairs = json.load(open(pj))['pairs']; ev = pickle.load(open(evp, 'rb')); EX = ev['EX']
tubes = [t for t in ev['tubes'].values()]
tails = [(p['seq5'], p['seq7']) for p in pairs]            # 管 k 的 (i5, i7)，与 R001.. 同序
KEEP = lambda nm: [i for i, x in enumerate(EX) if x[0] != nm]
exp_of = pdz.expand

def block(args):
    k, rows_new, ori_new, keep_idx = args       # rows_new: 新展开序列；ori_new: 'F'/'R'
    t5, t7 = tails[k]; tl = lambda o, s: (t5 if o == 'F' else t7) + s
    kept = [tl(EX[i][1], EX[i][2]) for i in keep_idx]; new = [tl(ori_new, s) for s in rows_new]
    Hm = []; Em = []; Ek = []
    for a in new:
        Hm.append([primer3.calc_heterodimer(a, b, temp_c=T, **KW).dg / 1000 for b in kept])
        Em.append([primer3.calc_end_stability(a, b, temp_c=T, **KW).dg / 1000 for b in kept])
        Ek.append([primer3.calc_end_stability(b, a, temp_c=T, **KW).dg / 1000 for b in kept])
    self_h = [primer3.calc_heterodimer(a, b, temp_c=T, **KW).dg / 1000 for a in new for b in new]
    self_e = [primer3.calc_end_stability(a, b, temp_c=T, **KW).dg / 1000 for a in new for b in new]
    return k, min(min(map(min, Hm)), min(self_h)), min(min(map(min, Em)), min(map(min, Ek)), min(self_e))

tier = lambda h, e: 'A' if h >= -7 and e >= -4.5 else 'B' if h >= -8.5 and e >= -5.5 else 'C'
N = len(EX); iu = np.triu_indices(N); res = {}
base_tiers = [tier(float(t['H'][iu].min()), float(t['E'].min())) for t in tubes]
print('基线', {x: base_tiers.count(x) for x in 'ABC'}, flush=True)
TARGETS = [('16S-A6-F.2', '16S-A2-R.1'), ('16S-A2-R.1', '16S-A6-F.2')]
if __name__ == '__main__':
    pool = Pool(4)
    for nm, partner in TARGETS:
        o = OL[nm]; keep_idx = KEEP(nm); part_idx = [i for i in keep_idx if EX[i][0] == partner]
        for vn, (seq, st, L) in variants(nm).items():
            key = f'{nm}|{vn}'; cov = site_cov(o['site'], (nm, seq, st, L)); tm = [round(cp.oligo_tm(x)[0], 1) for x in [seq]]
            new_rows = exp_of(seq)
            # 2) 只看与伙伴引物的相互作用
            jobs = [(k, new_rows, o['orient'], part_idx) for k in range(len(tubes))]
            out = pool.map(block, jobs)
            hmin = np.array([x[1] for x in out]); emin = np.array([x[2] for x in out])
            res[key] = dict(seq=seq, start=st, L=L, Tm=tm, cov=cov, nexp=len(new_rows), pair_E_lt45=int((emin < -4.5).sum()), pair_E_min=float(emin.min()), pair_H_lt7=int((hmin < -7).sum()), pair_H_min=float(hmin.min()))
            print(key, seq, res[key], flush=True)
    if FULL:
        sel = [k for k, v in res.items() if not k.endswith('原始') and v['pair_E_lt45'] <= 3]
        for key in sel:
            nm = key.split('|')[0]; o = OL[nm]; keep_idx = KEEP(nm); new_rows = exp_of(res[key]['seq'])
            jobs = [(k, new_rows, o['orient'], keep_idx) for k in range(len(tubes))]
            out = pool.map(block, jobs); ts = []
            for k, h, e in out:
                kept_h = tubes[k]['H'][np.ix_(keep_idx, keep_idx)]; kept_e = tubes[k]['E'][np.ix_(keep_idx, keep_idx)]
                ts.append(tier(min(float(kept_h[np.triu_indices(len(keep_idx))].min()), h), min(float(kept_e.min()), e)))
            res[key]['tiers_full'] = {x: ts.count(x) for x in 'ABC'}; print(key, res[key]['tiers_full'], flush=True)
    json.dump(res, open(outp, 'w'), ensure_ascii=False, indent=1)
