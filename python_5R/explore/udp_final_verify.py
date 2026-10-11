import sys, json, pickle, collections, itertools, re
import numpy as np
SC = sys.argv[1]; TSV = sys.argv[2]
sel48 = json.load(open(SC + '/selB20.json'))['res']['48']['sel']; sel96 = json.load(open(SC + '/sel96fix.json'))['res']['96']['sel']
U = {l.split()[0]: l.split() for l in open(TSV)}
def lev(a, b):
    p = list(range(len(b) + 1))
    for i, x in enumerate(a, 1):
        c = [i]
        for j, y in enumerate(b, 1): c.append(min(p[j] + 1, c[-1] + 1, p[j - 1] + (x != y)))
        p = c
    return p[-1]
rc = lambda s: s[::-1].translate(str.maketrans('ACGT', 'TGCA'))
R = pickle.load(open(SC + '/udp_eval_all371.pkl', 'rb')); iu = np.triu_indices(len(R['EX']))
T = {}
for d in R['tubes'].values():
    mh = float(d['H'][iu].min()); me = float(d['E'].min())
    T[d['name']] = dict(minH=mh, minE=me, tier='A' if (mh >= -7 and me >= -4.5) else 'B' if (mh >= -8.5 and me >= -5.5) else 'C', bad=float(d['bad']), hp=float(d['hg'].min()))
out = {}
for label, sel in (('48', sel48), ('96', sel96)):
    N = len(sel); i5 = [U[n][3] for n in sel]; i7 = [U[n][1] for n in sel]; r = {}
    r['数量'] = N; r['含 UDP0017'] = 'UDP0017' in sel
    r['i5 最小编辑距离'] = min(lev(a, b) for a, b in itertools.combinations(i5, 2)); r['i7 最小编辑距离'] = min(lev(a, b) for a, b in itertools.combinations(i7, 2))
    r['i5×i7 最小编辑距离（含同管）'] = min(lev(a, b) for a in i5 for b in i7)
    r['i5×rc(i7) 最小编辑距离'] = min(lev(a, rc(b)) for a in i5 for b in i7)
    gcs = [sum(c in 'GC' for c in s) for s in i5 + i7]; r['单个索引 GC(/10) 范围'] = (min(gcs), max(gcs))
    r['同聚物≥3'] = sum(max(len(m.group()) for m in re.finditer(r'(.)\1*', s)) >= 3 for s in i5 + i7)
    r['分级'] = dict(collections.Counter(T[n]['tier'] for n in sel))
    r['全局二聚体最差'] = min(T[n]['minH'] for n in sel); r['3′端最差'] = min(T[n]['minE'] for n in sel); r['发夹最差'] = min(T[n]['hp'] for n in sel)
    cyc = {}
    for key, seqs in (('i5', i5), ('i7', i7)):
        rows = []
        for p in range(10):
            c = collections.Counter(s[p] for s in seqs); gc = c['G'] + c['C']; red = c['A'] + c['C']; green = c['A'] + c['T']
            rows.append(dict(cycle=p + 1, A=c['A'], C=c['C'], G=c['G'], T=c['T'], GC=round(100 * gc / N, 1), AC=round(100 * red / N, 1), AT=round(100 * green / N, 1),
                             ok=bool(0.4 * N <= gc <= 0.6 * N and 0.4 * N <= red <= 0.6 * N and 0.4 * N <= green <= 0.6 * N and min(c['A'], c['C'], c['G'], c['T']) >= 0.15 * N)))
        cyc[key] = rows
    r['每周期'] = cyc; r['每周期全部达标'] = all(x['ok'] for k in cyc for x in cyc[k])
    out[label] = r
    print(label, {k: v for k, v in r.items() if k != '每周期'})
out['48⊂96'] = set(sel48) <= set(sel96); print('48 是 96 的子集', out['48⊂96'])
json.dump(dict(out=out, tubes={n: T[n] for n in sel96}, sel48=sel48, sel96=sel96), open(SC + '/final_verify.json', 'w'), ensure_ascii=False, indent=1)
