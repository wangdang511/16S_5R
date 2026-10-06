"""汇总 offtarget.py 的输出 pkl：每条引物在目标基因组上 ≤2/≤3 错配的位点数，以及线粒体命中。
用法：python3 summarize_offtarget.py <offtarget_hits_pool36.pkl> <offtarget_primers.json> <输出csv>
需要在仓库根目录用 PYTHONPATH=python_5R 运行。"""
import sys, pickle, json, collections, csv
import primer_design as pdz
pkl, js, out = sys.argv[1:4]
res, allh = pickle.load(open(pkl, "rb")); sets = json.load(open(js))
sname = next(iter(sets)); hits = allh[sname]; mt = res[(sname, next(k[1] for k in res if k[0] == sname))]
c2 = collections.Counter(h[3].split("#")[0] for h in hits if h[4] <= 2)
c3 = collections.Counter(h[3].split("#")[0] for h in hits if h[4] <= 3)
m3 = collections.Counter(h[3].split("#")[0] for h in mt if h[4] <= 3)
with open(out, "w", newline="", encoding="utf-8-sig") as f:
    w = csv.writer(f); w.writerow(["名称", "序列", "展开数", "≤2错配位点", "每展开", "≤3错配位点", "线粒体≤3错配"])
    for o in sets[sname]:
        e = len(pdz.expand(o["seq"])); n = o["name"]
        w.writerow([n, o["seq"], e, c2[n], round(c2[n] / e, 1), c3[n], m3[n]])
print("写入", out)
