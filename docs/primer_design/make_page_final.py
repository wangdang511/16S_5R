"""生成 docs/primer_compact_primers.html"""
import re, html, json
import numpy as np, pandas as pd
HERE = __file__.rsplit("/", 1)[0]; DOCS = HERE.rsplit("/", 1)[0]
old = open(DOCS + "/5R_SMURF_pipeline.html", encoding="utf-8").read()
head = re.sub(r"<title>.*?</title>", "<title>最少引物设计</title>", old[:old.index("</style>") + 8], 1)
E = pd.read_csv(f"{HERE}/final_compact_eval.csv"); O = json.load(open(f"{HERE}/final_oligos.json"))
NUMRE = re.compile(r"^[\d.%/ –→+-]+$")
def cell(v, seq=False): return '<td class="' + ("seq" if seq else ("n" if NUMRE.match(str(v)) else "")) + '">' + html.escape(str(v)) + "</td>"
def tbl(rows, heads, seqcol=None): return '<div class="tbl"><table><thead><tr>' + "".join(f"<th>{c}</th>" for c in heads) + "</tr></thead><tbody>" + "".join("<tr>" + "".join(cell(v, i == seqcol) for i, v in enumerate(r)) + "</tr>" for r in rows) + "</tbody></table></div>"
pc = lambda x: f"{x * 100:.0f}%"; p1 = lambda x: f"{x * 100:.1f}%"
rows = []
order = [("5 扩增子", 0.95), ("5 扩增子", 0.9), ("4 扩增子", 0.95), ("4 扩增子", 0.9), ("3 扩增子", 0.95), ("3 扩增子", 0.9)]
for sch, T in order:
    r = E[(E.scheme == sch) & (E.target == T)]
    if len(r) == 0:
        rows.append([sch, f"{T:.0%}", "不可行：V1·V2 的反向引物 ≤6 条寡核苷酸达不到 95%", "", "", "", "", "", "", "", ""]); continue
    r = r.iloc[0]
    rows.append([sch, f"{T:.0%}", int(r.oligos), int(r.expansions), f"{r.Tm_min:.0f}–{r.Tm_max:.0f}", int(r.severe), f"{pc(r.SILVA_p_all)} / {pc(r.GG_p_all)}", p1(r.ideal_acc), p1(r.realistic_acc),
                 r.SILVA_amp_cov, r.GG_amp_cov])
t_cmp = tbl(rows, ["扩增子数", "每个位点覆盖率目标", "寡核苷酸数", "展开后的序列数", "Tm 范围 °C", "严重二聚体", "全部扩增子都扩出：SILVA / GG", "属准确率（理想）", "属鉴定正确率（含扩增失败）", "各扩增子覆盖 %（SILVA）", "各扩增子覆盖 %（GG）"])
def oligo_tbl(key):
    rows = [[o["site"], o["pos"], o["seq"], o["length"], o["expansions"], o["GC"], f'{o["Tm"]}（{o["Tm_range"]}）', o["worst_dimer_dG"]] for o in O[key]["oligos"]]
    return tbl(rows, ["位点", "E. coli 位置", "序列 5′→3′", "nt", "展开数", "GC", "Tm °C（展开范围）", "与池内最差异二聚体 ΔG kcal/mol"], seqcol=2)
r4 = E[(E.scheme == "4 扩增子") & (E.target == 0.95)].iloc[0]; r3 = E[(E.scheme == "3 扩增子") & (E.target == 0.95)].iloc[0]
body = f"""
<div class="wrap"><main style="grid-column:1/-1;width:100%;max-width:980px;margin-inline:auto">
<header>
  <div class="eyebrow">16S_5R · 引物设计</div>
  <h1>覆盖率目标下的最少引物：5、4、3 个扩增子的比较</h1>
  <p class="lede">沿用“每个位点用简并引物集合覆盖序列”的做法：要求每个位点覆盖率达到目标，找最少的寡核苷酸，再均衡 Tm、用热力学计算检查二聚体。覆盖率目标分别取 90%（你要求的）和 95%（之前的）。</p>
  <p>这是计算机模拟，没有实验验证。数据库是 Greengenes 13_8 和 SILVA 128，不是最新版本。</p>
</header>

<h2>一、结论</h2>
<ul>
<li><b>90% 并没有比 95% 省引物。</b>4 个扩增子两个目标下都是 10 条寡核苷酸，3 个扩增子都是 8 条；95% 的设计覆盖更好，属鉴定正确率也更高（4 个扩增子 {p1(r4.realistic_acc)}，对 87.0%）。90% 的最小集合倾向于选覆盖刚好达标但不稳的位点，在验证集上 V4 扩增子（A3）在 Greengenes 里只有 50%。</li>
<li><b>5 个扩增子不能用最少引物实现 95%</b>：V1·V2 扩增子的反向引物（A1-R，243–260 附近）用 ≤6 条寡核苷酸也达不到 95%；90% 下需要 3 条，5 个扩增子合计 14 条，但“全部扩增子都扩出”的序列只有 42%（SILVA）/ 27%（GG）。</li>
<li><b>推荐 4 个扩增子、95% 目标、10 条寡核苷酸</b>（V3、V4、V6·V7、V8·V9）：属准确率（理想）{p1(r4.ideal_acc)}，把扩增失败算进去 {p1(r4.realistic_acc)}。比之前每位点 2–3 条的 4 扩增子设计（24 条寡核苷酸，90.8%）少了一半多的引物，鉴定正确率还略高。</li>
<li><b>极简版 3 个扩增子（V3、V4、V8·V9）</b>：8 条寡核苷酸，{p1(r3.ideal_acc)} / {p1(r3.realistic_acc)}，与 4 个扩增子只差约 0.3 / 0.7 个百分点。</li>
<li>剩下的问题：A5-F（V8·V9 的正向引物）的两条寡核苷酸 Tm 只有 55.1 和 50.8 °C，是整个池里最低的，我没有找到不损失覆盖率的办法抬高它们（见第四节）。</li>
</ul>

<h2>二、方法</h2>
<ul>
<li><b>覆盖率</b>：12 个主要细菌门等权的平均覆盖率；序列与引物最多 1 个错配，3′ 端 3 个碱基必须匹配。设计用 Greengenes（半数）+ SILVA（1.2 万条）；每个库的覆盖率都要达标，不再混合后平均。验证用没参与设计的 SILVA（10,443 条）和 Greengenes（3,151 条）。</li>
<li><b>每个位点的最小集合</b>：3′ 端位置允许 ±2 nt，长度 16–24 nt；贪心加入简并引物（每条 ≤4 倍简并且 ≤2 个简并位置，或 ≤8 倍且 ≤3 个），取满足目标的最少寡核苷酸数。</li>
<li><b>池优化</b>：各位点候选配置里做坐标下降，目标是严重二聚体最少、Tm 靠近 60.5 °C、寡核苷酸和展开序列最少；扩增子 150–290 bp、同管不重叠、间隔 ≥20 nt。</li>
<li><b>热力学检查</b>用 primer3（50 mM Na⁺、2 mM Mg²⁺、0.2 mM dNTP、250 nM 引物，60 °C），对所有简并展开序列两两检查：异二聚体 ΔG ≤ −9 kcal/mol 或 3′ 端 ΔG ≤ −6 kcal/mol 算严重；发夹 ΔG ≤ −3 kcal/mol 算有问题。</li>
<li><b>Tm 精调</b>：把 Tm 低于 56 °C 的寡核苷酸沿模板在 5′ 端延长（新增碱基取它所覆盖序列的多数碱基），只在覆盖率不掉时采用。A3-F 的第 2 条从 50.4 升到 56.5 °C，覆盖率基本不变，采用；A5-F 延长后覆盖率掉 5 个百分点，不采用。</li>
</ul>

<h2>三、比较</h2>
{t_cmp}
<p>“属鉴定正确率（含扩增失败）”用 2,801 条没参与设计的序列：序列只在它能扩出的扩增子上有数据，没扩出的当作没测到。绝对值比理想情况低很多，只有相对排序有参考意义。</p>
<ul>
<li>90% 的设计里，V4 扩增子在 Greengenes 验证集上只有 50%（A3-F 从 556 起，Greengenes 里这一段有长度多态，之前发现过）；设计集里它的覆盖率是达标的，我没能完全解释设计集和验证集之间的这个差距。95% 的设计把 V4 正向引物的第 1 条放在 558，验证集上正常；第 2 条为抬高 Tm 向 5′ 端延长了 2 nt，所以从 556 起，V4 扩增子在 Greengenes 验证集上有 91%。</li>
<li>所有 5 个设计的相邻外侧引物跨区产物最短都是 V3 正向 + V4 反向，约 464–468 bp，没有短于 460 bp 的跨区产物。</li>
</ul>

<h2>四、推荐设计：4 个扩增子，95%，10 条寡核苷酸</h2>
<p>扩增子：341–527（187 bp，V3）、556–804（249 bp，V4）、906–1192（287 bp，V6·V7）、1226–1510（285 bp，V8·V9），都 ≤290 bp，需要 2×150 读长。</p>
{oligo_tbl("4 扩增子|0.95")}
<ul>
<li>池内没有严重二聚体，没有发夹。最差的异二聚体是 A5-F 的第 1 条（ΔG −8.4 kcal/mol），接近阈值，需要用 OligoAnalyzer 复核。</li>
<li><b>Tm 不均匀</b>：50.8–63.5 °C。最低的是 A5-F 的两条（55.1、50.8）和 A3-F 的两条（57.4、56.5）。A5-F 在 95% 下只有 1 个满足覆盖率的配置，没有别的选择；延长 5′ 端会掉覆盖率。可行的办法是退火温度取偏低值（约 55 °C），或者对 A5-F 用高于 16 nt 的更长引物并接受覆盖率下降，或者换到 90% 目标下的 A5-F 配置（Tm 约 58–62 °C，但设计集覆盖率只有 93%–94%）。</li>
<li>简并碱基很多：10 条寡核苷酸展开成 56 种序列，每种序列在混合物里的有效浓度更低，需要提高对应寡核苷酸的浓度。</li>
</ul>
<h3>3 个扩增子的极简版</h3>
<p>去掉 V6·V7 扩增子（906–1192），其余 8 条寡核苷酸不变，所以序列与上表相同，只是去掉 A4-F、A4-R 两个位点。{p1(r3.ideal_acc)} / {p1(r3.realistic_acc)}。</p>

<h2>五、局限</h2>
<ul>
<li>属准确率是代理指标，只测到属；Greengenes 属标签偏向培养过的属。</li>
<li>覆盖率是序列匹配规则的结果，不是实际 PCR 效率；热力学检查是基于 primer3 的近似。没有检查对人基因组和线粒体的特异性。</li>
<li>池里的 A5-F 是已知的弱点；V1·V2 在这套方案里被去掉了，所以没有 V1·V2 的信息。</li>
<li>数据库不是最新版本，GTDB 没评估；k-mer 数据库需要用新引物重建，SMURF 的区域数（<code>nR = 5</code>）要改。</li>
</ul>
<footer>脚本：<code>python_5R/compact_primers.py</code>、<code>python_5R/explore/ 下的 run_compact.py、pool_opt.py、refine.py、final_eval.py</code>；数据：<code>docs/primer_design/final_compact_eval.csv</code>、<code>final_oligos.json</code>。</footer>
</main></div>
"""
open(DOCS + "/primer_compact_primers.html", "w", encoding="utf-8").write(head + body)
print("ok")
