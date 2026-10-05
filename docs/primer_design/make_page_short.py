"""生成 docs/primer_shorter_amplicons.html"""
import re, html
import pandas as pd
HERE = __file__.rsplit("/", 1)[0]; DOCS = HERE.rsplit("/", 1)[0]
old = open(DOCS + "/5R_SMURF_pipeline.html", encoding="utf-8").read()
head = re.sub(r"<title>.*?</title>", "<title>缩短 V6V7 和 V8V9</title>", old[:old.index("</style>") + 8], 1)
NUMRE = re.compile(r"^[\d.%/ –→+-]+$")
def cell(v, seq=False): return '<td class="' + ("seq" if seq else ("n" if NUMRE.match(str(v)) else "")) + '">' + html.escape(str(v)) + "</td>"
def tbl(rows, heads, seqcol=None): return '<div class="tbl"><table><thead><tr>' + "".join(f"<th>{c}</th>" for c in heads) + "</tr></thead><tbody>" + "".join("<tr>" + "".join(cell(v, i == seqcol) for i, v in enumerate(r)) + "</tr>" for r in rows) + "</tbody></table></div>"
E = pd.read_csv(f"{HERE}/short_design_eval.csv"); O = pd.read_csv(f"{HERE}/order_list_5amp_short.csv"); H = pd.read_csv(f"{HERE}/short_offtarget_per_oligo.csv").set_index("oligo")
E5 = E[E.design.str.startswith("5")]
t1 = tbl([[r.design.replace("5 扩增子 · ", ""), int(r.oligos), int(r.expansions), f"{r.Tm_min}–{r.Tm_max}", r.spans, r.SILVA_amp, r.GG_amp, f"{r.SILVA_all * 100:.1f}% / {r.GG_all * 100:.1f}%", f"{r.ideal * 100:.2f}%", f"{r.abs60 * 100:.1f} / {r.abs300 * 100:.1f} / {r['frac0.8'] * 100:.1f}%"] for _, r in E5.iterrows()],
         ["方案", "寡核苷酸", "展开数", "Tm °C", "扩增子范围（bp）", "扩增子覆盖 SILVA", "扩增子覆盖 GG", "全部扩出 SILVA / GG", "理想准确率", "属准确率（三种规则）"])
new = O[O.site.str.startswith(("A4", "A5-R"))]
names = {"SMURF5-A4-F.1": "A4-F_0", "SMURF5-A4-F.2": "A4-F_1", "SMURF5-A4-R": "A4-R_0", "SMURF5-A5-R.1": "A5-R_0", "SMURF5-A5-R.2": "A5-R_1"}
G = pd.read_csv(f"{HERE}/gap_design_eval.csv")
t3 = tbl([["原（V6·V7 290 / V8·V9 288）", "10 bp", "8 bp", "否", "否"], ["V6·V7 缩短（229 / 288）", "0", "8 bp", "是", "否"], ["V6·V7 缩短 + V8 缩短（229 / 185）", "0", "0", "是", "是"], ["V6·V7 缩短 + V8·V9 反向 1483（229 / 279）", "0", "0", "是", "是"], ["V6·V7 缩短 + V8·V9 反向 1485（229 / 281）", "0", "1 bp", "是", "否"]],
         ["方案", "V6·V7 缺口", "V8·V9 缺口", "V6·V7 测通", "V8·V9 测通"])
t2 = tbl([[r["name"], r.amplicon, r.seq, r.nt, r.pos, f"{r.Tm}（{r.Tm_lo}–{r.Tm_hi}）", int(r.expansions), f"{H.loc[names[r['name']]].loci_per_expansion:g}"] for _, r in new.iterrows()],
         ["订购名称", "区域", "序列 5′→3′", "nt", "E. coli 位置", "Tm °C（展开范围）", "展开数", "人基因组位点 / 展开序列（≤2 错配）"], seqcol=2)
body = f"""
<div class="wrap"><main style="grid-column:1/-1;width:100%;max-width:980px;margin-inline:auto">
<header>
  <div class="eyebrow">16S_5R · 引物设计</div>
  <h1>缩短 V6·V7（290 bp）和 V8·V9（288 bp）</h1>
  <p class="lede">可以缩短，而且不影响整体质量：V6·V7 缩到 229 bp，V8·V9 缩到 185 bp，覆盖率和属准确率没有下降，代价是多 2 条寡核苷酸，V8·V9 不再包含大部分 V9。这是计算机模拟，没有实验验证。</p>
</header>

<h2>一、结论</h2>
<ul>
<li><b>V6·V7：906–1195（290 bp）→ 967–1195（229 bp）</b>：正向引物移到 967F 附近的保守位点（3′ 端在 985），反向引物不动。理想准确率 95.33% → 95.31%，属准确率（算入扩增失败）持平或略高，需要 2 条正向引物（原来 1 条）。</li>
<li><b>V8·V9：1223–1510（288 bp）→ 1223–1407（185 bp）</b>：反向引物从 1492 附近移到 1389–1407（1391R 附近的保守位点）。覆盖率反而更好（A5 扩增子 SILVA 76% → 79%，GG 94% → 96%），理想准确率 95.31% → 95.19%（−0.12 个百分点），需要 2 条反向引物（原来 1 条）。代价是 1408–1491 不再被测到，V9 区（约 1435–1465）基本丢掉，所以这个扩增子实际是 V8（可算 V8 加一点 V9 的边）。</li>
<li>如果不想丢 V9，只做 V6·V7 的缩短即可（第一行和第二行的差别）。V8·V9 的另一个选择是把反向引物只往里移几个碱基（1485–1491 附近，覆盖率 93–96%），只能缩短不到 25 bp，意义不大。</li>
<li>理想准确率对扩增子内部几乎不敏感：在 928–1176 内部去掉一半，损失 ≤0.15 个百分点，因为其他 4 个扩增子已经把大部分序列分开了；真正的约束是缩短后的引物位点必须仍然保守。</li>
</ul>

<h2>二、比较（5 个扩增子，留出集）</h2>
{t1}
<p>属准确率的三列是三种“扩增失败也算”的比较规则；理想列是假设所有扩增子都扩出。这些指标对小差别不稳（之前有过 4 扩增子理想准确率反而升高的例子），差在 0.1–0.2 个百分点以内的，不要当作真实差别。</p>

<h2>三、改动的引物（其余 10 个位点的引物不变）</h2>
{t2}
<ul>
<li>池内所有寡核苷酸（含展开）两两检查，没有严重二聚体或发夹；平均 Tm 54.8–66.4 °C（原来 54.8–65.8 °C）。A5-R.2 最热，平均 66.4 °C。</li>
<li>人基因组（hg19）：新引物每个展开序列只有 1–3 个近似位点（≤2 个错配）；整池 ≤2 个错配有 2,946 个位点、1 个潜在产物（原来 0 个），≤3 个错配 62 个潜在产物（原来 56 个）；线粒体命中不变，没有潜在产物。</li>
<li>混合比例仍按“每个位点总量相同、位点内按单独覆盖率分”：A4-F 60% / 40%，A5-R 50% / 50%；完整订购清单在 <code>docs/primer_design/order_list_5amp_short.csv</code>。</li>
</ul>

<h2>五、内联 UDP（每端 10 bp）+ 2×150 读长的核对</h2>
<p>每条读段先是 10 bp UDP，再是引物，再是扩增子序列，所以每条读段能读到的扩增子长度是 140 bp（含引物）。两条读段合起来最多覆盖 280 bp；扩增子比这长，中间就有测不到的缺口，缺口长度 = 扩增子长度 − 280。读段对之间没有重叠，所以也不能靠合并读段补上。</p>
{t3}
<ul>
<li><b>原方案</b>：V6·V7（290 bp）缺口 10 bp（1046–1055），V8·V9（288 bp）缺口 8 bp（1363–1370）；其他 3 个扩增子（258 / 186 / 249 bp）都能测通。<b>V6·V7 的缺口在 V6 和 V7 之间的保守间隔区</b>（1044–1116，缺口 10 个碱基平均熵只有 0.15 bit，只有 1 个位点有变异）；<b>V8·V9 的缺口不是保守区</b>（平均熵 0.68 bit，62% 的位点有变异，在 V8 和 V9 之间），更正：我之前把这两个缺口都说成保守区，V8·V9 那个说错了。</li>
<li>把这些缺口位点从评估里去掉，理想准确率 95.33% → 95.32%，属准确率（三种规则）也基本不变（87.9 / 88.4 / 91.0%），所以缺口对属鉴定的影响可以忽略；真正的风险是下游流程需要合并读段（会丢掉这两个扩增子的读段对）。SMURF 把 R1 和反向互补的 R2 直接拼接，不需要重叠。</li>
<li><b>V6·V7 缩到 229 bp</b>：没有缺口，两条读段重叠约 51 bp；<b>V8·V9 缩到 185 bp</b>：没有缺口，重叠约 95 bp，但丢掉 V9。</li>
<li><b>想保留 V9 又要 ≤280 bp</b>：试了把 V8·V9 的反向引物移到 1483（扩增子 279 bp）和 1485（281 bp，缺口 1 bp），但覆盖率明显下降（A5 扩增子 SILVA 74–75%，GG 88–92%，对比原来的 76% / 94%），不推荐。</li>
<li>所有扩增子都比 140 bp 长，所以读段不会读穿到另一端的 UDP 和接头。</li>
</ul>
<h2>六、限制</h2>
<ul>
<li>覆盖率按“最多 1 个错配、3′ 端 3 个碱基匹配”的序列规则算；数据库是 Greengenes 13_8 和 SILVA 128。</li>
<li>缩短后各扩增子长度是 258 / 186 / 249 / 229 / 185 bp，混合比例没有校正长度差异引起的扩增偏倚，需要用已知组成的样本实测。</li>
<li>A5-R 两条的 Tm 偏高（64–66 °C），退火温度要照顾这个范围。</li>
</ul>
<footer>脚本：<code>python_5R/explore/short*.py</code>；数据：<code>docs/primer_design/short_*.csv</code>、<code>order_list_5amp_short.csv</code>。</footer>
</main></div>
"""
open(DOCS + "/primer_shorter_amplicons.html", "w", encoding="utf-8").write(head + body)
print("ok")
