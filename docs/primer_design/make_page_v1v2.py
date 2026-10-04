"""生成 docs/primer_v1v2.html"""
import re, html
import pandas as pd
HERE = __file__.rsplit("/", 1)[0]; DOCS = HERE.rsplit("/", 1)[0]
old = open(DOCS + "/5R_SMURF_pipeline.html", encoding="utf-8").read()
head = re.sub(r"<title>.*?</title>", "<title>V1V2 引物单独优化</title>", old[:old.index("</style>") + 8], 1)
NUMRE = re.compile(r"^[\d.%/ –→+-]+$")
def cell(v, seq=False): return '<td class="' + ("seq" if seq else ("n" if NUMRE.match(str(v)) else "")) + '">' + html.escape(str(v)) + "</td>"
def tbl(rows, heads, seqcol=None): return '<div class="tbl"><table><thead><tr>' + "".join(f"<th>{c}</th>" for c in heads) + "</tr></thead><tbody>" + "".join("<tr>" + "".join(cell(v, i == seqcol) for i, v in enumerate(r)) + "</tr>" for r in rows) + "</tbody></table></div>"
p1 = lambda x: f"{x * 100:.1f}%"
R = pd.read_csv(f"{HERE}/v1v2_scan_R_all.csv"); F = pd.read_csv(f"{HERE}/v1v2_scan_F_all.csv")
def best(D, k, fold, side):
    b = D[(D.k == k) & (D.fold == fold)].sort_values("cov_min", ascending=False).iloc[0]
    pos = f"{b.start}–{b.start + b.L - 1}"
    return [f"{k}", f"≤{fold}", int(b.total_fold), pos, f"{b.cov_min * 100:.1f}%", b.tops]
rowsR = [best(R, k, f, "R") for k, f in ((1, 16), (2, 8), (2, 16), (3, 8), (3, 16), (4, 8), (4, 16), (6, 16))]
rowsF = [best(F, k, f, "F") for k, f in ((1, 8), (2, 8), (3, 8), (4, 16))]
A = pd.read_csv(f"{HERE}/v1v2_amplicon_eval.csv")
rowsA = [[r.F, r.A1R.split(" + ", 1)[1], int(r.oligos), int(r.exp), f"{r.Tm_min:.0f}–{r.Tm_max:.0f}", f"{r.SILVA_amp}% / {r.GG_amp}%", p1(r.held_amp)] for _, r in A.iterrows() if r.F != "F17(2条)"]
D = pd.read_csv(f"{HERE}/v1v2_design_eval.csv"); RB = pd.read_csv(f"{HERE}/robust_realistic.csv")
rowsD = [[r.design, p1(r.abs60), p1(r.abs300), p1(r["frac0.8"]), p1(r.ideal_on_set)] for _, r in RB.iterrows() if r.design in ("S1", "紧凑 4 扩增子 95%", "紧凑 4 扩增子 95% + V1V2")]
rowsD = [rowsD[0], rowsD[1]] + [[r.design.replace("紧凑 4 扩增子 95% + V1V2", "4 扩增子 + V1V2"), p1(r.abs60), p1(r.abs300), p1(r["frac0.8"]), p1(r.ideal_on_set)] for _, r in D.iterrows() if r.design.startswith("紧凑") and ("F17(1条) + 2 条/32" in r.design or "F17(1条) + 4 条/64" in r.design)]
rowsD += [[r.design.replace("S1 + V1V2", "S1 + V1V2"), p1(r.abs60), p1(r.abs300), p1(r["frac0.8"]), p1(r.ideal_on_set)] for _, r in D.iterrows() if r.design.startswith("S1") and "F17(1条) + 4 条/64" in r.design]
tR = tbl(rowsR, ["引物条数", "单条展开上限", "展开后总数", "位置（E. coli）", "位点覆盖率（Greengenes / SILVA 取较低）", "序列 5′→3′（顶链方向）"], seqcol=5)
tF = tbl(rowsF, ["引物条数", "单条展开上限", "展开后总数", "位置", "位点覆盖率", "序列 5′→3′"], seqcol=5)
tA = tbl(rowsA, ["A1-F", "A1-R", "寡核苷酸数", "展开后总数", "Tm °C", "扩增子覆盖：SILVA / GG", "留出集覆盖"])
tD = tbl(rowsD, ["方案", "共同位点 ≥60", "≥300", "≥80% 自身位点", "理想（全部扩出）"])
body = f"""
<div class="wrap"><main style="grid-column:1/-1;width:100%;max-width:980px;margin-inline:auto">
<header>
  <div class="eyebrow">16S_5R · 引物设计</div>
  <h1>V1·V2 反向引物单独优化：值不值得把 V1·V2 加回去</h1>
  <p class="lede">单独为 V1·V2（A1）的反向引物扫描位置、长度、条数和简并度，并把正向引物 A1-F 也一起重做。位点覆盖率能提到 92–95%，但整个扩增子的覆盖率几乎不动，所以不建议把 V1·V2 加回推荐方案。</p>
  <p>这是计算机模拟，没有实验验证。数据库是 Greengenes 13_8 和 SILVA 128，不是最新版本。</p>
</header>

<h2>一、结论</h2>
<ul>
<li><b>A1-R 的位点覆盖率可以明显提高</b>：原来 3 条（展开 16）的组合，改成 2 条（展开 32，259–275）约 92%，3 条（展开 48，244–260）约 93%，4 条（展开 64，243–259）约 95%。位点覆盖率是 Greengenes 和 SILVA 中较低的那个。</li>
<li><b>但整个 V1·V2 扩增子的覆盖率只提高了几个百分点</b>：SILVA 从 74% 到 77–78%，Greengenes 从 87% 到 90–92%，而且 7–8 条寡核苷酸才到顶。正向引物 A1-F 一侧本身就有缺口，SILVA 里 A1-F 单独只有约 80% 的序列匹配（全部序列，不限于主要门）。</li>
<li><b>A1-F 缩短到 17 nt（8–24）更好</b>：1 条 <code>AGRGTTTGATYMTGGCT</code> 位点覆盖率 95%，比原来的 19 nt（8–26）高；扩增子覆盖率只多约 1 个百分点。</li>
<li><b>把 V1·V2 加回去，属准确率没有净收益</b>：理想情况（全部扩出）提高约 0.7–0.9 个百分点；把扩增失败算进去后，4 个扩增子方案加 V1·V2 在三种比较规则下是 90.0 / 90.1 / 91.4%，不加是 91.4 / 90.9 / 91.5%，持平或略差。</li>
<li><b>推荐不变</b>：4 个扩增子（V3、V4、V6·V7、V8·V9）、95%、10 条寡核苷酸。V1·V2 如果一定要，A1-F 用 17 nt 单条，A1-R 用 4 条（展开 64）；代价是再多 6 条寡核苷酸，Tm 下限 53 °C 左右。</li>
</ul>

<h2>二、更正：上一轮“加 V2 变差”的结论有一部分是评估假象</h2>
<ul>
<li><b>数据缺失</b>：约 20% 的带标签序列在 8–19 位没有数据，位于 8–26 的 A1-F 会被错判成没扩出。现在只用 8–1510 都有数据的序列评估（带标签 5,046 条，留出集 2,196 条、348 个属；SILVA 8,434 条，Greengenes 2,523 条）。</li>
<li><b>比较规则很敏感</b>：“扩增失败算进去”的准确率取决于两条序列至少要有多少共同位点才可比。规则改变，S1 加 V1·V2 的差距从 −7 个百分点变到 0（见第五节）。所以早先页面里“每种方案 87.0%、90.8%、84.7%”之类的数字，只能当作粗略排序，不要当作精确差距。</li>
<li>修正后的结论是：V1·V2 加回去的收益在误差范围内，不是明显有害，但也不值得。</li>
</ul>

<h2>三、A1-R 扫描（起点 243–280，长度 17/19/21/23，设计集 Greengenes 一半 + SILVA 12,000 条）</h2>
{tR}
<p>位点覆盖率在 4 条以后基本到顶（6 条 96.7%，8 条 96.9%）。Tm（引物池内最低/最高）见第四节。</p>

<h2>四、A1-F 扫描，以及两端组合后的扩增子覆盖率</h2>
{tF}
<p>下面是 A1-F × A1-R 的组合，在“8–1510 都有数据”的序列上评估：</p>
{tA}
<p>3 条/22 和 4 条/30 的最低 Tm 只有 50 °C，4 条/64 约 53 °C，没有严重二聚体。</p>

<h2>五、加回 V1·V2 的属准确率（留出集 2,196 条，属级最近邻）</h2>
{tD}
<p>三列是同一份数据、三种“可比”规则；理想列是假设所有扩增子都扩出。差距随规则变化，说明扩增失败的处理方式对结论的影响比引物本身更大。</p>

<h2>六、局限</h2>
<ul>
<li>引物覆盖按“最多 1 个错配、3′ 端 3 个碱基匹配”算，真实 PCR 更复杂。没有检查对人基因组和线粒体的特异性。</li>
<li>属准确率只测到属，Greengenes 属标签偏向培养过的属；扩增失败的处理是我设计的简化。</li>
<li>设计集和留出集虽然分开，但都来自同两个数据库；GTDB 没评估。</li>
</ul>
<footer>脚本：<code>python_5R/explore/v1v2_scan_R.py</code>、<code>v1v2_scan_F.py</code>、<code>v1v2_eval.py</code>；数据：<code>docs/primer_design/v1v2_*.csv</code>、<code>robust_realistic.csv</code>。</footer>
</main></div>
"""
open(DOCS + "/primer_v1v2.html", "w", encoding="utf-8").write(head + body)
print("ok")
