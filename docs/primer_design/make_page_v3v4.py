"""生成 docs/primer_v3v4_design.html"""
import re, json, html
import numpy as np, pandas as pd
HERE = __file__.rsplit("/", 1)[0]; DOCS = HERE.rsplit("/", 1)[0]
old = open(DOCS + "/5R_SMURF_pipeline.html", encoding="utf-8").read()
head = re.sub(r"<title>.*?</title>", "<title>V3+V4 单管引物设计</title>", old[:old.index("</style>") + 8], 1)
VR = pd.read_csv(f"{HERE}/vregion_acc.csv"); W = pd.read_csv(f"{HERE}/window_scan_clean.csv"); AD = pd.read_csv(f"{HERE}/acc_designs_clean.csv"); AF = pd.read_csv(f"{HERE}/acc_final5b.csv")
ST = pd.read_csv(f"{HERE}/final5b_sites.csv"); EV = pd.read_csv(f"{HERE}/final5b_eval.csv"); TM = pd.read_csv(f"{HERE}/final5b_tm.csv").set_index("site")
DIM = json.load(open(f"{HERE}/final5b_dimers.json")); BASE = pd.read_csv(f"{HERE}/base5R_valsets.csv")
NUMRE = re.compile(r"^[\d.%/ –→-]+$")
def cell(i, v, seqcol):
    cls = "seq" if i == seqcol else ("n" if NUMRE.match(str(v)) else "")
    return '<td class="' + cls + '">' + html.escape(str(v)) + "</td>"
def tbl(rows, heads, seqcol=None):
    h = "".join(f"<th>{c}</th>" for c in heads)
    b = "".join("<tr>" + "".join(cell(i, v, seqcol) for i, v in enumerate(r)) + "</tr>" for r in rows)
    return f'<div class="tbl"><table><thead><tr>{h}</tr></thead><tbody>{b}</tbody></table></div>'
pct = lambda x: f"{x * 100:.0f}%"

t_win = tbl([[f"{int(r.start)}–{int(r.end)}", f"{r.genus_acc * 100:.1f}%", f"±{r.se * 100:.1f}"] for _, r in W.iterrows()], ["250 nt 窗口（E. coli 位置）", "属水平最近邻准确率", "标准误 %"])
t_vr = tbl([[r.region, int(r.length), f"{r.genus_acc * 100:.1f}%", f"{r.bits_per_nt:.2f}"] for _, r in VR.sort_values("genus_acc", ascending=False).iterrows()], ["区域（只用该区域自身的碱基）", "长度 nt", "属水平最近邻准确率", "每个碱基平均熵 bit"])
t_accd = tbl([[r.design, int(r.bases), f"{r.genus_acc * 100:.1f}%", f"±{r.genus_se * 100:.1f}", f"{r.family_acc * 100:.1f}%"] for _, r in AD.iterrows()], ["方案（假设所有扩增子都扩出）", "测序碱基数", "属准确率", "标准误 %", "科准确率"])
t_accf = tbl([[r.design, int(r.bases), f"{r.genus_acc * 100:.1f}%", f"±{r.genus_se * 100:.1f}", f"{r.family_acc * 100:.1f}%"] for _, r in AF.iterrows()], ["方案（假设所有扩增子都扩出）", "测序碱基数", "属准确率", "标准误 %", "科准确率"])

amps = [("A1", "11–260", 250, "V1 · V2"), ("A2", "337–530", 194, "V3"), ("A3", "558–802", 245, "V4"), ("A4", "910–1195", 286, "V6 · V7"), ("A5", "1222–1509", 288, "V8 · V9")]
t_amp = tbl([[a, p, l, v] for a, p, l, v in amps], ["扩增子", "E. coli 位置", "长度 bp", "完整覆盖的可变区"])
rows = []
for _, r in ST.iterrows():
    if r.config == "single":
        t = TM.loc[r.site]
        rows.append([r.site, r.pos, r.seqs, t.len, t.fold, t.GC, t.Tm, pct(r["SILVA(验证)_mean"]), f"{pct(r['SILVA(验证)_min'])}（{r['SILVA(验证)_worst']}）", pct(r["GG干净(验证)_mean"]), f"{pct(r['GG干净(验证)_min'])}（{r['GG干净(验证)_worst']}）"])
t_single = tbl(rows, ["引物", "位置", "序列 5′→3′", "nt", "简并数", "GC", "Tm °C", "SILVA 主要门平均", "SILVA 最差门", "GG 主要门平均", "GG 最差门"], seqcol=2)
rows = []
for _, r in ST.iterrows():
    if r.config == "multi":
        rows.append([r.site, r.pos, r.seqs, r.fold, pct(r["SILVA(验证)_mean"]), f"{pct(r['SILVA(验证)_min'])}（{r['SILVA(验证)_worst']}）", pct(r["GG干净(验证)_mean"]), f"{pct(r['GG干净(验证)_min'])}（{r['GG干净(验证)_worst']}）"])
t_multi = tbl(rows, ["引物", "位置", "每位点多条引物（| 分隔；第 1 条即上表引物）", "总简并数", "SILVA 主要门平均", "SILVA 最差门", "GG 主要门平均", "GG 最差门"], seqcol=2)
cmp_rows = []
for _, r in BASE.iterrows():
    cmp_rows.append([r.db, r.scheme, r.per, r.worst, f"{r.regions:.2f}", pct(r.ge3), pct(r.all5)])
for _, r in EV.iterrows():
    cmp_rows.append([r.db, "新方案 5 扩增子，" + r.primers, r.per, r.worst, f"{r.regions:.2f}", pct(r.ge3), pct(r.all5)])
cmp_rows.sort(key=lambda x: (x[0], 0 if x[1].startswith("5R 现有") else (1 if x[1].startswith("5R 位点") else 2)))
t_cmp = tbl(cmp_rows, ["验证集", "方案", "各扩增子覆盖 %", "各扩增子最差门 %", "平均扩出区域数", "≥3 个区域", "全部区域"])
dim_txt = "；".join(f"{a} × {b}（{k} nt）" for a, b, k in DIM)

body = f"""
<div class="wrap"><main style="grid-column:1/-1;width:100%;max-width:980px;margin-inline:auto">
<header>
  <div class="eyebrow">16S_5R · 引物设计</div>
  <h1>V3 + V4 单管设计，扩增子上限 290 bp</h1>
  <p class="lede">回答三个问题：为什么常规方案用 V4、信息密度怎么算才合理、V3 和 V4 能不能放进同一个管。最后给出一套 5 个扩增子、同时包含 V3 和 V4、每个都不超过 290 bp 的单管方案。</p>
  <p>这是计算机模拟，没有做实验。数据库是 Greengenes 13_8 和 SILVA 128，不是最新版本；GTDB 本环境下载不了，没有评估。</p>
</header>

<h2>更正：我之前的做法有三处错误</h2>
<ul>
<li><b>“信息量”用熵来算不合适。</b>熵是对整个种群逐位点的统计，和“能不能区分属”不是一回事。在下面的窗口扫描里，熵最高的窗口（V1–V2）属准确率只有 78%，而覆盖 V4 的窗口最高，约 90%。我之前基于熵的设计会把 V4 排除掉，是这个原因。</li>
<li><b>556–557 位附近的引物评估在两个数据库里不一致，我之前的解释不完整。</b>我曾把原因归为参照序列里的一个 <code>N</code>；把 <code>N</code> 修正后，这个位点仍然有问题：Greengenes 里 44% 的序列在起点 554–557 的 18 nt 窗口里有 1 个碱基的缺口（起点 ≥558 只有 1%），同一条引物在 SILVA 上主要门平均 91%，在 Greengenes 上只有 35%。这个位点附近有长度多态性，两个库的比对对它的处理方式也不同，我没法确定哪个更接近真实扩增。起点 ≥558 的引物不跨过这个位置，在两个库里的结果一致，所以最终方案把 V4 正向引物放在 558。</li>
<li><b>“每位点追加多条引物”的算法有缺陷。</b>追加引物时，起始序列用的是全体序列的多数碱基，没有跟着“尚未被覆盖的序列”走，所以没法为某个类群设计专用引物。修复后，之前给 5R 的多引物结果略有提高（每位点 3 条：主要门平均 91.5% → 92.4%，最差门 67% → 74%），页面 <code>primer_5R_iteration.html</code> 已更新。</li>
</ul>

<h2>一、为什么常规单扩增子方案用 V4</h2>
<p>我用 3,908 条带属标签的近全长细菌序列（482 个属，Greengenes 97% 代表序列，用修正后的参照做的干净投影），对每个 250 nt 窗口做留一最近邻属分类：窗口里的序列碱基完全按错配数找最近的另一条序列，看它是不是同一个属。</p>
{t_win}
<ul>
<li><b>信息密度最高的窗口是覆盖 V4 的窗口</b>：460–709、510–759、560–809、610–859 的属准确率都在 89.5%–90.4%。V1–V2（60–309）是 78.4%，V3 窗口（310–559）是 79.6%，V6–V8（860–1209）是 76%–78%。</li>
<li>熵和准确率的相关性弱：窗口的熵与属准确率的 Spearman 相关系数只有 0.33。熵最高的窗口（60–309，239 bit）准确率只有 78%，准确率最高的窗口熵只有 184–195 bit。</li>
</ul>
<p>只取每个可变区自身的碱基（不含两侧）：</p>
{t_vr}
<ul>
<li><b>V4 单独 107 nt 就有 87.2% 的属准确率，长度相近的 V2（106 nt）只有 71.7%，而 V4 的每碱基熵（0.97 bit）在 9 个可变区里排最后。</b>熵衡量的是“变化多”，包括同一个属内部的变化；区分属需要的是“属之间不同、属内部一致”的变化，V4 的变异更符合这一点。这是对数据的解释，没有用其他数据集验证。</li>
<li>V1（31 nt）和 V9（31 nt）最短，准确率只有 44%–47%，说明两端的可变区单独用信息很少，主要靠和其他区域组合。</li>
<li>引物通用性是另一个原因：在 Greengenes 85% 集上（门等权），515F-Y 和 806R 的覆盖率分别约 92% 和 97%，对古菌也有约 97% 的覆盖；926R、1392R 也有类似的广谱性。其他原因（V4 扩增子约 290 bp，适合 2×150 或 2×250；参考数据库积累最多；EMP 等大型项目标准化）是常见说法，本次没有验证。</li>
</ul>
<p><b>这个计算是否合理？</b>比熵合理，但仍有局限：①只测了“能不能区分到属”，没有测种；②Greengenes 的属标签稀疏，只有有标签的 3,908 条参与，偏向培养过的属；③最近邻按序列错配数找，和 SMURF 的实际算法不同；④标准误约 0.5%，窗口之间差 1% 以内不能当作有差别。</p>

<h2>二、V3 和 V4 能不能兼顾</h2>
<p>你问的是“是不是不存在”：<b>没有一个扩增子能同时完整包含 V3 和 V4（≤290 bp，引物质量合格），但把它们放成两个相邻扩增子是可行的。</b></p>
<ul>
<li><b>单个扩增子</b>：V3 从 433 开始，V4 到 682 结束，中间的内部碱基至少 250 nt，加上两条 18 nt 引物，最短 285 bp。≤290 bp 的话，两条引物必须紧贴着 V3 起点（415–432）和 V4 终点（683–704），这两个位点不保守：搜索出的最好组合，最差主要门覆盖只有 3%–9%（Bacteroidetes 和 Spirochaetes 几乎扩不出）。</li>
<li><b>两个相邻扩增子</b>：V3 扩增子 337–530（194 bp）+ V4 扩增子 558–802（245 bp），两个扩增子区间不相交、间隔 27 nt。两个数据库各自独立搜索，得到的都是约 337–530 + 556–802 的组合；我把 V4 正向引物后移到 558，避开上面说的长度多态位点。</li>
<li>V4 正向引物（558–575）是关键：多数门是 <code>GGGCGTAAA</code>，Bacteroidetes 是 <code>GGGTTTAAA</code>，用一条引物最多容许 1 个错配，Bacteroidetes 只有 3%–8%。放 2 条引物（一条通用、一条 Bacteroidetes 专用）后，所有 11–12 个主要门在两个库的验证集上都 ≥79%。</li>
</ul>
<p>在干净的带标签集上，用两个扩增子的属准确率：</p>
{t_accd}

<h2>三、最终方案：5 个扩增子，单管，上限 290 bp，读长 2×150</h2>
{t_amp}
<p>覆盖 V1、V2、V3、V4、V6、V7、V8、V9，没有 V5。长度 194–288 bp，所以必须用 2×150：A4 和 A5（286、288 bp）在 2×126 下中间会有 34–36 nt 测不到。</p>
<p>属分类准确率（假设所有扩增子都扩出，干净带标签集，3,908 条）：</p>
{t_accf}
<ul>
<li>最终方案 94.4%，现有 5R 93.4%，标准误约 0.4%，差 1 个百分点，只有大约 2 个标准误，不要过度解读。</li>
<li>去掉 V3 变成 93.6%，去掉 V4 变成 91.6%：V4 贡献的区分力比 V3 大。</li>
<li>之前我用熵算“信息量”，得到新设计比 5R 多 69% 的结论；用属准确率（假设所有扩增子都扩出）看，差别只有约 1 个百分点。</li>
</ul>

<h3>扩增覆盖，对比现有 5R</h3>
<p>验证集：SILVA 其余序列的近全长部分（10,443 条）和没参与设计的干净 Greengenes 一半（3,151 条）。设计只用了另一半 Greengenes 和 SILVA 的 1.2 万条。</p>
{t_cmp}
<p>现有 5R 在 R1–R5 上各扩增子的最差门覆盖只有 0%–5%，每个扩增子都有至少一个门几乎扩不出来。新方案每位点放多条引物后，每个扩增子的最差门覆盖是 67%–89%；只放 1 条时，V4 扩增子（A3）的最差门只有 2%–8%（Bacteroidetes）。</p>

<h3>引物（单引物版本）</h3>
{t_single}
<p>3′ 端互补 ≥5 nt 的引物对：{dim_txt}。</p>
<h3>每位点多条引物版本（推荐）</h3>
<p>第 1 条是上表的引物，其余只为还没被覆盖的序列追加（各 ≤4 倍简并）。</p>
{t_multi}

<h2>四、局限和风险</h2>
<ul>
<li><b>扩增子长度不均匀（194–288 bp）。</b>在示例数据里，长扩增子在管里占比明显低（现有 5R 里 244 bp 的 R3 只有 4%–11%），这个方案的 A4、A5 接近 290 bp，很可能占比偏低，FFPE 降解样本更严重。需要用标准菌群调引物浓度。</li>
<li><b>总简并数高。</b>每位点多条引物版本一个位点最多 16 种寡核苷酸，每种有效浓度更低，非特异性风险也更高，没有检查对人基因组和线粒体的特异性。</li>
<li><b>Tm 和二聚体。</b>各引物 Tm 差异较大（见表），个别引物 3′ 端有互补，需要在 Primer3 或 OligoAnalyzer 里复核。</li>
<li><b>首尾引物弱。</b>A1 和 A5 引物对某些门（如 Bacteroidetes、Tenericutes）覆盖偏低。</li>
<li><b>数据库不是最新，且分类信息稀疏。</b>属准确率只测了有属标签的类群；SILVA 的门标签是我推断的。</li>
<li>k-mer 数据库必须用新引物重建（<code>smurf5r.build_region_db_from_fasta</code> 已支持每位点多条引物）；SMURF 的区域数要从 5 的硬编码里解放（<code>nR = 5</code>）。</li>
<li>本次没有重算 5 个方案之外的结果；<code>primer_info_design.html</code>（熵目标）和 <code>primer_redesign.html</code>（两库设计）里的结论请以本页为准。</li>
</ul>
<footer>脚本：<code>python_5R/primer_design.py</code>、<code>info_design.py</code>、<code>iterate_5R.py</code>；数据：<code>docs/primer_design/</code>（window_scan_clean.csv、acc_*.csv、final5b_*.csv）。</footer>
</main></div>
"""
open(DOCS + "/primer_v3v4_design.html", "w", encoding="utf-8").write(head + body)
print("ok")
