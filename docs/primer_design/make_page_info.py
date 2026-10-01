"""由 info_*.csv / designs_info.json 生成 docs/primer_info_design.html"""
import re, json, html
import numpy as np, pandas as pd
HERE = __file__.rsplit("/", 1)[0]; DOCS = HERE.rsplit("/", 1)[0]
old = open(DOCS + "/5R_SMURF_pipeline.html", encoding="utf-8").read()
head = re.sub(r"<title>.*?</title>", "<title>5R 信息量优化设计</title>", old[:old.index("</style>") + 8], 1)
ID = pd.read_csv(f"{HERE}/info_design.csv"); CMP = pd.read_csv(f"{HERE}/info_compare.csv"); SITES = pd.read_csv(f"{HERE}/info_sites.csv")
V4 = pd.read_csv(f"{HERE}/info_v4.csv"); TM = pd.read_csv(f"{HERE}/info_single_tm.csv"); DIM = json.load(open(f"{HERE}/info_single_dimers.json"))
DI = json.load(open(f"{HERE}/designs_info.json"))
SITES4 = pd.read_csv(f"{HERE}/info_sites_v4.csv"); TM4 = pd.read_csv(f"{HERE}/info_single_tm_v4.csv"); DIM4 = json.load(open(f"{HERE}/info_single_dimers_v4.json"))
DI4 = json.load(open(f"{HERE}/designs_info_v4.json")); EV4 = pd.read_csv(f"{HERE}/info_eval_v4.csv"); PROF = pd.read_csv(f"{HERE}/profile_bact.csv", index_col=0)
NUMRE = re.compile(r"^[\d.%/ –→-]+$")
def cell(i, v, seqcol):
    cls = "seq" if i == seqcol else ("n" if NUMRE.match(str(v)) else "")
    return '<td class="' + cls + '">' + html.escape(str(v)) + "</td>"
def tbl(rows, heads, seqcol=None):
    h = "".join(f"<th>{c}</th>" for c in heads)
    b = "".join("<tr>" + "".join(cell(i, v, seqcol) for i, v in enumerate(r)) + "</tr>" for r in rows)
    return f'<div class="tbl"><table><thead><tr>{h}</tr></thead><tbody>{b}</tbody></table></div>'
pct = lambda x: f"{x * 100:.0f}%"

# 1. 不限条件下的最优（K=3..6），阈值 0.4/0.6
base = ID[ID.thr == "5R现有"].iloc[0]
free = ID[(ID.thr == "0.4/0.6") & (ID.read_len.isin([126, 150]))]
rows = [["5R 现有", "2×126", 5, f"{base.silva_bits_w:.0f}", f"{base.gg_bits_w:.0f}", f"{base.silva_regions:.1f}", base.worst_silva, base.amplicons]]
for _, r in free.iterrows():
    rows.append([f"优化 K≤{r.K}", f"2×{r.read_len}", r.K, f"{r.silva_bits_w:.0f}", f"{r.gg_bits_w:.0f}", f"{r.silva_regions:.1f}", r.worst_silva, r.amplicons])
t_free = tbl(rows, ["方案", "读长", "扩增子数", "SILVA 期望信息量 bit", "GG 留出 期望信息量 bit", "SILVA 平均扩出区域数", "各扩增子最差门覆盖（SILVA）%", "扩增子（E. coli 位置，括号为长度）"])

# 2. 推荐设计
amps = DI["amplicons"]
t_amp = tbl([[a["name"], a["lo"], a["hi"], a["len"], a["nbases"], f"{a['bits']:.0f}"] for a in amps], ["扩增子", "起点", "终点", "长度 bp", "测序碱基数", "信息量 bit"])
VREG = {"A1": "V1 · V2", "A2": "V3", "A3": "V5", "A4": "V6 · V7", "A5": "V9"}
t_amp = tbl([[a["name"], a["lo"], a["hi"], a["len"], VREG[a["name"]], a["nbases"], f"{a['bits']:.0f}"] for a in amps], ["扩增子", "起点", "终点", "长度 bp", "完整覆盖的可变区", "测序碱基数", "信息量 bit"])
tm = TM.set_index("site")
srows = []
for _, r in SITES.iterrows():
    if r.config == "single":
        srows.append([r.site, r.pos, r.seqs, tm.loc[r.site, "len"], tm.loc[r.site, "fold"], tm.loc[r.site, "GC"], tm.loc[r.site, "Tm"], pct(r["SILVA(验证)_mean"]), pct(r["SILVA(验证)_min"]), pct(r["GG留出_mean"]), pct(r["GG留出_min"])])
t_single = tbl(srows, ["引物", "E. coli 位置", "序列 5′→3′", "nt", "简并数", "GC", "Tm °C", "SILVA 主要门平均", "SILVA 最差门", "GG 主要门平均", "GG 最差门"], seqcol=2)
mrows = []
for _, r in SITES.iterrows():
    if r.config == "multi":
        mrows.append([r.site, r.pos, r.seqs, r.fold, pct(r["SILVA(验证)_mean"]), pct(r["SILVA(验证)_min"]), pct(r["GG留出_mean"]), pct(r["GG留出_min"])])
t_multi = tbl(mrows, ["引物", "E. coli 位置", "每位点 3 条（第 1 条即上表引物，| 分隔）", "总简并数", "SILVA 主要门平均", "SILVA 最差门", "GG 主要门平均", "GG 最差门"], seqcol=2)

# 3. 对照
crow = []
for _, r in CMP.iterrows():
    crow.append([r.db, r.scheme, f"{r.bits_w:.0f}", f"{r.bases_exp:.0f}", f"{r.regions:.2f}", pct(r.ge3), pct(r.all_), r.bases_all, f"{r.median_nn:.0f}", f"{r.nn_le2 * 100:.1f}%"])
t_cmp = tbl(crow, ["验证集", "方案", "期望信息量 bit", "期望测序碱基数", "平均扩出区域数", "≥3 个区域", "全部区域", "全部扩出时的测序碱基数", "最近邻距离中位数（错配数）", "有 ≤2 错配邻居的序列"])

# 4. V4
v4rows = []
for _, r in V4[V4.K.isin([4, 5])].iterrows():
    v4rows.append([f"2×{r.read_len}", r.K, r["mode"], f"{r.silva_bits:.0f}", f"{r.gg_bits:.0f}", r.amplicons, r.vregions.replace(":", " ").replace("V", "V")])
t_v4 = tbl(v4rows, ["读长", "K", "条件", "SILVA 期望信息量", "GG 期望信息量", "扩增子", "各可变区被测到的比例 %"])
def single_multi(SITES_, TM_, VREG_, DI_):
    tm_ = TM_.set_index("site")
    sr = [[r.site, r.pos, r.seqs, tm_.loc[r.site, "len"], tm_.loc[r.site, "fold"], tm_.loc[r.site, "GC"], tm_.loc[r.site, "Tm"], pct(r["SILVA(验证)_mean"]), pct(r["SILVA(验证)_min"]), pct(r["GG留出_mean"]), pct(r["GG留出_min"])] for _, r in SITES_.iterrows() if r.config == "single"]
    mr = [[r.site, r.pos, r.seqs, r.fold, pct(r["SILVA(验证)_mean"]), pct(r["SILVA(验证)_min"]), pct(r["GG留出_mean"]), pct(r["GG留出_min"])] for _, r in SITES_.iterrows() if r.config == "multi"]
    ta_ = tbl([[a["name"], a["lo"], a["hi"], a["len"], VREG_[a["name"]], a["nbases"], f"{a['bits']:.0f}"] for a in DI_["amplicons"]], ["扩增子", "起点", "终点", "长度 bp", "完整覆盖的可变区", "测序碱基数", "信息量 bit"])
    ts_ = tbl(sr, ["引物", "E. coli 位置", "序列 5′→3′", "nt", "简并数", "GC", "Tm °C", "SILVA 主要门平均", "SILVA 最差门", "GG 主要门平均", "GG 最差门"], seqcol=2)
    tm2_ = tbl(mr, ["引物", "E. coli 位置", "每位点 3 条（第 1 条即上表引物，| 分隔）", "总简并数", "SILVA 主要门平均", "SILVA 最差门", "GG 主要门平均", "GG 最差门"], seqcol=2)
    return ta_, ts_, tm2_

VREG4 = {"A1": "V1 · V2", "A2": "V4", "A3": "V6 · V7", "A4": "V8 · V9"}
t4a, t4s, t4m = single_multi(SITES4, TM4, VREG4, DI4)
dim4_txt = "；".join(f"{a} × {b}（{k} nt）" for a, b, k in DIM4)
VR_ = {"V1": (69, 99), "V2": (137, 242), "V3": (433, 497), "V4": (576, 682), "V5": (822, 879), "V6": (986, 1043), "V7": (1117, 1173), "V8": (1243, 1294), "V9": (1435, 1465)}
ent_rows = []
for v, (a, b) in VR_.items():
    e = PROF.entropy.reindex(range(a, b + 1)).fillna(0)
    ent_rows.append([v, b - a + 1, f"{e.sum():.1f}", f"{e.mean():.2f}"])
t_ent = tbl(sorted(ent_rows, key=lambda r: -float(r[3])), ["可变区", "长度 nt", "总熵 bit", "每个碱基平均熵 bit"])
ev4 = EV4.copy()
t_ev4 = tbl([[r.db, r.primers, f"{r.bits_w:.0f}", f"{r.regions:.2f}", pct(r.ge3), pct(r.all_), r.per, r.worst] for _, r in ev4.iterrows()],
            ["验证集", "每位点引物数", "期望信息量 bit", "平均扩出区域数", "≥3 个区域", "4 个全部", "各扩增子覆盖 %", "各扩增子最差门 %"])
dim_txt = "；".join(f"{a} × {b}（{k} nt）" for a, b, k in DIM)

body = f"""
<div class="wrap"><main style="grid-column:1/-1;width:100%;max-width:980px;margin-inline:auto">
<header>
  <div class="eyebrow">16S_5R · 引物设计</div>
  <h1>以信息量为目标：4–6 个扩增子、延伸到 16S 两端</h1>
  <p class="lede">仍是单管。目标改成“每条序列平均测到的信息量”：被测序碱基位置的保守度熵之和（bit），再乘以该序列被扩出的概率。用动态规划在“同管扩增子区间互不相交、间隔 ≥20 nt”的约束下求精确最优，允许边界延伸到 11 和 1509。</p>
  <p>我把“4–6 条引物”理解为 4–6 个扩增子（正反引物各一对）。这是计算机模拟，没有做实验；信息量是熵的代理指标，不等于种水平鉴定率。</p>
<p style="border:1px solid var(--warn);background:var(--surface);padding:10px 14px;border-radius:8px"><b>更正（以 <code>primer_v3v4_design.html</code> 为准）：</b>本页用“位点熵”作为信息量，这个指标和“能否区分属”相关性弱（Spearman 0.33）：熵最高的窗口（V1–V2）属准确率 78%，覆盖 V4 的窗口最高约 90%。因此本页“最优方案不含 V4”和“信息量 +69%”的结论不可靠；此外本页的多引物结果受一个已修复的算法缺陷影响（起始序列未按未覆盖序列加权）。</p>
  </header>

<h2>一、结论</h2>
<ul>
<li>相同的验证集（SILVA 4.8 万条里的近全长序列）上，现有 5R 的期望信息量是 <b>360 bit</b>。新设计 5 个扩增子、每位点 1 条引物是 <b>608 bit（+69%）</b>；每位点放 3 条引物是 <b>646 bit（+80%）</b>。Greengenes 留出集上是 387 → 643 → 689 bit。</li>
<li>如果只给现有 5R 的位点加到每位点 3 条引物，是 535 bit。所以新设计相对“迭代后的 5R”仍多 21%。</li>
<li><b>6 个扩增子不比 5 个更好</b>：5 个扩增子之间剩下的空隙最长只有 135 nt（1196–1331），放不下 ≥180 bp 的第 6 个；要放第 6 个就得压缩其他扩增子，期望信息量反而下降，优化器没有选它。4 个扩增子是 517 bit，5 个是 608 bit。</li>
<li>信息量增加主要来自三件事：①延伸到两端，加入 V1、V2、V9；②每个扩增子都接近 250 bp，两条 126 nt 的读段几乎不重叠，没有浪费碱基（现有 R4 只有 160 bp，重叠约 90 nt）；③引物更通用，平均每条序列能扩出的区域数从 2.8 增加到 3.6。</li>
<li><b>代价</b>：最优方案<b>没有 V4 和 V8</b>；两端的引物较弱（A1 最差门覆盖 44%，每位点放 3 条后提高到 62%）。</li>
</ul>

<h2>二、方法</h2>
<ul>
<li>熵：Greengenes 13_8 细菌序列在每个 E. coli 位置上的碱基分布。全长 16S 总熵约 1,098 bit，现有 5R 理想情况下覆盖约 671 bit（61%）。</li>
<li>一个扩增子的信息量：2×126 的读段（含引物）去掉引物后，被测到的位置的熵之和，重叠区只算一次。</li>
<li>期望信息量 = Σ（扩增子信息量 × 正反引物都匹配的概率），主要门等权。它对扩增子可加，所以可以用动态规划精确求解。</li>
<li>候选引物位点沿用前面的扫描（Greengenes 与 SILVA 两库的缺口比例都低；主要门平均覆盖率 ≥60%，最差门 ≥40%）。</li>
<li>设计集：Greengenes 85% OTU 集 + SILVA 1.2 万条里的近全长序列（位置 20–1510 都有数据）。验证集：SILVA 其余序列（10,443 条）和 Greengenes 留出集（3,099 条）。</li>
</ul>

<h2>三、不同扩增子数与读长下的最优</h2>
{t_free}
<p>读长 2×150 时扩增子可以放宽到 290 bp，信息量更高（5 个扩增子 665 bit，并且 V8 也被覆盖），但末端扩增子（1222–1509）的最差门覆盖只有 47%。</p>

<h2>四、推荐设计：5 个扩增子，2×126</h2>
{t_amp}
<p>引物位置（单引物版本）：</p>
{t_single}
<p>3′ 端互补 ≥5 nt 的引物对：{dim_txt}。Tm 范围 51–67 °C，不均匀，下单前需要重新平衡。</p>
<h3>每位点 3 条引物版本</h3>
<p>第 1 条是上表的引物，另外两条只为还没被覆盖的序列追加（各 ≤4 倍简并）。这样每个位点的覆盖不会比单引物差。</p>
{t_multi}

<h2>五、与现有 5R 的对照</h2>
{t_cmp}
<p>“最近邻距离”：在所有扩增子都扩出的序列里随机取最多 1,500 条，量每条序列与其他序列拼接签名的最小错配数。SMURF 把 ≤2 个错配当噪声，所以有 ≤2 错配邻居的序列无法被区分。各方案参与比较的序列不完全相同，只能作相对参考。</p>

<h2>六、V4 是高变区，为什么强制包含会降低信息量</h2>
<p>V4 确实是高变区，但在这个熵指标下它并不是信息最密的区域。每个碱基的平均熵：</p>
{t_ent}
<p>V4 总熵（103 bit）比 V3、V5 大，只是因为它长（107 nt）；按每个碱基算，V4 排在后面。强制 V4 之后信息量下降，原因是位置安排的机会成本，不是 V4 信息少：</p>
<ul>
<li><b>2×126</b>：V4 两侧的保守位点（515 和 785–806）相距 290 bp，放不进 250 bp。能放进去的是 514–709（196 bp），它自己贡献 103 bit，期望 104 bit。但它占用了原来 V3（338–576）和 V5（692–925）两个扩增子的位置（期望合计 236 bit），后面的扩增子被迫下移到 780–985（期望 87 bit）和 1050–1238（期望 75 bit），而自由最优的 947–1196 期望 133 bit。合计 607 → 504 bit。</li>
<li><b>2×150</b>：516–805（290 bp）可以被完整测到，所以 4 个扩增子的自由最优本来就包含 V4（620 bit），强制 V4 没有任何损失。之前我写的“降低约 8%”是拿它和 5 个扩增子且没有 V4 的 665 bit 比；含 V4 的方案放不下第 5 个扩增子，因为扩增子都是 286–290 bp。</li>
<li>熵是对整个种群逐位点的统计，不能反映“哪些位点能区分近缘种”。V4 在文献里常用于属水平分类，熵指标没有捕捉到这一点。如果你更看重 V4，这个指标会低估它的价值。</li>
</ul>
<p>两个读长下强制 V4 的结果：</p>
{t_v4}

<h2>七、包含 V4 的推荐版本（2×150，4 个扩增子）</h2>
<p>读长 2×150，扩增子 11–260 和 3 个 286–290 bp 的长扩增子，完整覆盖 V1、V2、V4、V6、V7、V8、V9（没有 V3 和 V5）。</p>
{t4a}
<p>期望信息量 614 bit（单引物）/ 649 bit（每位点 3 条），而现有 5R 是 360 bit，迭代后的 5R 是 535 bit。</p>
{t_ev4}
<p>引物（单引物版本）：</p>
{t4s}
<p>3′ 端互补 ≥5 nt 的引物对：{dim4_txt}。A2-F 的 Tm 高达 68–77 °C，与其他引物差 10 °C 以上。</p>
<h3>每位点 3 条引物版本</h3>
{t4m}
<ul>
<li>这版没有 V3 和 V5；V4 扩增子（516–805）和另外两个长扩增子都是 286–290 bp，在 FFPE 降解样本里扩增效率会下降，也可能在各区域 reads 占比里偏低（现有 5R 里 244 bp 的 R3 只有 4%–11%）。</li>
<li>首个扩增子（11–260）和末个扩增子（1222–1509）的引物较弱：最差门覆盖 44% 和 47%（单引物）/ 62% 和 71%（每位点 3 条）。</li>
<li>读长必须至少是 2×150；如果只有 2×126，290 bp 的扩增子中间会有 38 nt 测不到。</li>
</ul>

<h2>八、局限</h2>
<ul>
<li>熵是代理指标：位点保守度高低不等于能区分到种。真正检验要用已知组成的标准菌群或带种级分类的参考集。</li>
<li>数据库是 Greengenes 13_8 和 SILVA 128，不是最新版本；两端（1–100 位、1450 位以后）序列数据少，熵和覆盖率的估计不如中间可靠。</li>
<li>设计和验证的序列限定为近全长，可能偏向容易测通的类群。</li>
<li>首个扩增子的两条引物对主要门的最差覆盖只有 44%（单条）/ 62%（3 条），古菌几乎检不到。</li>
<li>5 个扩增子都接近 250 bp，按示例数据里“长扩增子占比低”的规律，各区域的 reads 可能比现有 5R 更均衡也可能更偏，没有实验数据，需要用标准菌群验证。</li>
<li>Tm 不均匀，部分引物有 3′ 端互补；没有检查对人基因组和线粒体的特异性。</li>
<li>k-mer 数据库必须用新引物重建；SMURF 需要把区域数从 5 的硬编码里解放出来（<code>scott_format_newer_func.m</code> 的 <code>nR = 5</code>）。</li>
</ul>
<footer>脚本：<code>python_5R/info_design.py</code>、<code>python_5R/iterate_5R.py</code>；数据：<code>docs/primer_design/info_*.csv</code>、<code>designs_info.json</code>。</footer>
</main></div>
"""
open(DOCS + "/primer_info_design.html", "w", encoding="utf-8").write(head + body)
print("ok")
