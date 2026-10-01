"""由 designs_single.json 等数据生成 docs/primer_redesign.html（单管设计；用 Greengenes 13_8 与 SILVA 128 交叉验证）"""
import json, re, html
import numpy as np, pandas as pd

HERE = __file__.rsplit("/", 1)[0]
DOCS = HERE.rsplit("/", 1)[0]
OUT = DOCS + "/primer_redesign.html"
old = open(DOCS + "/5R_SMURF_pipeline.html", encoding="utf-8").read()
head = re.sub(r"<title>.*?</title>", "<title>5R 引物重设计</title>", old[:old.index("</style>") + 8], 1)
D = json.load(open(f"{HERE}/designs_single.json"))
prof = pd.read_csv(f"{HERE}/profile_bact.csv", index_col=0)
RB = pd.read_csv(f"{HERE}/region_bias.csv")
EV = pd.read_csv(f"{HERE}/eval_silva_fixed.csv")
HO = pd.read_csv(f"{HERE}/holdout.csv")
VR = {"V1": (69, 99), "V2": (137, 242), "V3": (433, 497), "V4": (576, 682), "V5": (822, 879),
      "V6": (986, 1043), "V7": (1117, 1173), "V8": (1243, 1294), "V9": (1435, 1465)}
FIVE = [("R1", 103, 332), ("R2", 338, 536), ("R3", 685, 927), ("R4", 944, 1104), ("R5", 1175, 1391)]
NAMES = {"S1_V4full_V6_V8": "S1", "S2_V3_V6V7_le250": "S2", "S3_V4full_V6V7_V9": "S3", "S4_V4short_V6V7_V9": "S4"}
pct = lambda x: f"{x * 100:.0f}%"


def track_svg():
    X0, X1 = 105, 885
    sx = lambda p: X0 + (p - 1) / 1541 * (X1 - X0)
    s = prof.top_freq.reindex(range(1, 1543)).fillna(0).rolling(9, center=True, min_periods=1).mean()
    pts = " ".join(f"{sx(p):.1f},{150 - 120 * (v - 0.4) / 0.6:.1f}" for p, v in s.items())
    rows = [("5R 现有", [(n, a, b) for n, a, b in FIVE])]
    for k, lab in [("S1_V4full_V6_V8", "方案 S1"), ("S2_V3_V6V7_le250", "方案 S2"), ("S3_V4full_V6V7_V9", "方案 S3"), ("S4_V4short_V6V7_V9", "方案 S4")]:
        rows.append((lab, [(a["amplicon"], a["start"], a["end"]) for a in D[k]["amplicons"]]))
    o = ['<svg viewBox="0 0 900 355" width="100%" style="min-width:680px" role="img" aria-label="16S 保守性与单管引物平铺">']
    for v, (a, b) in VR.items():
        o.append(f'<rect x="{sx(a):.1f}" y="20" width="{sx(b) - sx(a):.1f}" height="140" fill="var(--accent-soft)"/>'
                 f'<text x="{(sx(a) + sx(b)) / 2:.1f}" y="14" font-size="11" text-anchor="middle" fill="var(--muted)">{v}</text>')
    for y, lab in [(150, "40%"), (90, "70%"), (30, "100%")]:
        o.append(f'<line x1="{X0}" x2="{X1}" y1="{y}" y2="{y}" stroke="var(--line)"/><text x="{X0 - 6}" y="{y + 4}" font-size="10" text-anchor="end" fill="var(--muted)">{lab}</text>')
    o.append(f'<polyline points="{pts}" fill="none" stroke="var(--ink)" stroke-width="1.2"/>')
    o.append(f'<text x="{X0}" y="178" font-size="11" fill="var(--muted)">每个位点最常见碱基的频率（Greengenes 细菌，9 nt 滑动平均）</text>')
    y = 196
    for lab, items in rows:
        o.append(f'<text x="{X0 - 6}" y="{y + 11}" font-size="11" text-anchor="end" fill="var(--ink)">{lab}</text>')
        for n, a, b in items:
            o.append(f'<rect x="{sx(a):.1f}" y="{y}" width="{max(sx(b) - sx(a), 2):.1f}" height="14" rx="2" fill="var(--accent)"/>'
                     f'<text x="{sx(a) + 3:.1f}" y="{y + 11}" font-size="10" fill="var(--surface)">{n}</text>')
        y += 24
    o.append(f'<text x="{X0}" y="{y + 16}" font-size="11" fill="var(--muted)">E. coli 16S 位置 1–1542。所有方案都是单管，扩增子区间互不相交。</text></svg>')
    return "\n".join(o)


def tbl(df, cols, heads, num=()):
    h = "".join(f"<th>{c}</th>" for c in heads)
    b = "".join("<tr>" + "".join(f'<td class="{"seq" if c == "primer" else ("n" if c in num else "")}">{html.escape(str(r[c]))}</td>' for c in cols) + "</tr>" for _, r in df.iterrows())
    return f'<div class="tbl"><table><thead><tr>{h}</tr></thead><tbody>{b}</tbody></table></div>'


def plan(name):
    P = pd.DataFrame(D[name]["primers"]); A = pd.DataFrame(D[name]["amplicons"])
    for c in ["GG_key_mean", "GG_key_min", "SILVA_key_mean", "SILVA_key_min", "SILVA_arch"]:
        P[c] = P[c].map(pct)
    P["SILVA_min"] = P.SILVA_key_min + "（" + P.SILVA_worst + "）"
    ta = tbl(A, ["amplicon", "start", "end", "length", "regions", "F_site", "R_site"],
             ["扩增子", "起点", "终点", "长度 bp", "完整覆盖", "正向位点", "反向位点"], num=("start", "end", "length"))
    tp = tbl(P, ["name", "site", "primer", "length", "fold", "GC", "Tm", "GG_key_mean", "SILVA_key_mean", "SILVA_min", "SILVA_arch"],
             ["引物", "E. coli 位置", "序列 5′→3′", "nt", "简并数", "GC", "Tm °C", "GG 主要门平均", "SILVA 主要门平均", "SILVA 最差主要门", "SILVA 古菌"],
             num=("length", "fold", "GC", "Tm", "GG_key_mean", "SILVA_key_mean", "SILVA_min", "SILVA_arch"))
    d = D[name]["dimers"]
    dt = "；".join(f"{x['pair']}（{x['nt']} nt）" for x in d) if d else "无"
    return ta + tp + f"<p>3′ 端互补 ≥5 nt 的引物对：{dt}。</p>"


def summary_row(lab, sel):
    g = lambda c: np.mean([x for x in sel[c]])
    return f"<tr><td>{lab}</td>" + "".join(f'<td class="n">{pct(v)}</td>' for v in sel) + "</tr>"


five_sv = EV[(EV.set == "5R current") & (EV.db == "SILVA128")]
five_gg = EV[(EV.set == "5R current") & (EV.db == "GG")]
five_ho = HO[HO.set == "5R current"]
rows_cmp = []
vals = {"5R 现有": (five_gg.key_mean.mean(), five_gg.key_min.mean(), five_ho.key_mean.mean(), five_ho.key_min.mean(), five_sv.key_mean.mean(), five_sv.key_min.mean(), five_sv.arch.mean())}
for k, n in NAMES.items():
    P = pd.DataFrame(D[k]["primers"])
    vals[n] = (P.GG_key_mean.mean(), P.GG_key_min.mean(), P.GGholdout_key_mean.mean(), P.GGholdout_key_min.mean(),
               P.SILVA_key_mean.mean(), P.SILVA_key_min.mean(), P.SILVA_arch.mean())
cmp_html = "".join(f"<tr><td>{n}</td>" + "".join(f'<td class="n">{v * 100:.1f}%</td>' for v in vs) + "</tr>" for n, vs in vals.items())
regs = {"S4": "V4 · V6+V7 · V9", "S1": "V4 · V6 · V8", "S2": "V3 · V6+V7", "S3": "V4 · V6+V7 · V9", "5R 现有": "V2 · V3 · V5 · V6 · V8"}

rb = RB.copy()
rb_tab = rb.pivot(index="region", columns="sample", values="observed_share")
rb_exp = rb.pivot(index="region", columns="sample", values="expected_equal")
L = {"R1": 232, "R2": 193, "R3": 244, "R4": 160, "R5": 217}
bias_rows = "".join(
    f'<tr><td>{r}（{L[r]} bp）</td>' + "".join(f'<td class="n">{rb_tab.loc[r, s] * 100:.1f}%</td>' for s in rb_tab.columns) +
    "".join(f'<td class="n">{rb_exp.loc[r, s] * 100:.0f}%</td>' for s in rb_exp.columns) + "</tr>" for r in rb_tab.index)
cols = list(rb_tab.columns)

body = f"""
<div class="wrap"><main style="grid-column:1/-1;width:100%;max-width:980px;margin-inline:auto">
<header>
  <div class="eyebrow">16S_5R · 引物设计</div>
  <h1>5R 引物重设计：单管、加入 V4、用 SILVA 交叉验证</h1>
  <p class="lede">和 5R 一样，所有引物放在同一管里扩增。在“扩增子 180–250 bp、引物位于保守区、同管引物位点互不相交”的约束下，用 Greengenes 13_8 设计候选，再用 SILVA 128 的 6 万条序列和 Greengenes 的 1.2 万条留出序列做交叉验证。</p>
  <p>这是计算机模拟的设计，没有做实验验证。SILVA 的门标签是我用 Greengenes 训练的分类器推断的，不是 SILVA 官方标签；GTDB 在本环境下载不了，没有评估。</p>
</header>

<h2>一、现有 5R：引物位置与示例数据里各区域的 reads</h2>
<p>10 条引物按 <i>E. coli</i> 坐标：R1 103–332、R2 338–536、R3 685–927、R4 944–1104、R5 1175–1391。相邻扩增子的引物位点最近间隔 5 nt（R1/R2）和 16 nt（R3/R4），没有重叠，所以不会出现短产物抢引物。</p>
<p>但示例数据里 5 个区域的 reads 差得很大（进入模型的 reads 占比）：</p>
<div class="tbl"><table><thead><tr><th>区域（扩增子长度）</th>{"".join(f"<th>{c} 实测</th>" for c in cols)}{"".join(f"<th>{c} 等效率期望</th>" for c in cols)}</tr></thead><tbody>{bias_rows}</tbody></table></div>
<ul>
<li>最多的区域与最少的区域相差 <b>3.7 倍</b>（RDB123）和 <b>12.2 倍</b>（RDB1）。如果各区域扩增效率相同，重建出的菌群应使每个区域各占约 19–22%。</li>
<li>偏倚与<b>扩增子长度</b>明显相关：reads 占比与长度的 Spearman 相关系数是 −0.9 和 −0.7（只有 5 个区域、2 个样本，只能说方向一致，不能证明因果）。最长的 R3（244 bp）在 RDB1 里只有 4%，最短的 R2、R4 反而占得多。</li>
<li>引物 0 错配的比例解释不了它（Spearman 约 0.3）：R5 在 RDB123 里 0 错配的菌群只有 1%，reads 占比却有 23%。引物浓度、GC 含量、测序簇生成对短片段的偏好也可能起作用，这里没法区分。</li>
<li>SMURF 模型假设各区域效率相同，所以这个偏倚会直接进入丰度估计。设计新引物时应把扩增子长度尽量拉近，并留出按引物对调浓度的余地。</li>
</ul>

<h2>二、数据来源和方法</h2>
<ul>
<li><b>Greengenes 13_8</b>：4,797 条 85% OTU 代表序列做设计，另取 12,000 条不在其中的 97% OTU 序列做留出验证。</li>
<li><b>SILVA 128</b>：SEPP 参考包 <code>sepp-refsilva128</code>（bioconda）里的 99% OTU 全长比对，共 395,440 条，随机取 60,000 条（其中古菌 2,815 条），直接用它自己的比对列换算到 <i>E. coli</i> 坐标。门由 Greengenes 97% 代表序列训练的 8-mer 朴素贝叶斯分类器推断，99.96% 的序列后验概率 &gt; 0.9。</li>
<li><b>覆盖率</b>：最多 1 个错配、3′ 端 3 个碱基必须匹配。主要门指 12 个（Firmicutes、Proteobacteria、Bacteroidetes、Actinobacteria、Fusobacteria、Verrucomicrobia、Tenericutes、Spirochaetes、Cyanobacteria、Chloroflexi、Acidobacteria、Planctomycetes），每个门权重相同。</li>
<li><b>候选位点</b>：18–20 nt，最多 3 个简并位置、总简并数 ≤ 8，在两个数据库里各自取最优引物，取覆盖率较低者作为保守估计，并要求位点处两库的缺口比例都很低（见第六节）。</li>
<li><b>单管约束</b>：扩增子 180–250 bp（含引物；要完整包含 V4 时放宽到 292 bp），任意两个扩增子的区间不相交且间隔 ≥ 20 nt。用动态规划求完整覆盖可变区最多的组合。Tm 把 5′ 端调到 17–24 nt，靠近 62 °C。</li>
</ul>

<h2>三、保守性与单管能放进去的东西</h2>
<figure class="figure">{track_svg()}
<figcaption>灰蓝色底纹为可变区。引物必须落在曲线接近 100% 的位置。</figcaption></figure>
<p>两个数据库一致的结论：</p>
<ul>
<li><b>要在单管里完整放进 V4，有两条路</b>：①正向引物在 515 附近（经典 515F 位点）、反向引物在 785–804，扩增子 290 bp，引物质量好（每条引物对最差主要门的覆盖率约 86%–96%）；②正向引物在 520 附近、反向引物在 692–709，扩增子只有 190 bp，但反向引物落在保守性差的位置，对 Tenericutes 等门覆盖只有约 70%（SILVA）、62%（Greengenes）。我把候选位点的最低覆盖率要求从 80% 一路放宽到 30% 搜索，除了这两条没有找到更好的单个 V4 扩增子；第二条要把要求降到 60% 才会出现。</li>
<li>把 V4 压到 250 bp 的另一种做法（正向引物在 556 附近）只在 Greengenes 上看起来可行，在 SILVA 上这条引物几乎不匹配（主要门平均只有 23%），原因见第七节，不要采用。</li>
<li>所以单管下要在“引物质量”和“扩增子长度均衡”之间取舍：S1 用 290 bp 的 V4 换来好引物；S4 用 190 bp 的 V4，三个扩增子长度都在 180–226 bp，但 V4 反向引物和 V9 正向引物较弱。</li>
</ul>

<h2>四、单管方案与现有 5R 的比较</h2>
<div class="tbl"><table><thead><tr><th>方案</th><th>完整覆盖的可变区</th><th>扩增子</th><th>GG 设计集 主要门平均</th><th>GG 设计集 最差门</th><th>GG 留出集 主要门平均</th><th>GG 留出集 最差门</th><th>SILVA 主要门平均</th><th>SILVA 最差门</th><th>SILVA 古菌</th></tr></thead><tbody>
<tr><td>5R 现有</td><td>{regs["5R 现有"]}</td><td class="n">5</td>{"".join(f'<td class="n">{v * 100:.1f}%</td>' for v in vals["5R 现有"])}</tr>
<tr><td><b>S1（推荐）</b></td><td>{regs["S1"]}</td><td class="n">3（290/178/227）</td>{"".join(f'<td class="n">{v * 100:.1f}%</td>' for v in vals["S1"])}</tr>
<tr><td>S2</td><td>{regs["S2"]}</td><td class="n">2（190/226）</td>{"".join(f'<td class="n">{v * 100:.1f}%</td>' for v in vals["S2"])}</tr>
<tr><td>S3</td><td>{regs["S3"]}</td><td class="n">3（290/226/180）</td>{"".join(f'<td class="n">{v * 100:.1f}%</td>' for v in vals["S3"])}</tr>
<tr><td>S4</td><td>{regs["S4"]}</td><td class="n">3（190/226/180）</td>{"".join(f'<td class="n">{v * 100:.1f}%</td>' for v in vals["S4"])}</tr>
</tbody></table></div>
<p>每个数字是该方案全部引物的平均值；“最差门”是每条引物在 12 个主要门里最低覆盖率的平均。S1 的数字在 Greengenes 设计集、Greengenes 留出集和 SILVA 上一致（94.7%、97.8%、97.3%），说明不是对设计集过拟合。古菌覆盖率只在 SILVA 评估（2,815 条）：S1 的 V4 引物对约 55%–98%，V6 对 62%–91%，V8 正向引物为 0%。</p>
<p>与现有 5R 相比，S1 的得失：<b>得到</b> V4，引物更通用；<b>失去</b> V2、V3、V5（现有 5R 完整覆盖 V2、V3、V5、V6、V8）。扩增子从 5 个减到 3 个。S1 的 V4 扩增子（290 bp）是所有方案里最长的，长度偏倚的影响见第五节。</p>

<h2>五、推荐方案 S1：V4 + V6 + V8</h2>
{plan("S1_V4full_V6_V8")}
<ul>
<li><b>V4 扩增子 290 bp 比其他扩增子长很多</b>（178 bp 和 227 bp）。按第一节的规律，长扩增子在同一管里会占比偏低（现有 5R 里 244 bp 的 R3 只有 4%–11%），FFPE 降解样本更严重。建议用已知组成的标准菌群（mock community）把三对引物的浓度做配比，让长扩增子的引物浓度更高。</li>
<li>同管里最近的跨区域产物是 F1+R2（515–1082，568 bp）和 F2+R3（905–1406，502 bp），都很长，在降解样本里很难扩出来。</li>
<li>A2-F 的 GC 只有 28%、Tm 47–53 °C，明显低于其他引物，A1-F 的 Tm 又高达 69–74 °C。这是因为 3′ 端必须落在保守位置，只能调 5′ 端。下单前需要在 Primer3 或 IDT 工具里重新做 Tm 平衡，并微调位点 1–3 nt。</li>
<li>A1-F 与自身有 8 nt 的 3′ 互补（简并展开版本之间），A3-F 与自身有 6 nt。需要用 OligoAnalyzer 核对实际最差组合。</li>
</ul>

<h2>六、备选方案</h2>
<p><b>怎么选</b>：样本 DNA 较完整、能用 mock community 做引物配比，选 S1；样本降解严重（如 FFPE），或者更看重各区域 reads 均衡，选 S4（扩增子 180–226 bp）或 S2（没有 V4）。S4 的短 V4 反向引物和 V9 正向引物覆盖较差，实验前要重点验证。</p>
<h3>S2：所有扩增子 ≤250 bp，没有 V4</h3>
<p>适合降解严重的样本。两个扩增子长度接近（190 和 226 bp），引物覆盖略低于 S1。</p>
{plan("S2_V3_V6V7_le250")}
<h3>S4：短 V4 + V6/V7 + V9（三个扩增子都 ≤226 bp）</h3>
<p>V4 扩增子 190 bp，只覆盖 V4 本身，不含两侧。A1-R（692–709）GC 只有 36%、Tm 50–57 °C，对主要门最低覆盖 62%（GG）/ 70%（SILVA），古菌几乎检不到；A3-F（V9）对主要门最低 61%（GG）/ 66%（SILVA）。</p>
{plan("S4_V4short_V6V7_V9")}
<h3>S3：V4 + V6/V7 + V9</h3>
<p>在 S1 的基础上用 V6+V7 的长扩增子（226 bp）和 V9（180 bp）。V9 正向引物 A3-F 对主要门最低只有 66%（SILVA），GG 留出集 63%，覆盖明显不如其他引物，古菌覆盖也低。V9 在 Greengenes 里序列被大量截短，Greengenes 的估计不可靠。</p>
{plan("S3_V4full_V6V7_V9")}

<h2>七、这次验证里发现的两个问题</h2>
<ul>
<li><b>Greengenes 投影的伪影</b>：我最初在 Greengenes 上还找到过一个 237 bp 的 V4 组合（555–574 正向引物），在 Greengenes 上主要门最低覆盖率 66%。换到 SILVA 上，这条引物的主要门平均只有 23%。原因是我的参照序列在 558 位有一个简并碱基 <code>N</code>，导致 Greengenes 里 44% 的序列在该位置被比对软件插入假缺口，而评估时有缺口的序列被排除在分母之外，覆盖率被高估。SILVA 用的是它自己的比对列，缺口只有 0.2%。这个问题还可能出现在 642 和 1138 位。现在设计时要求引物位点处两库的缺口比例都很低（Greengenes &lt; 3%、SILVA &lt; 5%），上面三个方案的位点都满足。</li>
<li><b>两库各自的最优引物不同</b>：同一个位点，Greengenes 上最优的简并引物和 SILVA 上最优的可能不是同一条。所以我只用 Greengenes 设计的引物去评估 SILVA，不把两库各自最优引物的覆盖率当成同一条引物的覆盖率。</li>
</ul>

<h2>八、局限</h2>
<ul>
<li>没有评估 GTDB。本环境的网络策略下载不了 SILVA 官网和 GTDB；SILVA 128 是通过 conda 渠道里的 SEPP 参考包拿到的（2016 年发布的版本，较旧）。</li>
<li>SILVA 的门标签是推断的，继承了 Greengenes 的分类体系，所以门级的评估对分类体系差异不敏感。候选门（OP11、OD1、GN02）序列少，统计不稳。</li>
<li>评估采用 1 个错配、3′ 端 3 个碱基匹配的标准，真实 PCR 的扩增效率可能更低；区域间的效率差异（第一节）在模拟里完全没有体现。</li>
<li>下游 SMURF 需要用 <code>build_region_db_from_fasta</code> 重新生成 3 个区域的 k-mer 数据库，并修改 <code>scott_format_newer_func.m</code> 里写死的 <code>nR = 5</code>。</li>
<li>至少需要用一个已知组成的 mock community 做实验验证，并调整引物配比。</li>
</ul>

<footer>脚本：<code>python_5R/primer_design.py</code>、<code>python_5R/design_primers.py</code>、<code>python_5R/silva_ref.py</code>；数据：<code>docs/primer_design/</code>（候选位点 sites_key.csv / sites_silva.csv / sites_joint.csv、设计结果 designs_single.json、区域偏倚 region_bias.csv）。</footer>
</main></div>
"""
open(OUT, "w", encoding="utf-8").write(head + body)
print("written", OUT, len(head + body))
