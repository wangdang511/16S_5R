"""生成 docs/primer_order_list_v67short.html"""
import re, html
import pandas as pd
HERE = __file__.rsplit("/", 1)[0]; DOCS = HERE.rsplit("/", 1)[0]
old = open(DOCS + "/5R_SMURF_pipeline.html", encoding="utf-8").read()
head = re.sub(r"<title>.*?</title>", "<title>订购清单：V6V7 缩短版</title>", old[:old.index("</style>") + 8], 1)
NUMRE = re.compile(r"^[\d.%/ –→+-]+$")
def cell(v, seq=False): return '<td class="' + ("seq" if seq else ("n" if NUMRE.match(str(v)) else "")) + '">' + html.escape(str(v)) + "</td>"
def tbl(rows, heads, seqcol=None): return '<div class="tbl"><table><thead><tr>' + "".join(f"<th>{c}</th>" for c in heads) + "</tr></thead><tbody>" + "".join("<tr>" + "".join(cell(v, i == seqcol) for i, v in enumerate(r)) + "</tr>" for r in rows) + "</tbody></table></div>"
D = pd.read_csv(f"{HERE}/order_list_5amp_v67short.csv"); I = pd.read_csv(f"{HERE}/v67_info.csv"); P = pd.read_csv(f"{HERE}/v67_primer.csv")
t1 = tbl([[r["name"], r.amplicon, r.seq, r.nt, r.pos, f"{r.Tm}（{r.Tm_lo}–{r.Tm_hi}）", f"{r.GC * 100:.0f}%", int(r.expansions), f"{r.share * 100:.0f}%" if r.n_at_site > 1 else "单条", f"{r.uL_100uM_per_100uL_10x:.2f}"] for _, r in D.iterrows()],
         ["订购名称", "区域", "序列 5′→3′（IUPAC）", "nt", "E. coli 位置", "Tm °C（展开范围）", "GC", "展开数", "位点内比例", "100 µM 母液 µL / 100 µL 10× 混合液"], seqcol=2)
t2 = tbl([[r.region, int(r.bp), f"{r.alone * 100:.1f}%"] for _, r in I.iterrows()], ["区域（E. coli 位置）", "长度 bp", "单独用这一段的属准确率"])
t3 = tbl([[r.site, int(r.oligos), r.prim, int(r.expansions), r.Tm.replace("[", "").replace("]", ""), f"{r.SILVA_mean * 100:.1f}% / {r.SILVA_min * 100:.1f}%", f"{r.GG_mean * 100:.1f}% / {r.GG_min * 100:.1f}%"] for _, r in P.iterrows()],
         ["位点", "寡核苷酸数", "序列 5′→3′", "展开数", "Tm °C", "覆盖率 SILVA：主要门平均 / 最差门", "覆盖率 GG：主要门平均 / 最差门"], seqcol=2)
body = f"""
<div class="wrap"><main style="grid-column:1/-1;width:100%;max-width:980px;margin-inline:auto">
<header>
  <div class="eyebrow">16S_5R · 引物设计</div>
  <h1>订购清单：V6·V7 缩短（229 bp）、V8·V9 保持 288 bp</h1>
  <p class="lede">15 条寡核苷酸、102 个展开序列、10 个引物位点。V6·V7 缩短只去掉保守序列，信息密度反而升高，没有降低；代价是 V6·V7 的正向引物位点没有原来那么保守，需要多 1 条寡核苷酸，对个别菌门（支原体所在的 Tenericutes 等）的覆盖率较弱。这是计算机设计，没有实验验证。</p>
</header>

<h2>一、订购清单</h2>
{t1}
<ul>
<li>相对上一版订购清单，只有 V6·V7 的三条引物不同：A4-F 从 1 条变成 2 条（967–985），A4-R 第 10 位 C→Y；其余引物不变，A5-F 沿用已处理自二聚体的 `GCKACACRCRTGCTACAAT`。</li>
<li>序列、位置、Tm、10× 混合液配法的含义和上一版相同：序列是 IUPAC 简并码，按混合碱基合成；每个位点总量相同（终浓度 250 nM），位点内按单独覆盖率分；没有接头和 UDP。</li>
<li>每条寡核苷酸的 Tm 用 primer3（50 mM Na⁺、2 mM Mg²⁺、0.2 mM dNTP、250 nM 引物）算，括号是各展开序列的范围。池内两两检查：没有严重二聚体或发夹，平均 Tm 54.8–65.7 °C。</li>
<li>扩增子（含引物）：V1·V2 258、V3 186、V4 249、<b>V6·V7 229</b>、V8·V9 288 bp。内联 UDP 每端 10 bp、2×150 读长下，只有 V8·V9 缺 8 bp（1363–1370，不是保守区，平均熵 0.68 bit），其余测通；缺口对属鉴定的影响可忽略，SMURF 拼接读段不需要重叠。</li>
<li>人基因组（hg19）：≤2 个错配 2,940 个位点、1 个潜在产物；≤3 个错配 58 个潜在产物（和上一版 56 相近）；线粒体没有潜在产物。A1-F 仍是人基因组命中最多的引物。</li>
</ul>

<h2>二、V6·V7 缩短有没有降低信息密度？没有</h2>
{t2}
<ul>
<li>被去掉的 928–985 这 58 bp，平均熵只有 0.24 bit，变异位点（熵 >0.5）只占 19%；保留的 986–1176 平均熵 0.67 bit，变异位点占 52%。去掉的是保守序列。</li>
<li>V6（986–1043）和 V7（1117–1173）完整保留：V6 平均熵 1.15 bit，90% 的位点有变异；V7 平均熵 0.81 bit，63% 有变异。</li>
<li>信息密度（平均熵）从 0.57 bit/bp 升到 0.67 bit/bp（+18%），整体设计里每 bp 的边际准确率贡献从 0.0040 升到 0.0051 个百分点。</li>
<li>单独看 V6·V7，属准确率从 78.1% 降到 77.3%（−0.8 个百分点）；放进 5 个扩增子的整体设计，理想准确率只从 95.33% 降到 95.31%（−0.02 个百分点），因为其他区域已经把大部分序列分开了。去掉整个 V6·V7 则降到 94.33%（−1.0 个百分点），说明这个区域有价值，但缩短 61 bp 几乎不影响它。</li>
<li>局限：熵是粗略指标（之前发现它和属准确率相关性只有 0.33），属准确率是代理指标，只测到属。</li>
</ul>

<h2>三、有没有增加引物设计难度？有</h2>
{t3}
<ul>
<li><b>正向引物位点</b>：原来的 906–927（926F 附近）极保守，1 条 22 nt、4 个展开序列，主要门平均 97%，最差门 94%。新位点 967–985（967F 附近）保守度低一些，同样要求覆盖率 ≥93% 时需要 2 条、16 个展开序列，主要门平均降到 94% / 93%，最差门掉到 84%（SILVA）和 63%（Greengenes）。</li>
<li>最弱的是 Tenericutes（Greengenes 63%、SILVA 84%）、SILVA 里 Cyanobacteria 86%、Spirochaetes 90%；Tenericutes 的 A4-R 本来也只有 77–92%，两端都弱，所以这个菌门的 V6·V7 扩增子覆盖率会比较低。其他主要门基本在 92% 以上。</li>
<li>反向引物 A4-R 还是同一个位点（1177–1195），改了一个碱基（C→Y）让覆盖率更高：SILVA 主要门平均 94.3% → 96.1%，GG 92.3% → 94.8%，代价是展开数 4 → 8。</li>
<li>池的规模：寡核苷酸 14 → 15，展开序列 86 → 102（+19%），整体 Tm 范围基本不变；没有新的二聚体问题。</li>
<li>整体扩增子覆盖率：V6·V7 扩增子 SILVA 78% → 77%，GG 91% → 91%；5 个扩增子全部扩出 SILVA 56.2% → 56.0%，GG 68.7% → 68.4%，影响很小。属准确率（三种规则）87.9 / 88.1 / 90.8% → 88.8 / 88.6 / 90.9%，持平或略高。</li>
</ul>
<p>如果你的样本里支原体（Tenericutes）等菌门很重要，V6·V7 位点的覆盖要多留意；也可以保留原来的 906–927 正向引物，但那样扩增子是 290 bp，会有 10 bp 缺口（保守区），见 <code>primer_shorter_amplicons.html</code> 第五节。</p>

<h2>四、限制</h2>
<ul>
<li>覆盖率按“最多 1 个错配、3′ 端 3 个碱基匹配”的序列规则算；数据库是 Greengenes 13_8 和 SILVA 128，不是最新版本；没有实验验证。</li>
<li>各扩增子长度不同（186–288 bp），混合比例没有校正长度引起的扩增偏倚，需要用已知组成的样本实测。</li>
<li>A1-F 在人基因组里有很多近似位点，人 DNA 多的样本会浪费读段。</li>
</ul>
<footer>脚本：<code>python_5R/explore/v67_*.py</code>、<code>order_d1.py</code>；数据：<code>docs/primer_design/order_list_5amp_v67short.csv</code>、<code>v67_*.csv</code>。</footer>
</main></div>
"""
open(DOCS + "/primer_order_list_v67short.html", "w", encoding="utf-8").write(head + body)
print("ok")
