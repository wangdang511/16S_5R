"""生成 docs/primer_order_list_full.html"""
import re, html
import pandas as pd
HERE = __file__.rsplit("/", 1)[0]; DOCS = HERE.rsplit("/", 1)[0]
old = open(DOCS + "/5R_SMURF_pipeline.html", encoding="utf-8").read()
head = re.sub(r"<title>.*?</title>", "<title>完整引物订购清单</title>", old[:old.index("</style>") + 8], 1)
NUMRE = re.compile(r"^[\d.%/ –→+-]+$")
def cell(v, seq=False): return '<td class="' + ("seq" if seq else ("n" if NUMRE.match(str(v)) else "")) + '">' + html.escape(str(v)) + "</td>"
def tbl(rows, heads, seqcol=None): return '<div class="tbl"><table><thead><tr>' + "".join(f"<th>{c}</th>" for c in heads) + "</tr></thead><tbody>" + "".join("<tr>" + "".join(cell(v, i == seqcol) for i, v in enumerate(r)) + "</tr>" for r in rows) + "</tbody></table></div>"
D = pd.read_csv(f"{HERE}/order_list_full.csv"); E = pd.read_csv(f"{HERE}/final_design_eval2.csv"); T = pd.read_csv(f"{HERE}/final_ten_eval2.csv"); S = pd.read_csv(f"{HERE}/final_offtarget_summary.csv").set_index("mm")
core = D[D.group != "备选"]
fm = lambda v: "—" if pd.isna(v) else f"{v:.2f}"
t1 = tbl([[r.order, r.recommend, r.region, r.seq, r.nt, r.pos, f"{r.Tm}（{r.Tm_lo}–{r.Tm_hi}）", int(r.expansions), f"{r.human_per_exp:g}", fm(r.uL_with_sup), fm(r.uL_no_sup)] for _, r in core.iterrows()],
         ["订购名称", "建议", "区域", "序列 5′→3′（IUPAC）", "nt", "E. coli 位置", "Tm °C（展开范围）", "展开数", "人基因组位点 / 展开序列（≤2 错配）", "µL 母液（含补充）", "µL 母液（不含补充）"], seqcol=3)
alt = D[D.group == "备选"]
t2 = tbl([[r.order, r.recommend, r.seq, r.nt, r.pos, f"{r.Tm}（{r.Tm_lo}–{r.Tm_hi}）", int(r.expansions), f"{r.human_per_exp:g}"] for _, r in alt.iterrows()],
         ["订购名称（备选）", "用途", "序列 5′→3′", "nt", "位置", "Tm °C", "展开数", "人基因组位点 / 展开序列"], seqcol=2)
e = E.set_index("design")
cols = ["核心（不含补充）", "推荐：核心 + 4 条支原体补充（延长后）", "核心 + 3 条补充（不含 A1-R）"]
t3 = tbl([[n, int(e.loc[n].oligos), int(e.loc[n].expansions), f"{e.loc[n].Tm_min}–{e.loc[n].Tm_max}", int(e.loc[n].severe), e.loc[n].SILVA_amp, e.loc[n].GG_amp, f"{e.loc[n].SILVA_all * 100:.1f}% / {e.loc[n].GG_all * 100:.1f}%", f"{e.loc[n].ideal * 100:.2f}%", f"{e.loc[n].abs60 * 100:.1f} / {e.loc[n].abs300 * 100:.1f} / {e.loc[n]['frac0.8'] * 100:.1f}%"] for n in cols],
         ["方案", "寡核苷酸", "展开数", "Tm °C", "严重二聚体", "扩增子覆盖 SILVA", "扩增子覆盖 GG", "全部扩出 SILVA / GG", "理想准确率", "属准确率（三种规则）"])
tt = T.set_index("design"); pc = lambda x: f"{x * 100:.0f}%"
t4 = tbl([[a, f"{pc(tt.loc[cols[0]]['GG_' + k])} → {pc(tt.loc[cols[1]]['GG_' + k])}", f"{pc(tt.loc[cols[0]]['SILVA_' + k])} → {pc(tt.loc[cols[1]]['SILVA_' + k])}"] for a, k in (("V1·V2", "A1"), ("V3", "A2"), ("V4", "A3"), ("V6·V7（906 正向）", "A4"), ("V8·V9", "A5"))] + [["至少扩出 3 个扩增子", f"{pc(tt.loc[cols[0]]['GG_>=3扩增子'])} → {pc(tt.loc[cols[1]]['GG_>=3扩增子'])}", f"{pc(tt.loc[cols[0]]['SILVA_>=3扩增子'])} → {pc(tt.loc[cols[1]]['SILVA_>=3扩增子'])}"]],
         ["支原体（Tenericutes）扩增子覆盖：核心 → 含补充", "Greengenes", "SILVA"])
body = f"""
<div class="wrap"><main style="grid-column:1/-1;width:100%;max-width:1040px;margin-inline:auto">
<header>
  <div class="eyebrow">16S_5R · 引物设计</div>
  <h1>完整引物订购清单（推荐版 + 备选）</h1>
  <p class="lede">推荐 18 条寡核苷酸（14 条核心 + 4 条支原体补充），展开序列 114 个，覆盖 V1·V2、V3、V4、V6·V7、V8·V9 五个区域。其中 3 条补充“建议”订，1 条“可选”。V3 推荐长版（338–535 bp 198），V6·V7 推荐保留原来的 906 正向引物（290 bp）。这是计算机设计，没有实验验证。</p>
</header>

<h2>一、推荐订购清单</h2>
{t1}
<ul>
<li>序列用 IUPAC 简并码，按混合碱基合成（R=A/G，Y=C/T，W=A/T，M=A/C，K=G/T，S=G/C，D=A/G/T，B=C/G/T，V=A/C/G）；反向引物的序列已经是 5′→3′ 的引物序列，位置是它在 E. coli 16S（顶链）上占据的范围；只含基因特异序列，没有接头和 UDP，带接头后的 Tm 和二聚体没有评估。</li>
<li><b>“必订”</b>14 条是核心；<b>“建议”</b>3 条是支原体补充（A2-F.2、A4-R.2、A5-F.3），人基因组命中低、Tm 合适，对其他菌门没有害处（覆盖率反而略升）；<b>“可选”</b>1 条（A1-R.4）只在 SILVA 里把 V1·V2 的支原体覆盖提高约 6 个百分点，Greengenes 里没有提高。</li>
<li>混合比例：每个引物位点（共 10 个）总量相同，终浓度 250 nM；位点内有多条时，核心按单独覆盖率分（A1-R 45/35/20，A3-F 90/10，A5-F 50/50），补充各占 10%。表中是“100 µL 的 10× 混合液里要加多少 µL 的 100 µM 母液”，含补充时总量 25 µL（加 75 µL 水），不含补充时同样总量 25 µL（有补充的位点，核心按原比例分满）。每 25 µL PCR 加 2.5 µL。</li>
<li>人基因组（hg19，整池）：≤2 个错配 {int(S.loc[2].loci):,} 个位点、{int(S.loc[2].products)} 个潜在产物；≤3 个错配 {int(S.loc[3].products)} 个潜在产物；线粒体没有潜在产物。人基因组命中最多的是 A1-F（17 nt，291 个/展开序列）；其余都在 30 个以内。</li>
</ul>

<h2>二、备选（选一个版本时才订）</h2>
{t2}
<ul>
<li><b>V3：长版 vs 短版（推荐长版）</b>。长版用 A2-R 516–535（20 nt，2 个展开序列），扩增子 338–535（198 bp），理想准确率 95.38%（短版 95.33%）；对支原体长版好得多：V3 反向位点 Greengenes 92%、SILVA 97%，短版是 73% / 86%（需要再加一条补充才能到 97% / 98%）。代价：长版 Tm 66.9 °C（短版 63.8 °C），是池里最高的。短版更省事、Tm 低，但要加一条补充（`TGCTGGCACATAGTTWGYY`，19 nt，人基因组位点 8.9 个/展开序列）才对支原体有同样覆盖。选短版就不要订 A2-R（长版）。</li>
<li><b>V6·V7：906 正向（推荐）vs 967 缩短版</b>。推荐 906 正向（290 bp）：对支原体几乎全覆盖（96% / 96%），缺 10 bp，在 V6 和 V7 之间的保守间隔区（1046–1055，平均熵 0.15 bit），不丢可变区。缩短版（967–1195，229 bp）不缺口、不需要考虑重叠，但正向位点对支原体弱（74% / 88%），要 2 条寡核苷酸，再加一条补充（`RTACMCGAARAACCTTACC`，14 个/展开序列）。选缩短版就不要订 A4-F（906）。</li>
<li>V8·V9 保持 288 bp：内联 UDP（每端 10 bp）+ 2×150 读长下缺 8 bp（1363–1370），不是保守区（平均熵 0.68 bit），属准确率影响可忽略；需要合并读段的流程要注意。更短的 V8（1223–1407，185 bp）没有缺口，但丢掉 V9，没有放进这份清单。</li>
</ul>

<h2>三、整池评估</h2>
{t3}
{t4}
<ul>
<li>推荐的 18 条加上表中 3 条备选（共 21 条，含各自的展开序列）两两检查，没有严重二聚体和发夹（备注里提到的另外两条备选补充，在延长前的版本做过同样检查，也没有问题）；推荐池平均 Tm 54.8–66.9 °C，退火温度建议先做梯度（大致 52–60 °C）。</li>
<li>补充寡核苷酸对其他菌门没有害处：5 个扩增子全部扩出的序列 SILVA 56.8% → 57.7%，Greengenes 69.8% → 70.2%；属准确率持平或略高（三种规则：88.0 / 88.3 / 91.0% → 88.0 / 88.4 / 91.0%）。</li>
<li>对支原体（Greengenes / SILVA）：补充后至少扩出 3 个扩增子的序列 83% → 90%，78% → 83%。</li>
</ul>

<h2>四、其他备注和限制</h2>
<ul>
<li>各扩增子长度：V1·V2 258、V3 198（长版）、V4 249、V6·V7 290、V8·V9 288 bp。V3 与 V1·V2、V4 的位点之间分别留 ≥20 bp，没有重叠。</li>
<li>混合比例没有校正扩增子长度差异引起的扩增偏倚（V3 最短，V6·V7 最长），需要用已知组成的样本实测后调整各位点的量。</li>
<li>A1-F 对人基因组有很多近似位点，人 DNA 多的样本（比如 FFPE）会浪费读段，也没法延长（数据库里序列在 8 位之前没有数据，评估不了）。</li>
<li>补充寡核苷酸是针对支原体序列设计的，只在 Greengenes（202 条）和 SILVA（355 条，门由 Greengenes 推断）里评估，序列少，百分比误差有几个点；Greengenes 里 Phytoplasma 占了近四分之一，和你的样本可能不一致。</li>
<li>覆盖率按“最多 1 个错配、3′ 端 3 个碱基匹配”的序列规则算，不是实际 PCR 效率；数据库不是最新版本；没有实验验证。</li>
</ul>
<footer>脚本：<code>python_5R/explore/final*.py</code>；数据：<code>docs/primer_design/order_list_full.csv</code>、<code>final_*.csv</code>。</footer>
</main></div>
"""
open(DOCS + "/primer_order_list_full.html", "w", encoding="utf-8").write(head + body)
print("ok")
